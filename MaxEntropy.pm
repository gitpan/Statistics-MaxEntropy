package Statistics::MaxEntropy;

##---------------------------------------------------------------------------##
##  Author:
##      Hugo WL ter Doest       terdoest@cs.utwente.nl
##  Description:
##      Object-oriented implementation of
##      Generalised Iterative Scaling algorithm, 
##	Improved Iterative Scaling algorithm, and
##      Feature Induction algorithm
##      for inducing maximum entropy probability distributions
##  Keywords:
##      Maximum Entropy Modeling
##      Kullback-Leibler Divergence
##      Exponential models
##
##---------------------------------------------------------------------------##
##  Copyright (C) 1998 Hugo WL ter Doest terdoest@cs.utwente.nl
##
##  This library is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This library  is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU Library General Public 
##  License along with this program; if not, write to the Free Software
##  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
##	Globals
##---------------------------------------------------------------------------##
use vars qw($VERSION
	    @ISA
	    @EXPORT
	    $SPARSE
	    $VECTOR_PACKAGE

	    $debug
	    $SAMPLE_size
	    $NEWTON_max_it
	    $KL_max_it
	    $KL_min
	    $NEWTON_min
	    $cntrl_c_pressed
	    $cntrl_backslash_pressed
	    );


##---------------------------------------------------------------------------##
##	Require libraries
##---------------------------------------------------------------------------##
use strict;
use diagnostics -verbose;
$SPARSE = 1;
if ($SPARSE) {
    $VECTOR_PACKAGE = "Statistics::SparseVector";
    use Statistics::SparseVector;
}
else {
    $VECTOR_PACKAGE = "Bit::Vector";
    use Bit::Vector;
}
use POSIX;
use Carp;
use Data::Dumper;
require Exporter;
require AutoLoader;

@ISA = qw(Exporter AutoLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw($KL_min
	     $NEWTON_min
	     $debug
	     $nr_add
	     $KL_max_it
	     $NEWTON_max_it
	     $SAMPLE_size

	     new
	     DESTROY
	     write
	     scale
	     dump
	     undump
	     
	     fi
	     random_parameters
	     set_parameters_to
	     write_parameters
	     write_parameters_with_names
	     );

$VERSION = '0.8';


# default values for some configurable parameters
$NEWTON_max_it = 200;
$NEWTON_min = 0.001;
$KL_max_it = 100;
$KL_min = 0.001;
$debug = 0;
$SAMPLE_size = 250; # the size of MC samples
$cntrl_c_pressed = 0;
$cntrl_backslash_pressed = 0;
$SIG{INT} = \&catch_cntrl_c;
$SIG{QUIT} = \&catch_cntrl_backslash;


# interrrupt routine for control c
sub catch_cntrl_c {
    my($signame) = shift;

    $cntrl_c_pressed = 1;
    die "<CONTROL-C> pressed\n";
}


# interrrupt routine for control \ (originally core-dump request)
sub catch_cntrl_backslash {
    my($signame) = shift;

    $cntrl_backslash_pressed = 1;
}


# creates a new event space
# depending on the $arg parameter samples it or reads it from a file
sub new {
    my($this, $arg) = @_;

    # for calling $self->new($someth):
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;
    $self->{SCALER} = "gis"; # default
    $self->{SAMPLING} = "corpus"; # default
    $self->{NR_CLASSES} = 0;
    $self->{NR_EVENTS} = 0;
    $self->{NR_FEATURES} = 0;
    if ($arg) { # hey a filename
	$self->read($arg);
    }
    return($self);
}


# decides how to sample, "enum", "mc", or "corpus"
sub sample {
    my($self) = @_;

    my($sample);

    if ($self->{SAMPLING} eq "mc") {
	$sample = $self->new();
	$sample->{SCALER} = $self->{SCALER};
	$sample->{NR_FEATURES} = $self->{NR_FEATURES};
	# refer to the parameters of $self
	$sample->{PARAMETERS} = $self->{PARAMETERS};
	$sample->{CORRECTION_PARAMETER} = $self->{CORRECTION_PARAMETER};
	$sample->{E_REF} = $self->{E_REF};
	$sample->{THIS_IS_A_SAMPLE} = 1;
	$sample->mc($self);
	$self->prepare_model();
    }
    elsif ($self->{SAMPLING} eq "enum") {
	$sample = $self->new();
	$sample->{SCALER} = $self->{SCALER};
	$sample->{NR_FEATURES} = $self->{NR_FEATURES};
	$sample->{PARAMETERS} = $self->{PARAMETERS};
	$sample->{CORRECTION_PARAMETER} = $self->{CORRECTION_PARAMETER};
	$sample->{E_REF} = $self->{E_REF};
	$sample->{THIS_IS_A_SAMPLE} = 1;
	$sample->enum();
	$sample->prepare_model();
    }
    else { # taken to be "corpus"
	$sample = $self;
    }
    return($sample);
}


# makes sure that when prepare_model is called, everything is recomputed
sub clear {
    my($self) = @_;
    
    undef $self->{PARAMETERS_INITIALISED};
    $self->{PARAMETERS_CHANGED} = 1;
    $self->{CLASSES_CHANGED} = 1;
}



sub DESTROY {
    my($self) = @_;
    
    if ($cntrl_c_pressed) {
	$self->dump();
    }
}


# reads an events file, dies in case of inconsistent lines
# syntax first line: <name> <tab> <name> <tab> ..... <newline>
# syntax other lines: <freq> <bitvector> <newline>
sub read {
    my($self, $file) = @_;

    my($features,
       $feature_names);

    $feature_names = "";
    open(EVENTS, $file) ||
	$self->die("Could not open $file\n");
    print "Opened $file\n";

    # read the names of the features, skip comment
    do {
	$feature_names = <EVENTS>;
    } until ($feature_names !~ /\#.*/);
    chomp $feature_names;
    $self->{FEATURE_NAMES} = [split(/\t/, $feature_names)];

    # read the bitvectors
    while (<EVENTS>) {
	if (!/\#.*/) {
	    chomp;

	    ($self->{FREQ}[$self->{NR_CLASSES}], $features) = split;
	    if ($self->{FREQ} == 0) {
		$self->die("Class $self->{NR_CLASSES} has zero probability\n");
	    }
	    $self->{NR_EVENTS} += $self->{FREQ}[$self->{NR_CLASSES}];

	    # if first event set nr_features
	    if ($self->{NR_CLASSES} == 0) {
		$self->{NR_FEATURES} = length($features);
	    }
	    # else check nr of features for this event
	    else {
		if (length($features) != $self->{NR_FEATURES}) {
		    $self->die("Events file corrupt (class $self->{NR_CLASSES})\n");
		}
	    }
	    # create and initialise bit vector
	    $self->{CLASSES}[$self->{NR_CLASSES}] = 
	      $VECTOR_PACKAGE->new_Bin($self->{NR_FEATURES}, $features);
	    $self->{NR_CLASSES}++;
	}
    }
    close(EVENTS);

    print "Read $self->{NR_EVENTS} events, $self->{NR_CLASSES} classes, " . 
	"and $self->{NR_FEATURES} features\n";
    print "Closed $file\n";

    $self->{FILENAME} = $file;
    $self->{CLASSES_CHANGED} = 1;
    $self->{PARAMETERS_CHANGED} = 1;
}


# reads an initial distribution
# syntax: one parameter per line
sub read_parameters {
    my($self, $file) = @_;

    my($i);

    $i = 0;
    open(DISTR,$file) ||
	$self->die("Could not open $file\n");
    print "Opened $file\n";

    while (<DISTR>) {
	if (!/\#.*/) {
	    chomp;
	    $self->{PARAMETERS}[$i++] = $_;
	}
    }

    close(DISTR);
    if ($i != $self->{NR_FEATURES}) {
	$self->die("Initial distribution file corrupt\n");
    }
    print "Read $i parameters; closed $file\n";
    $self->{PARAMETERS_CHANGED} = 1;
}


# writes the the current parameters
# syntax: <parameter> <newline>
sub write_parameters {
    my($self, $file) = @_;

    my($i);

    open(DISTR,">$file") ||
	$self->die("Could not open $file\n");
    print "Opened $file\n";

    for ($i = 0; $i < $self->{NR_FEATURES}; $i++) {
	if ($self->{FEATURE_IGNORED}{$i}) {
	    print DISTR "IGNORED\n";
	}
	else {
	    print DISTR "$self->{PARAMETERS}[$i]\n";
	}
    }

    close(DISTR);
    print "Closed $file\n";
}


# writes the the current features with their parameters
# syntax first line: <$nr_features> <newline>
# syntax last line: <bitmask> <newline>
# syntax other lines: <name> <parameter> <newline>
sub write_parameters_with_names {
    my($self, $file) = @_;

    my($x,
       $bitmask);

    open(DISTR,">$file") ||
	$self->die("Could not open $file\n");
    print "Opened $file\n";

    print DISTR "$self->{NR_FEATURES}\n";
    $bitmask = "";
    for ($x = 0; $x < $self->{NR_FEATURES}; $x++) {
	print DISTR "$self->{FEATURE_NAMES}[$self->{NR_FEATURES} - $x - 1]\t" .
	    "$self->{PARAMETERS}[$x]\n";
	if ($self->{FEATURE_IGNORED}{$x}) {
	    $bitmask .= "0";
	}
	else {
	    $bitmask .= "1";
	}
    }
    print DISTR "$bitmask\n";

    close(DISTR);
    print "Closed $file\n";
}


# generate random parameters
sub random_parameters {
    my($self) = @_;

    my($f);

    srand();
    for ($f = 0; $f < $self->{NR_FEATURES}; $f++) {
	$self->{PARAMETERS}[$f] = rand() + 1;
    }
    if ($self->{SCALER} eq "gis") {
	$self->{CORRECTION_PARAMETER} = rand();
    }
    $self->{PARAMETERS_CHANGED} = 1;
}


# sets parameters to $val
sub set_parameters_to {
    my($self, $val) = @_;

    my($f);

    for ($f = 0; $f < $self->{NR_FEATURES}; $f++) {
	$self->{PARAMETERS}[$f] = $val;
    }
    if ($self->{SCALER} eq "gis") {
	$self->{CORRECTION_PARAMETER} = $val;
    }
    $self->{PARAMETERS_CHANGED} = 1;
}


# initialise if !$self->{PARAMETERS_INITIALISED}; subsequent calls 
# of scale (by fi) should not re-initialise parameters
sub init_parameters {
    my($self) = @_;

    if (!$self->{PARAMETERS_INITIALISED}) {
	if ($self->{SCALER} eq "gis") {
	    $self->set_parameters_to(1);
	}
	else {
	    $self->random_parameters();
#	    $self->set_parameters_to(1);
	}
	$self->{PARAMETERS_INITIALISED} = 1;
    }
}


# make sure \tilde{p} << q_0
# constant feature functions are forbidden: that is why
# we check whether for all features \sum_x f(x) > 0
# and \sum_x f(x) != $corpus_size
sub check {
    my($self) = @_;

    my ($x,
	$f,
	$sum);

    $sum = 0;
    for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
	if ($self->{CLASS_PROBS}[$x] == 0) {
	    print "Initial distribution not ok; event $x\n";
	}
	$sum += $self->{CLASS_PROBS}[$x];
    }
    if ($debug && ($sum != 1)) {
	print "Sum of the probabilities: $sum\n";
    }

    for ($f = 0; $f < $self->{NR_FEATURES}; $f++) {
	$sum = 0;
	for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
	    $sum += $self->{CLASSES}[$x]->bit_test($f);
	}
	if (!$sum || ($sum == $self->{NR_CLASSES})) {
	    print "Feature ", $f + 1, " is constant ($sum), and will be ignored\n";
	    $self->{FEATURE_IGNORE}{$f} = 1;
	}
    }
}


# writes events to a file 
# usefull in case new features have been added
# syntax: same as input events file
sub write {
    my($self, $file) = @_;

    my($x, $f);

    # prologue
    open(EVENTS,">$file") ||
	$self->die("Could not open $file\n");
    print "Opened $file\n";

    # write a line with the feature names
    print EVENTS join("\t", $self->{FEATURE_NAMES}),"\n";
    # write the events themselves
    for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
	print EVENTS $self->{FREQ}[$x],"\t";
	print EVENTS $self->{CLASSES}[$x]->to_Bin(), "\n";
    }

    # close the file and tell you did that
    close EVENTS;
    print "Wrote $self->{NR_EVENTS} events, $self->{NR_CLASSES} classes, " . 
	"and $self->{NR_FEATURES} features\n";
    print "Closed $file\n";
}


# makes enum strings from Bit::Vector internal (C) format if necessary
sub perldata {
    my($self) = @_;

    @{$self->{CLASSES}} = map {$_->to_Enum()} @{$self->{CLASSES}};
}


# makes Bit::Vector internal format from bitstrings if necessary
sub cdata {
    my($self) = @_;

    @{$self->{CLASSES}} = map 
    {$VECTOR_PACKAGE->new_Enum($self->{NR_FEATURES}, $_)} 
    @{$self->{CLASSES}};
}


# reads a dump, and evaluates it into an object
sub undump {
    my($class, $file) = @_;

    my($x,
       $VAR1);

    # open, slurp, and close file
    open(UNDUMP, "$file") ||
	croak "Could not open $file\n";
    print "Opened $file\n";
    undef $/;
    $x = <UNDUMP>;
    $/ = "\n";
    close(UNDUMP);

    # and undump
    eval $x;
    if (!$SPARSE) {
	$VAR1->cdata();
    }
    print "Undumped $VAR1->{NR_EVENTS} events, $VAR1->{NR_CLASSES} classes, " . 
	"and $VAR1->{NR_FEATURES} features\n";
    print "Closed $file\n";
    return($VAR1);
}


# makes dump of the event space using Data::Dumper
sub dump {
    my($self, $file) = @_;

    my(@bitvecs,
       $dump,
       %features,
       $f);

    if (!$file) {
	$file = POSIX::tmpnam();
    }
    open(DUMP, ">$file") ||
	croak "Could not open $file\n";
    print "Opened $file\n";

    # build something that we can sort
    # ONLY FOR CORPUS!
    if (!$self->{THIS_IS_A_SAMPLE}) {
    for ($f = 0; $f < $self->{NR_FEATURES}; $f++) {
        $features{$self->{FEATURE_NAMES}[$self->{NR_FEATURES} - $f - 1]} = 
	    $self->{PARAMETERS}[$f];
    }
    if ($self->{NEED_CORRECTION_FEATURE} && ($self->{SCALER} eq "gis")) {
        $features{"correction$self->{M}"} = 
	    $self->{CORRECTION_PARAMETER};
    }
    # and print it into $self
    $self->{FEATURE_SORTED} = join(' > ',
				   sort {
				       if ($features{$b} == $features{$a}) {
					   return($b cmp $a)} 
				       else {
					   return ($features{$b} <=> $features{$a})
					   }
				   }
				   keys(%features));
}

    # save classes
    if (!$SPARSE) {
	@bitvecs = @{$self->{CLASSES}};
    }
    $dump = Data::Dumper->new([$self]);
    # perldata makes bitstrings
    if (!$SPARSE) {
	$dump->Freezer('perldata');
    }
    print DUMP $dump->Dump();
    # restore classes
    if (!$SPARSE) {
	@{$self->{CLASSES}} = @bitvecs;
    }

    print "Dumped $self->{NR_EVENTS} events, $self->{NR_CLASSES} classes, " . 
	"and $self->{NR_FEATURES} features\n";

    close(DUMP);
    print "Closed $file\n";
}


# $msg is logged, the time is logged, a dump is created, and the
# program dies with $msg
sub die {
    my($self, $msg) = @_;

    $self->log($msg);
    $self->log(time());
    $self->dump();
    croak $msg;
}


# prints a msg to STDOUT, and appends it to $self->{LOG}
# so an emergency dump will contain some history information
sub log {
    my($self, $x) = @_;

    $self->{LOG} .= $x;
    print $x;
}


# computes f_# for alle events; results in @sample_nr_feats_on
# computes %$sample_m_feats_on; a HOL from m 
sub active_features {
    my($self) = @_;

    my($i,
       $j);

    if ($self->{CLASSES_CHANGED}) {
	# M is needed for both gis and iis
	# NEED_CORRECTION_FEATURE is for gis only
	# NR_FEATURES_ACTIVE for iis only
	$self->{M} = 0;
	$self->{NEED_CORRECTION_FEATURE} = 0;
	undef $self->{NR_FEATURES_ACTIVE};
	for ($i = 0; $i < $self->{NR_CLASSES}; $i++) {
	    $self->{NR_FEATURES_ACTIVE}[$i] = 0;
	    for ($j = 0; $j < $self->{NR_FEATURES}; $j++) {
		$self->{NR_FEATURES_ACTIVE}[$i] += 
		    $self->{CLASSES}[$i]->bit_test($j);
	    }
	    if (!$self->{M}) { 
		# is undefined!
		$self->{M} = $self->{NR_FEATURES_ACTIVE}[$i];
	    }
	    elsif ($self->{NR_FEATURES_ACTIVE}[$i] > $self->{M}) {
		# higher nr_features_active found
		$self->{M} = $self->{NR_FEATURES_ACTIVE}[$i];
		$self->{NEED_CORRECTION_FEATURE} = 1;
	    }
	    if ($debug) {
		print "f_#($i) = $self->{NR_FEATURES_ACTIVE}[$i]\n";
	    }
	}
	if ($debug) {
	    print "M = $self->{M}\n";
	}
	# set up a hash from m to classes HOL; and the correction_feature
	# M_FEATURES_ACTIVE IS FOR iis
	# CORRECTION_FEATURE FOR gis
	undef $self->{M_FEATURES_ACTIVE};
	for ($i = 0; $i < $self->{NR_CLASSES}; $i++) {
	    push @{$self->{M_FEATURES_ACTIVE}{$self->{NR_FEATURES_ACTIVE}[$i]}}, 
	         $i;
	    $self->{CORRECTION_FEATURE}[$i] = $self->{M} - 
		$self->{NR_FEATURES_ACTIVE}[$i];
	    # changed $self->{M} into $self->{NR_FEATURES}
#	    $self->{CORRECTION_FEATURE}[$i] = $self->{NR_FEATURES} - 
#		$self->{NR_FEATURES_ACTIVE}[$i];
	}
	if ($debug) {
	    print "M = $self->{M}\n";
	}
	# observed feature expectations
	if (!$self->{THIS_IS_A_SAMPLE}) {
	    for ($j = 0; $j < $self->{NR_FEATURES}; $j++) {
		# observed feature expectations; gis and iis
		$self->E_reference($j);
	    }
	    $self->E_reference_correction();
	}
	undef $self->{CLASSES_CHANGED};
    }
}


# compute the class probabilities according to the parameters
sub prepare_model {
    my($self) = @_;

    my ($x, 
	$f, 
	$sum);

    $self->active_features();
    if ($self->{PARAMETERS_CHANGED}) {
	$self->{Z} = 0;
	for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {    
	    $sum = 0;
	    for ($f = 0; $f < $self->{NR_FEATURES}; $f++) {
		if (!$self->{FEATURE_IGNORE}{$f} && 
		    $self->{CLASSES}[$x]->bit_test($f)) {
		    $sum += $self->{PARAMETERS}[$f];
		}
	    }
	    if ($self->{NEED_CORRECTION_FEATURE} && ($self->{SCALER} eq "gis")) {
		$sum += $self->{CORRECTION_FEATURE}[$x] * 
		    $self->{CORRECTION_PARAMETER};
	    }
	    $self->{CLASS_WEIGHTS}[$x] = exp($sum);
	    $self->{Z} += $self->{CLASS_WEIGHTS}[$x];
	}
	# normalise
	for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {    
	    $self->{CLASS_PROBS}[$x] = $self->{CLASS_WEIGHTS}[$x] / $self->{Z};
	}
	# expectations
	for ($f = 0; $f < $self->{NR_FEATURES}; $f++) {
	    $self->E_loglinear($f);
	}
	$self->E_loglinear_correction();
	# A_{mj}
	if ($self->{SCALER} eq "iis") {
	    $self->A();
	}
	if (!$self->{THIS_IS_A_SAMPLE}) {
	    $self->entropies();
	}
	$self->check();
	undef $self->{PARAMETERS_CHANGED};
    }
}


sub E_loglinear {
    my($self, $i) = @_;

    my($x,
       $sum);

    $sum = 0;
    for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {    
	if ($self->{CLASSES}[$x]->bit_test($i)) {
	    $sum += $self->{CLASS_PROBS}[$x];
	}
    }
    $self->{E_LOGLIN}[$i] = $sum;
}


sub E_loglinear_correction {
    my($self) = @_;

    my($x,
       $sum);

    if ($self->{NEED_CORRECTION_FEATURE} && ($self->{SCALER} eq "gis")) {
	$sum = 0;
	for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {    
	    $sum += $self->{CORRECTION_FEATURE}[$x] * $self->{CLASS_PROBS}[$x];
	}
	$self->{E_LOGLIN}[$self->{NR_FEATURES}] = $sum;
    }
}


sub E_reference {
    my($self, $i) = @_;

    my($x,
       $sum);

    $sum = 0;
    for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
	if ($self->{CLASSES}[$x]->bit_test($i)) {
	    $sum += $self->{FREQ}[$x];
	}
    }
    $self->{E_REF}[$i] = $sum / $self->{NR_EVENTS};
}


sub E_reference_correction {
    my($self) = @_;

    my($x,
       $sum);

    if (($self->{SCALER} eq "gis") && ($self->{NEED_CORRECTION_FEATURE})) {
	$sum = 0;
	for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
	    $sum += $self->{CORRECTION_FEATURE}[$x] * 
		$self->{FREQ}[$x];
	}
	$self->{E_REF}[$self->{NR_FEATURES}] = $sum / $self->{NR_EVENTS};
    }
}


# compute several entropies
sub entropies {
    my($self) = @_;

    my ($i, 
	$p,
	$log_p,
	$p_ref,
	$log_p_ref);

    $self->{H_p} = 0;
    $self->{H_cross} = 0;
    $self->{H_p_ref} = 0;
    $self->{KL} = 0;
    for ($i = 0; $i < $self->{NR_CLASSES}; $i++) {
	$p = $self->{CLASS_PROBS}[$i];
	# we don't know whether $p > 0
	$log_p = ($p > 0) ? log($p) : 0;
	$p_ref = $self->{FREQ}[$i] / $self->{NR_EVENTS};
	# we know that $p_ref > 0
	$log_p_ref = log($p_ref);
	$self->{H_p} -= $p * $log_p;
	$self->{H_cross} -= $p_ref * $log_p;
	$self->{KL} += $p_ref * ($log_p_ref - $log_p);
	$self->{H_p_ref} -= $p_ref * $log_p_ref;
	if ($p == 0) {
	    $self->log("entropies: skipping event $i (p^n($i) = 0)\n");
	}
    }
    $self->{L} = -$self->{H_cross};
}


# for GIS, if f_# is not a constant function
sub correction_feature {
    my($self, $bitvec) = @_;

    my($i,
       $sum);

    $sum = 0;
    for ($i = 0; $i < $self->{NR_FEATURES}; $i++) {
	$sum += $bitvec->bit_test($i);
    }
    return($self->{M} - $sum);
}


# unnormalised density
sub weight {
    my($self, $bitvec) = @_;

    my ($f, 
	$sum);

    $sum = 0;
    for ($f = 0; $f < $self->{NR_FEATURES}; $f++) {
	if (!$self->{FEATURE_IGNORE}{$f} &&
	    $bitvec->bit_test($f)) {
	    $sum += $self->{PARAMETERS}[$f];
	}
    }
    if ($self->{NEED_CORRECTION_FEATURE}) {
	$sum += $self->correction_feature($bitvec) *
	    $self->{CORRECTION_PARAMETER};
    }
    return(exp($sum));
}


# computes the probability of a bitvector
sub prob {
    my($self, $bitvec) = @_;

    return($self->weight($bitvec) / $self->{Z});
}


# computes the `a' coefficients of 
# \sum_{m=0}^{M} a_{m,j}^{(n)} e^{\alpha^{(n)}_j m}
# according to the current distribution
sub A {
    my($self) = @_;

    my($j, 
       $m,
       $class);

    for  ($j = 0; $j < $self->{NR_FEATURES}; $j++) {
	for ($m = 0; $m <= $self->{M}; $m++) {
	    if ($m == 0) {
		$self->{A}[0][$j] = -$self->{E_REF}[$j];
	    }
	    else {
		$self->{A}[$m][$j] = 0;
		for $class (@{$self->{M_FEATURES_ACTIVE}{$m}}) {
		    if ($self->{CLASSES}[$class]->bit_test($j)) {
			$self->{A}[$m][$j] += $self->{CLASS_PROBS}[$class];
		    }
		}
		if ($debug) {
		    print "a[$m][$j] = $self->{A}[$m][$j]\n";
		}
	    }
	}
    }
}


# enumerates complete Omega; do not use if $nr_features is large!
sub enum {
    my($self) = @_;

    my($f,
       $vec);

    $vec = $VECTOR_PACKAGE->new($self->{NR_FEATURES});
    for ($f = 0; $f < 2 ** $self->{NR_FEATURES}; $f++) {
	push @{$self->{CLASSES}}, $vec;
	push @{$self->{FREQ}}, 1;
	$vec++;
    }
    $self->{NR_CLASSES} = 2 ** $self->{NR_FEATURES};
    $self->{NR_EVENTS} = $f;
    $self->{PARAMETERS_CHANGED} = 1;
    $self->{CLASSES_CHANGED} = 1;
}


#
# Monte Carlo sampling with the Metropolis update
#

# returns heads up with probability $load 
sub loaded_die {
    my($load) = @_;

    (rand() <= $load) ? 1 : 0;
}


# samples from the probability distribution of $other to create $self
# we use the so-called Metropolis update R = h(new)/h(old)
# Metropolis algorithm \cite{neal:probabilistic}
sub mc {
    my($self, $other, $type) = @_;

    my($R,
       $weight,
       $state,
       $old_weight,
       $k,
       %events
       );

    srand();
    # take some class from the sample space as initial state
    $state = $VECTOR_PACKAGE->new($self->{NR_FEATURES});
    # make sure there are no constant features!
    $state->Fill();
    $events{$state->to_Bin()}++;
    $state->Empty();
    $weight = 1;
    # iterate 
    $k = 0;
    do {
	$old_weight = $weight;
	if ($state->bit_flip($k)) {
	    $weight += $self->{PARAMETERS}[$k];
	}
	else {
	    $weight -= $self->{PARAMETERS}[$k];
	}
	$R = exp($weight - $old_weight);
	if (!loaded_die(1 < $R ? 1 : $R)) { # stay at the old state
	    $state->bit_flip($k);
	    $weight = $old_weight;
	}
	$events{$state->to_Bin()}++;
	# next component
	$k = ($k + 1) % $self->{NR_FEATURES};
    } until (scalar(keys(%events)) == $SAMPLE_size);
    for (keys(%events)) {
	push @{$self->{CLASSES}}, $VECTOR_PACKAGE->new_Bin($self->{NR_FEATURES}, $_);
    }
    $self->{NR_CLASSES} = scalar(keys(%events)) - 1;
	
    $self->{CLASSES_CHANGED} = 1;
    $self->{PARAMETERS_CHANGED} = 1;
}


#
# IIS
#
sub a_func {
    my($self, $j, $x, $e_x) = @_;

    my($m,
       $sum_func,
       $sum_deriv,
       $a_x_m);

    $sum_func = $self->{"A"}[0][$j];
    $sum_deriv = 0;
    for ($m = 1; $m <= $self->{M}; $m++) {
	if ($self->{"A"}[$m][$j] != 0) {
	    if ($e_x) {
		$a_x_m = $self->{"A"}[$m][$j] * ($x ** $m);
	    }
	    else {
		$a_x_m = $self->{"A"}[$m][$j] * exp($x * $m);
	    }
	    $sum_func += $a_x_m;
	    $sum_deriv += $m * $a_x_m;
	}
    }
    return($sum_func, $sum_deriv);
}


# solves \alpha from 
# \sum_{m=0}^{M} a_{m,j}^{(n)} e^{\alpha^{(n)}_j m}=0
sub iis_estimate_with_newton {
    my($self, $i, $e_x) = @_;

    my($x, 
       $old_x,
       $deriv_res,
       $func_res,
       $k);

    $x = 1;
    $k = 0;

    # do newton's method
    do {
	# save old x
	$old_x = $x;
	# compute new x
	($func_res, $deriv_res) = $self->a_func($i, $x, $e_x);
	if (($deriv_res eq 'NaN') || ($deriv_res == 0) ||
	    ($deriv_res eq 'Infinity')) {
	    print "a_deriv($i, $x) = $deriv_res\n";
	}
	if (($func_res eq 'NaN') || ($func_res eq 'Infinity') ||
	    ($func_res == 0)) {
	    print "a_func($i, $x) = $func_res\n";
	}
	$x -= ($func_res / $deriv_res);
	if ($x eq 'NaN') {
	    print "$func_res / $deriv_res = NaN\n";
	    print "feature $i will be ignored in future iterations\n";
	    $self->{FEATURE_IGNORE}{$i} = 1;
	    return(0);
	}
    } until ((abs($x - $old_x) <= $NEWTON_min) ||
	     ($k++ > $NEWTON_max_it));
    if ($debug) {
	print "Estimated gamma_$i with Newton's method: $x\n";
    }
    return($e_x ? log($x) : $x);
}


# the iterative scaling algorithms
sub scale {
    my($self, $sampling, $scaler) = @_;

    my($k,
       $i,
       $kl,
       $old_kl,
       $diff,
       $sample,
       $old_correction_parameter,
       @old_parameters);

    if ($sampling) {
	$self->{SAMPLING} = $sampling;
    }
    if ($scaler) {
	$self->{SCALER} = $scaler;
    }
    $self->init_parameters();
    $self->prepare_model();
    if ($self->{SAMPLING} eq "enum") {
	$sample = $self->sample();
    }
    $self->log("($self->{SCALER}, $self->{SAMPLING}): H(p_ref)=$self->{H_p_ref}\nit.\tD(p_ref||p)\t\tH(p)\t\t\tL(p_ref,p)\t\ttime\n0\t$self->{KL}\t$self->{H_p}\t$self->{L}\t" . time() . "\n");
    $k = 0;
    $kl = 1e99;
    do {
	# store parameters for reverting if converging stops
	@old_parameters = @{$self->{PARAMETERS}};
	$old_correction_parameter = $self->{CORRECTION_PARAMETER};
	if ($self->{SAMPLING} ne "enum") { 
	    if ($sample) {
		$sample->DESTROY();
	    }
	    $sample = $self->sample();
	}
	$sample->prepare_model();
	for ($i = 0; $i < $self->{NR_FEATURES}; $i++) {
	    if (!$self->{FEATURE_IGNORE}{$i} && !$sample->{FEATURE_IGNORE}{$i}) {
		if ($self->{SCALER} eq "gis") {
		    $self->{PARAMETERS}[$i] *= ($self->{E_REF}[$i] / 
			$sample->{E_LOGLIN}[$i]) ** (1 / $self->{M});
		}
		else {
		    $self->{PARAMETERS}[$i] += 
			$sample->iis_estimate_with_newton($i);
		}
	    }
	}
	if (($self->{SCALER} eq "gis") && ($self->{NEED_CORRECTION_FEATURE})) {
	    $self->{CORRECTION_PARAMETER} *=
		($self->{E_REF}[$self->{NR_FEATURES}] / 
		 $sample->{E_LOGLIN}[$self->{NR_FEATURES}]) ** 
		     (1 / $self->{M});
	}
	$self->{PARAMETERS_CHANGED} = 1;
	$self->prepare_model();
	$diff = $kl - $self->{KL};
	$kl = $self->{KL};

	$k++;
	$self->log("$k\t$self->{KL}\t$self->{H_p}\t$self->{L}\t" . time() . "\n");
	if ($debug) {
	    $self->check();
	}
	if ($diff < 0) {
	    $self->log("Scaling is not converging (anymore); will revert parameters!\n");
	    # restore old parameters
	    $self->{PARAMETERS} = \@old_parameters;
	    $self->{CORRECTION_PARAMETER} = $old_correction_parameter;
	    $self->{PARAMETERS_CHANGED} = 1;
	    $self->prepare_model();
	}
	if ($cntrl_backslash_pressed) { 
	    $self->dump();
	    $cntrl_backslash_pressed = 0;
	}
    } until ($diff <= $KL_min || 
	     ($k > $KL_max_it) ||
	     ($diff < 0));
}


#
# Field Induction Algorithm
#

# add feature $g to $self
sub add_feature {
    my($self, $candidates, $g) = @_;

    my($i);

    $self->{NR_FEATURES}++;
    for ($i = 0; $i < $self->{NR_CLASSES}; $i++) {
	$self->{CLASSES}[$i]->Interval_Substitute($candidates->{CANDIDATES}[$i],
						  $self->{CLASSES}[$i]->Size(),
						  0, $g, 1);
    }
    if ($self->{SCALER} eq "gis") {
	$self->{PARAMETERS}[$self->{NR_FEATURES} - 1] = 1;
    }
    else {
	$self->{PARAMETERS}[$self->{NR_FEATURES} - 1] = $candidates->{ALPHA}[$g];
    }
    unshift @{$self->{FEATURE_NAMES}}, $candidates->{CANDIDATE_NAMES}[$g];
    $self->{PARAMETERS_CHANGED} = 1;
    $self->{CLASSES_CHANGED} = 1;
    $self->prepare_model();
}


# remove feature $g
sub remove_feature {
    my($self, $g) = @_;

    my($i
       );

    for ($i = 0; $i < $self->{NR_CLASSES}; $i++) {
	# substitute offset $g length 1 by nothing
	$self->{CLASSES}[$i]->Interval_Substitute($self->{CLASSES}[$i],
						  $g, 1, 0, 0);
    }
    splice(@{$self->{PARAMETERS}}, $g, 1);
    splice(@{$self->{FEATURE_NAMES}}, $self->{NR_FEATURES} - 1 - $g, 1);
    $self->{NR_FEATURES}--;
    $self->{PARAMETERS_CHANGED} = 1;
    $self->{CLASSES_CHANGED} = 1;
    $self->prepare_model();
}


# checks for $event, if not there adds it, otherwise increases its {FREQ}
sub add_event {
    my($self, $event) = @_;

    my($i,
       $found);

    $found = 0;
    for ($i = 0; $i < $self->{NR_CLASSES}; $i++) {
	$found = ($event->Compare($self->{CLASSES}[$i]) == 0);
	if ($found) {
	    $self->{FREQ}[$i]++;
	    last;
	}
    }
    if (!$found) {
	$self->{CLASSES}[$self->{NR_CLASSES}] = $event->Clone();
	$self->{FREQ}[$self->{NR_CLASSES}] = 1;
	$self->{NR_CLASSES}++;
    }
    $self->{NR_EVENTS}++;
}


# computes the gain for all $candidates
sub gain {
    my($self, $candidates) = @_;

    my($c,
       $x,
       $kl,
       $below,
       $above,
       $sum_p_ref,
       $sum_p);

    $candidates->{MAX_GAIN} = 0;
    $candidates->{BEST_CAND} = 0;
    for ($c = 0; $c < $candidates->{NR_CANDIDATES}; $c++) {
	if (!$candidates->{ADDED}{$c}) {
	    $sum_p_ref = 0;
	    $sum_p = 0;
	    for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
		if ($candidates->{CANDIDATES}[$x]->bit_test($c)) {
		    $sum_p += $self->{CLASS_PROBS}[$x];
		    $sum_p_ref += $self->{FREQ}[$x];
		}
	    }
	    $sum_p_ref /= $self->{NR_EVENTS};
	    $above = $sum_p_ref * (1 - $sum_p);
	    $below = $sum_p * (1 - $sum_p_ref);
	    if ((($above > 0) && ($below > 0)) || (($above < 0) && ($below < 0))) {
		$candidates->{ALPHA}[$c] = log(abs($above)) - log(abs($below));
	    }
	    else {
		$self->die("Cannot take log of negative/zero value: $above / $below\n");
	    }
	    
	    # temporarily add feature to classes and compute $gain
	    $kl = $self->{KL};
	    $self->add_feature($candidates, $c);
	    $candidates->{GAIN}[$c] = $kl - $self->{KL};
	    $self->log("G($c, $candidates->{ALPHA}[$c]) = $candidates->{GAIN}[$c]\n");
	    if (($candidates->{MAX_GAIN} <= $candidates->{GAIN}[$c])) {
		$candidates->{MAX_GAIN} = $candidates->{GAIN}[$c];
		$candidates->{BEST_CAND} = $c;
	    }
	    # remove the feature
	    $self->remove_feature($self->{NR_FEATURES} - 1);
	}
    }
}


# adds the $n best candidates
sub fi {
    my($self, $scaler, $candidates, $n, $sample) = @_;

    my ($i,
	$kl);

    $self->log("(fi, $scaler, $sample, $n)\n");
    if ($scaler) {
	$self->{SCALER} = $scaler;
    }
    if ($sample) {
	$self->{SAMPLING} = $sample;
    }

    if ($self->{NR_CLASSES} != $candidates->{NR_CLASSES}) {
	$self->die("Candidates have the wrong number of events\n");
    }

    $self->scale();
    $kl = $self->{KL};
    $n = ($n > $candidates->{NR_CANDIDATES}) ? $candidates->{NR_CANDIDATES} : $n;
    for ($i = 0; $i < $n; $i++) {
	$self->gain($candidates);
	$self->add_feature($candidates, $candidates->{BEST_CAND});
	$candidates->{ADDED}{$candidates->{BEST_CAND}} = 1;
	$self->log("Adding candidate $candidates->{BEST_CAND}\n");
	$self->scale();
	$self->log("Actual gain: " . ($self->{KL} - $kl) . "\n");
	$kl = $self->{KL};
    }
    return(1);
}


1;

__END__


# Below is the stub of documentation for your module. You better edit it!

=head1 NAME

MaxEntropy - Perl5 module for Maximum Entropy Modeling and Feature Induction

=head1 SYNOPSIS

  use Statistics::MaxEntropy;

  # debugging messages; default 0
  $Statistics::MaxEntropy::debug = 0;

  # maximum number of iterations for IIS; default 100
  $Statistics::MaxEntropy::NEWTON_max_it = 100;

  # minimal distance between new and old x for Newton's method; 
  # default 0.001
  $Statistics::MaxEntropy::NEWTON_min = 0.001;

  # maximum number of iterations for Newton's method; default 100
  $Statistics::MaxEntropy::KL_max_it = 100;

  # minimal distance between new and old x; default 0.001
  $Statistics::MaxEntropy::KL_min = 0.001;

  # the size of Monte Carlo samples; default 1000
  $Statistics::MaxEntropy::SAMPLE_size = 1000;

  # creation of a new event space from an events file
  $events = Statistics::MaxEntropy::new($file);

  # Generalised Iterative Scaling, "corpus" means no sampling
  $events->scale("corpus", "gis");

  # Improved Iterative Scaling, "mc" means Monte Carlo sampling
  $events->scale("mc", "iis");

  # Feature Induction algorithm, also see Statistics::Candidates POD
  $candidates = Statistics::Candidates->new($candidates_file);
  $events->fi("iis", $candidates, $nr_to_add, "mc");

  # writing new events, candidates, and parameters files
  $events->write($some_other_file);
  $events->write_parameters($file);
  $events->write_parameters_with_names($file);

  # dump/undump the event space to/from a file
  $events->dump($file);
  $events->undump($file);


=head1 DESCRIPTION

This module is an implementation of the Generalised and Improved
Iterative Scaling (GIS, IIS) algorithms and the Feature Induction (FI)
algorithm as defined in (B<Darroch and Ratcliff 1972>) and (B<Della
Pietra et al. 1997>). The purpose of the scaling algorithms is to find
the maximum entropy distribution given a set of events and
(optionally) an initial distribution. Also a set of candidate features
may be specified; then the FI algorithm may be applied to find and add
the candidate feature(s) that give the largest `gain' in terms of
Kullback Leibler divergence when it is added to the current set of
features.

Events are specified in terms of a set of feature functions
(properties) f_1...f_k that map each event to {0,1}: an event is a
string of bits. In addition of each event its frequency is given. We
assume the event space to have a probability distribution that can be
described by

=begin roff

    p(x) = 1/Z e^{sum_i alpha_i f_i(x)}

=end roff

=begin text

    p(x) = 1/Z e^{sum_i alpha_i f_i(x)}

=end text

=begin latex

\begin{equation*}
    p(x) = \frac{1}{Z} \exp[\sum_i \alpha_i f_i(x)]
\end{equation*}
where $Z$ is a normalisation factor given by
\begin{equation*}
    Z = \sum_x \exp[\sum_i \alpha_i f_i(x)]
\end{equation*}

=end latex

=begin roff

where Z is a normalisation factor. The purpose of the IIS algorithm is
the find alpha_1..alpha_k such that D(p~||p), defined by

    D(p~||p) = 
       sum_x p~ . log(p~(x) / p(x)),

is minimal under the condition that p~[f_i] = p[f_i], for all i.

=end roff

=begin latex

The purpose of the scaling algorithms IIS GIS is
the find $\alpha_1..\alpha_k$ such that $D(\tilde{p}||p)$, defined by
\begin{equation*}
    D(\tilde{p}||p) = 
       \sum_x \tilde{p} \log (\frac{\tilde{p}(x)}{p(x)}),
\end{equation*}
is minimal under the condition that for all $i$
$\tilde{p}[f_i]=p[f_i]$.

=end latex

The module requires the C<Bit::SparseVector> module by Steffen Beyer and the
C<Data::Dumper> module by Gurusamy Sarathy. Both can be obtained from
CPAN just like this module.


=head2 CONFIGURATION VARIABLES

=over 4

=item C<$Statistics::MaxEntropy::debug>

If set to C<1>, lots of debug information, and intermediate results will be
output. Default: C<0>

=item C<$Statistics::MaxEntropy::NEWTON_max_it>

Sets the maximum number of iterations in Newton's method. Newton's
method is applied to find the new parameters \alpha_i of the features
C<f_i>. Default: C<100>.

=item C<$Statistics::MaxEntropy::NEWTON_min>

Sets the minimum difference between x' and x in Newton's method (used for
computing parameter updates in IIS); if either the maximum number of
iterations is reached or the difference between x' and x is small enough,
the iteration is stopped. Default: C<0.001>. Sometimes features have
Infinity or -Infinity as a solution; these features are excluded from future
iterations.

=item C<$Statistics::MaxEntropy::KL_max_it>

Sets the maximum number of iterations applied in the IIS
algorithm. Default: C<100>.

=item C<$Statistics::MaxEntropy::KL_min>

Sets the minimum difference between KL divergences of two
distributions in the IIS algorithm; if either the maximum number of
iterations is reached or the difference between the divergences is
enough, the iteration is stopped. Default: C<0.001>.

=item C<$Statistics::MaxEntropy::SAMPLE_size>

Determines the number of (unique) events a sample should contain. Only
makes sense if for sampling "mc" is selected (see below). Its default
is C<1000>.

=back


=head2 METHODS

=over 4

=item C<new>

 $events = Statistics::MaxEntropy::new($events_file);

A new event space is created, and the events are read from C<$file>. The
events file is required, its syntax is described in 
L<FILE SYNTAX>.

=item C<write>

 $events->write($file);

Writes the events to a file. Its syntax is described in 
L<FILE SYNTAX>.

=item C<scale>

 $events->scale($sample, $scaler);

If C<$scaler> equals C<"gis">, the Generalised Iterative Scaling algorithm
(B<Darroch and Ratcliff 1972>) is applied on the event space; C<$scaler>
equals C<"iis">, the Improved Iterative Scaling Algorithm (B<Della Pietra et
al. 1997>) is used. If C<$sample> is C<"corpus">, there is no sampling done
to re-estimate the parameters (the events previously read are considered a
good sample); if it equals C<"mc"> Monte Carlo (Metropolis-Hastings)
sampling is performed to obtain a random sample; if C<$sample> is C<"enum">
the complete event space is enumerated.

=item C<fi>

 fi($scaler, $candidates, $nr_to_add, $sampling);

Calls the Feature Induction algorithm. The parameter C<$nr_to_add> is for
the number of candidates it should add. If this number is greater than the
number of candidates, all candidates are added. Meaningfull values for
C<$scaler> are C<"gis"> and C<"iis">; default is C<"gis"> (see previous
item). C<$sampling> should be one of C<"corpus">, C<"mc">, C<"enum">.
C<$candidates> should be in the C<Statistics::Candidates> class:

 $candidates = Statistics::Candidates->new($file);

See L<Statistics::Candidates>.

=item C<write_parameters>

 $events->write_parameters($file);

=item C<write_parameters_with_names>

 $events->write_parameters_with_names($file);

=item C<dump>

 $events->dump($file);

C<$events> is written to C<$file> using C<Data::Dumper>.

=item C<undump>

 $events = Statistics::MaxEntropy->undump($file);

The contents of file C<$file> is read and eval'ed into C<$events>.

=back


=head1 FILE SYNTAX

Lines that start with a C<#> and empty lines are ignored.

Below we give the syntax of in and output files.


=head2 EVENTS FILE (input/output)

Syntax of the event file (C<n> features, and C<m> events); the following
holds for features:

=over 4

=item *

each line is an event; 

=item *

each column represents a feature function; the co-domain of a feature
function is {0,1};

=item *

no space between feature columns; 

=item *

constant features (i.e. columns that are completely 0 or 1) are
forbidden;

=item *

2 or more events should be specified (this is in fact a consequence of
the previous requirement;

=back

The frequency of each event precedes the feature columns. Features are
indexed from right to left. This is a consequence of how
C<Bit::SparseVector> reads bit strings. Each C<f_ij> is a bit and C<freq_i>
an integer in the following schema:

    name_n <tab> name_n-1 ... name_2 <tab> name_1 <newline>
    freq_1 <white> f_1n ... f_13 f_12 f_11 <newline>
      .                     .
      .                     .
      .                     .
    freq_i <white> f_in ... f_i3 f_i2 f_i1 <newline>
      .                     .
      .                     .
      .                     .
    freq_m <white> f_mn ... f_m3 f_m2 f_m1

(C<m> events, C<n> features) The feature names are separated by tabs,
not white space. The line containing the feature names will be split
on tabs; this implies that (non-tab) white space may be part of the
feature names.


=head2 PARAMETERS FILE (input/output)

Syntax of the initial parameters file; one parameter per line:

    par_1 <newline>
     .
     .
     .
    par_i <newline>
     .
     .
     .
    par_n

The syntax of the output distribution is the same. The alternative
procedure for saving parameters to a file
C<write_parameters_with_names> writes files that have the following
syntax

    n <newline>
    name_1 <tab> par_1 <newline>
     .
     .
     .
    name_i <tab> par_i <newline>
     .
     .
     .
    name_n <tab> par_n <newline>
    bitmask

where bitmask can be used to tell other programs what features to use
in computing probabilities. Features that were ignored during scaling
or because they are constant functions, receive a C<0> bit.


=head2 DUMP FILE (input/output)

A dump file contains the event space (which is a hash blessed into
class C<Statistics::MaxEntropy>) as a Perl expression that can be
evaluated with eval.


=head1 BUGS

It's slow.


=head1 SEE ALSO

L<perl(1)>, L<Statistics::Candidates>, L<Statistics::SparseVector>,
L<Bit::Vector>, L<Data::Dumper>, L<POSIX>, L<Carp>.


=head1 DIAGNOSTICS

The module dies with an appropriate message if

=over 4

=item *

it cannot open a specified events file;

=item *

if you specified a constant feature function (in the events file or
the candidates file);

=item *

if the events file, candidates file, or the parameters file is not
consistent; possible causes are (a.o.): insufficient or too many
features for some event; inconsistent candidate lines; insufficient,
or to many event lines in the candidates file.

=back

The module captures C<SIGQUIT> and C<SIGINT>. On a C<SIGINT>
(typically <CONTROL-C> it will dump the current event space(s) and
die. If a C<SIGQUIT> (<CONTROL-BACKSLASH>) occurs it dumps the current
event space as soon as possible after the first iteration it finishes.


=head1 REFERENCES

=over 4

=item (Darroch and Ratcliff 1972) 

J. Darroch and D. Ratcliff, Generalised Iterative Scaling for
log-linear models, Ann. Math. Statist., 43, 1470-1480, 1972.

=item (Jaynes 1983)

E.T. Jaynes, Papers on probability, statistics, and statistical
physics. Ed.: R.D. Rosenkrantz. Kluwer Academic Publishers, 1983.

=item (Jaynes 1997) 

E.T. Jaynes, Probability theory: the logic of science, 1997,
unpublished manuscript.
C<URL:http://omega.math.albany.edu:8008/JaynesBook.html>

=item (Della Pietra et al. 1997) 

Stephen Della Pietra, Vincent Della Pietra, and John Lafferty,
Inducing features of random fields, In: Transactions Pattern Analysis
and Machine Intelligence, 19(4), April 1997.

=back


=head1 VERSION

Version 0.7.


=head1 AUTHOR

=begin roff

Hugo WL ter Doest, terdoest@cs.utwente.nl

=end roff

=begin latex

Hugo WL ter Doest, \texttt{terdoest\symbol{'100}cs.utwente.nl}

=end latex


=head1 COPYRIGHT

=begin roff

Copyright (C) 1998 Hugo WL ter Doest, terdoest@cs.utwente.nl
Univ. of Twente, Dept. of Comp. Sc., Parlevink Research, Enschede,
The Netherlands.

=end roff

=begin latex

Copyright (C) 1998 Hugo WL ter Doest, \texttt{terdoest\symbol{'100}cs.utwente.nl}
Univ. of Twente, Dept. of Comp. Sc., Parlevink Research, Enschede,
The Netherlands.

=end latex

C<Statistics::MaxEntropy> comes with ABSOLUTELY NO WARRANTY and may be copied
only under the terms of the GNU Library General Public License (version 2, or
later), which may be found in the distribution.

=cut
