package Statistics::MaxEntropy;

##---------------------------------------------------------------------------##
##  Author:
##      Hugo WL ter Doest       terdoest@cs.utwente.nl
##  Description: Statistics::MaxEntropy
##	Improved Iterative Scaling algorithm and 
##      Feature Induction algorithm
##  Keywords:
##      Maximum Entropy Modeling
##      Kullback-Leibler Divergence
##      Exponential models
##
##---------------------------------------------------------------------------##
##  Copyright (C) 1998 Hugo WL ter Doest terdoest@cs.utwente.nl
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
##	Require libraries
##---------------------------------------------------------------------------##
use strict;
use diagnostics -verbose;
use Bit::Vector;
# for floor and ceil functions
use POSIX;


##---------------------------------------------------------------------------##
##	Globals
##---------------------------------------------------------------------------##
use vars qw($VERSION 
	    @ISA 
	    @EXPORT 
	    @events 
	    $nr_events
	    $nr_features
	    @candidates
	    $nr_candidates
	    @pars
	    @p
	    @p_ref
	    $Z
	    $M
	    @nr_feats_on
	    @a
	    $normalise
	    @E_ref
	    $debug
	    $NEWTON_max_it
	    $KL_max_it
	    $KL_min
	    $NEWTON_min
	    %is_added
	    $rand_max_features
	    $rand_max_candidates
	    @correction_feature
	    $correction_parameter
	    );

use subs qw(GIS
	    IIS
	    );

require Exporter;
require AutoLoader;

@ISA = qw(Exporter AutoLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw(
	     $KL_min
	     $NEWTON_min
	     $debug
	     $normalise
	     $nr_add
	     $KL_max_it
	     $NEWTON_max_it
	     $rand_max_features
	     $rand_max_candidates
	     init
	     random_init
	     done
	     IIS
	     GIS
	     FI
);

$VERSION = '0.5';
$normalise = 1;
$NEWTON_max_it = 200;
$NEWTON_min = 0.001;
$KL_max_it = 100;
$KL_min = 0.001;
$rand_max_candidates = 20;
$rand_max_features = 20;
$debug = 0;

##---------------------------------------------------------------------------##
##	Private routines
##---------------------------------------------------------------------------##


##	I/O routines

# reads an events file
# dies in case of inconsistent lines
sub read_events {
    my($file) = shift;
    my $features;
    my $sum = 0;
    my $i;

    open(EVENTS,$file) ||
	die "Could not open $file\n";
    print "Opened $file\n";
    $nr_events = 0;
    $nr_features = 0;
    while (<EVENTS>) {
	if (!/\#.*/) {
	    chomp;

	    ($p_ref[$nr_events],$features) = split;
	    $sum += $p_ref[$nr_events];

	    # if first event set nr_features
	    if ($nr_events == 0) {
		$nr_features = length($features);
	    }
	    # else check nr of features for this event
	    else {
		if (length($features) != $nr_features) {
		    die "Events file corrupt (line $nr_events)\n";
		}
	    }
	    # create and initialise bit vector
	    $events[$nr_events] = Bit::Vector->new_Bin($nr_features,$features);
	    $nr_events++;
	}
    }
    close(EVENTS);

    # normalise, if necessary
    if ($normalise) {
	normalise_p_ref();
    }

    print "Read $nr_events events and $nr_features features;\n";
    print "closed $file\n";
}


# normalises the reference distribution
sub normalise_p_ref {
    my ($i, $sum);

    $sum = 0;
    for ($i = 0; $i < $nr_events; $i++) {
	$sum += $p_ref[$i];
	if ($debug) {
	    print "prob event $i: $p_ref[$i]\n";
	}
    }
    for ($i = 0; $i < $nr_events; $i++) {
	$p_ref[$i] /= $sum;
	if ($debug) {
	    print "prob event $i: $p_ref[$i]\n";
	}
    }
}


# reads a candidates file
# dies if insufficient events or inconsistent lines
sub read_candidates {
    my($file) = shift;
    my $features;
    my $sum = 0;
    my $event;

    open(CANDS,$file) ||
	die "Could not open $file\n";
    print "Opened $file\n";
    $nr_candidates = 0;
    $event = 0;
    while (<CANDS>) {
	if (!/\#.*/) {
	    chomp;
	    $features = $_;
	    if ($event == 0) {
		$nr_candidates = length($features);
	    }
	    else {
		if ($nr_candidates != length($features)) {
		    die "Candidate file corrupt ".
			"(line $event has insufficient features)\n";
		}
	    }
	    $candidates[$event++] = 
	      Bit::Vector->new_Bin($nr_candidates,$features);
	}
    }
    close(CANDS);

    # check the candidates for constant functions
    check_candidates();

    print "Read $nr_candidates candidates; closed $file\n";
    if ($nr_events != $event) {
	die "Candidate file has insufficient events\n";
    }
}


# writes events to file (in case of new feature(s))
sub write_candidates {
    my($file) = shift;
    my($x,$f);

    open(CANDIDATES,">$file") ||
	die "Could not open $file\n";
    print "Opened $file\n";
    
    # write candidates that were not added
    for ($x = 0; $x < $nr_events; $x++) {
	for ($f = 0; $f < $nr_candidates; $f++) {
	    if (!$is_added{$f}) {
		print CANDIDATES $candidates[$x]->bit_test($f);
	    }
	}
	print CANDIDATES "\n";
    }
    close CANDIDATES;
    print "Closed $file\n";
}


# initialise parameters randomly
sub random_parameters {
    my $i;

    srand();
    for ($i = 0; $i < $nr_features; $i++) {
	$pars[$i] = rand();
    }
    if ($debug) {
	print_distr();
    }
}


# reads an initial distribution
sub read_parameters {
    my($file) = shift;
    my $i = 0;

    open(DISTR,$file) ||
	die "Could not open $file\n";
    print "Opened $file\n";
    while (<DISTR>) {
	if (!/\#.*/) {
	    chomp;
	    $pars[$i++] = $_;
	}
    }
    close(DISTR);
    print "Read $nr_features features; closed $file\n";
    if ($i != $nr_features) {
	die "Initial distribution file corrupt\n";
    }
}


# make sure \tilde{p} << q_0
sub check_initial_distr {
    my $i;

    for ($i = 0; $i < $nr_events; $i++) {
	if ((p($i) == 0) && (p_ref($i) != 0)) {
	    die "Initial distribution not ok!\n";
	}
    }
}


# constant feature functions are forbidden: that is why
# we check whether for all features \sum_x f(x) > 0
# and \sum_x f(x) != $nr_events
sub check_features {
    my($x,$f);
    my $sum = 0;

    for ($f = 0; $f < $nr_features; $f++) {
	$sum = 0;
	for ($x = 0; $x < $nr_events; $x++) {
	    $sum += $events[$x]->bit_test($f);
	}
	if (!$sum || ($sum == $nr_events)) {
	    die "Feature ",$f+1, " is constant, remove it\n";
	}
    }
}


# check whether for all features f, \sum_x f(x) > 0
sub check_candidates {
    my($x,$f);
    my $sum = 0;

    for ($f = 0; $f < $nr_candidates; $f++) {
	$sum = 0;
	for ($x = 0; $x < $nr_events; $x++) {
	    $sum += $candidates[$x]->bit_test($f);
	}
	if (!$sum || ($sum == $nr_events)) {
	    die "Candidate feature ",$f+1, " is constant, remove it\n";
	}
    }
}


# writes the the current parameters to a file
sub write_parameters {
    my($file) = shift;

    open(DISTR,">$file") ||
	die "Could not open $file\n";
    print "Opened $file\n";
    for (@pars) {
	print DISTR "$_\n";
    }
    close(DISTR);
    print "Closed $file\n";
}


# writes events to a file 
# usefull in case new features have been added
sub write_events {
    my($file) = shift;
    my($x,$f);

    open(EVENTS,">$file") ||
	die "Could not open $file\n";
    print "Opened $file\n";
    for ($x = 0; $x < $nr_events; $x++) {
	print EVENTS $events[$x]->to_Bin();
	print EVENTS " $p_ref[$x]\n";
    }
    close EVENTS;
    print "Closed $file\n";
}


# prints the parameters of the current distr to the screen
sub print_distr {
    my $i = 0;

    for (@pars) {
	print "$i: $_\n";
	$i++
    }
}


# writes the added features to a file
# in the feature it might even more (sic) info...
sub write_info {
    my($file) = @_;

    my $f;

    open(INFO,">$file") ||
	die "Could not open $file\n";
    print "Opened $file\n";
    if ($nr_candidates) {
	print INFO "Added ";
	for ($f = 0; $f < $nr_candidates; $f++) {
	    if ($is_added{$f}) {
		print INFO "$f\t";
	    }
	}
    }
    print INFO "\n";
    print INFO "p(events) ", p_events(), "\n";
    print INFO "Z $Z\n", ;
    close(INFO);
    print "Closed $file\n";
}


##	Computation routines

# computes \sum_{x\in\Omega} p(x)g(x)
# functions $p and $g receive event number!
sub E {
    my($p) = shift;
    my($g) = shift;

    my $sum = 0;
    my $i;

    for ($i = 0; $i < $nr_events; $i++) {
	$sum += &$g($i) * &$p($i);
    }
    return($sum);
}


# returns p^{(n)}(x)
sub p {
    my($x) = shift;
    
    return $p[$x];
}


# (re)computes the current distribution p^{(n)}
sub curr_distr {
    my($corr) = @_;

    my ($x, $f);
    my $sum;

    $Z = 0;
    for ($x = 0; $x < $nr_events; $x++) {    
	$sum = 0;
	for ($f = 0; $f < $nr_features; $f++) {
	    $sum += $events[$x]->bit_test($f) * $pars[$f];
	}
	if ($corr) {
	    $sum += ($M - $nr_feats_on[$x]) * $correction_parameter;
	}
	$p[$x] = exp($sum);
	$Z += $p[$x];
    }
    # normalise
    for ($x = 0; $x < $nr_events; $x++) {    
	$p[$x] /= $Z;
    }
}


# gives reference probability for an event
sub p_ref {
    my($event) = shift;

    return($p_ref[$event]);
}


# compute a_{m,j}^{(n)}
sub a_mj {
    my($m) = shift;
    my($j) = shift;

    my $func;

    if ($m == 0) {
	return(-$E_ref[$j]);
    }
    else {
	$func = sub {
	    my($event) = shift;
	    
	    # f_j * \delta(m,f_#)
	    if ($m == $nr_feats_on[$event]) {
		return($events[$event]->bit_test($j));
	    }
	    else {
		return(0);
	    }
	};
	return(E(\&p,$func));
    }
}


# computes the `a' coefficients of 
# \sum_{m=0}^{M} a_{m,j}^{(n)} e^{\alpha^{(n)}_j m}
# according to the current distribution
sub a {
    my($j, $m);

    for  ($j = 0; $j < $nr_features; $j++) {
	for ($m = 0; $m <= $M; $m++) {
	    $a[$m][$j] = a_mj($m, $j);
	    if ($debug) {
		print "a[$m][$j] = $a[$m][$j]\n";
	    }
	}
    }
}


# computes f_# for alle events
# results in $nr_feats_on
sub compute_nr_feats_on {
    my($i, $j);

    for ($i = 0; $i < $nr_events; $i++) {
	$nr_feats_on[$i] = 0;
	for ($j = 0; $j < $nr_features; $j++) {
	    $nr_feats_on[$i] += $events[$i]->bit_test($j);
	}
	if ($debug) {
	    print "f_#($i) = $nr_feats_on[$i]\n";
	}
    }
}


# the maximum number of feature that some event has on: M
sub M {
    my $i;

    $M = 0;
    for ($i = 0;$i < $nr_events; $i++) {
	if ($nr_feats_on[$i] > $M) {
	    $M = $nr_feats_on[$i];
	}
    }
    if ($debug) {
	print "M = $M\n";
    }
}


# separate procedure; called if one candidate is added
sub E_ref_i {
    my($i) = @_;
    my $func;

    $func = sub {
	my ($event) = shift;
	
	return($events[$event]->bit_test($i));
    };
    $E_ref[$i] = E(\&p_ref,$func);
}


# fills in $E_ref of the features
sub E_ref {
    my($i);
    my $func;

    for $i (0..$nr_features-1) {
	E_ref_i($i);
    }
}


sub a_func {
    my($x) = shift;
    my($j) = shift;
    my($m, $sum);

    $sum = 0;
    for ($m = 0; $m <= $M; $m++) {
	$sum += $a[$m][$j] * exp($x * $m);
    }
    return($sum);
}


sub a_deriv {
    my($x) = shift;
    my($j) = shift;
    my($m, $sum);

    $sum = 0;
    for ($m = 1; $m <= $M; $m++) {
	$sum += $a[$m][$j] * $m * exp($x * $m);
    }
    return($sum);
}


# solves \alpha from 
# \sum_{m=0}^{M} a_{m,j}^{(n)} e^{\alpha^{(n)}_j m}=0
sub iis_estimate_with_newton {
    my($i) = shift;
    my($x, $old_x);
    my($deriv_res);
    my $diff;
    my $k = 0;

    # initialise x such that a_deriv(x) != 0
    $x = 10;

    # do newton's method
    do {
	# save old x
	$old_x = $x;
	# compute new x
	$deriv_res = a_deriv($x,$i);
	if ($deriv_res == 0) {
	    print STDERR "Derivative of $i is zero\n";
	    return(0);
	}
	$x = $x - a_func($x,$i) / $deriv_res;
	if ($debug) {
	    printf("x = %e\n",$x);
	}
	$diff = abs($x - $old_x);
	$k++;
    } until (($diff <= $NEWTON_min) ||
	     ($k > $NEWTON_max_it));
    if ($debug) {
	print "Estimated new alpha_$i with Newton's method\n";
    }
    return($x);
}


# determines Kullback Leibler divergence between reference distribution
# and current distribution
# takes two function arguments
sub KL {
    my($p1,$p2) = @_;

    my ($i, $sum);

    $sum = 0;
    for ($i = 0; $i < $nr_events; $i++) {
	$sum += &$p1($i) * log(&$p1($i) / &$p2($i));
    }
    return($sum);
}


# returns entropy of p
sub H {
    my ($i, $sum);

    $sum = 0;
    for ($i = 0; $i < $nr_events; $i++) {
	$sum += p($i) * log(p($i));
    }
    return(-$sum);
}


# cross entropy
sub cross_H {
    my ($i, $sum);

    $sum = 0;
    for ($i = 0; $i < $nr_events; $i++) {
	$sum += p_ref($i) * log(p($i));
    }
    return(-$sum);
}


# likelihood measure
sub L {
    return(-cross_H());
}


# sums the probabilities over all events (should be one)
sub check_p {
    my ($i, $sum);

    for ($i = 0; $i < $nr_events; $i++) {
	$sum += p($i)
    }
    print "Sum of final distr: $sum\n";
}


sub p_events {
    my ($i, $prod);

    $prod =1;
    for ($i = 0; $i < $nr_events; $i++) {
	$prod *= p($i);
    }
    return($prod);
}


#
# Field induction, given a set of candidates and the current 
# distribution p, the best candidate is found and added to the model
#

# computes the gain G_p(g) of a candidate feature
sub G {
    my($g, $is) = @_;

    my $func;
    my $i;
    my $new_kl;
    my $max;
    my $below;
    my $above;
    my $E_p_ref_g;
    my $E_p_g;

    $func = sub {
	my($event) = shift;

	return($candidates[$event]->bit_test($g));
    };
    $E_p_ref_g = E(\&p_ref,$func);
    $E_p_g = E(\&p,$func);
    $above = $E_p_ref_g * (1 - $E_p_g);
    $below = $E_p_g * (1 - $E_p_ref_g);
    if ((($above > 0) && ($below > 0)) ||
	(($above < 0) && ($below < 0))) {
	$max = log(abs($above)) - log(abs($below));
    }
    else {
	die "Cannot take log of negative or zero value: $above / $below\n";
    }

    # temporarily add feature to field
    $nr_features++;
    for ($i = 0; $i < $nr_events; $i++) {
	$events[$i]->Interval_Substitute($candidates[$i],
					       $events[$i]->Size(),
					       0,
					       $g,
					       1);
    }
    $pars[$nr_features - 1] = $max;

    # renormalise and compute KL
    curr_distr($is == \&GIS);
    $new_kl = KL(\&p_ref,\&p);

    # restore old field and renormalise
    for ($i = 0; $i < $nr_events; $i++) {
	$events[$i]->Resize($events[$i]->Size() - 1);
    }
    undef $pars[$nr_features - 1];
    $nr_features--;
    curr_distr($is == \&GIS);

    # return gain
    return(KL(\&p_ref,\&p) - $new_kl, $max);
}


#
# Improved Iterative Scaling Algorithm
# 
sub iis {
    my $k = 0;
    my $i;
    my $kl = 1e99;
    my $old_kl;
    my $diff;
    
    print "IIS:\n";
    print "it.\tD(p_ref||p)\tL(p)\t\tp(events)\n";
    printf "0\t%e\t%e\t%e\n", KL(\&p_ref,\&p), L(), p_events();
    do {
	a();
	for ($i = 0; $i < $nr_features; $i++) {
	    $pars[$i] += iis_estimate_with_newton($i);
	}
	$old_kl = $kl;
	curr_distr(0);
	$kl = KL(\&p_ref,\&p);
	$k++;
	printf "%u\t%e\t%e\t%e\n", $k, $kl, L(), p_events();
	if ($debug) {
	    print_distr();
	}
	$diff = $old_kl - $kl;
	if ($diff < 0) {
	    die "Something is very wrong; IIS is not converging!\n";
	}
    } until ($diff <= $KL_min || 
	     ($k > $KL_max_it));
    if ($debug) {
	curr_distr(0);
	check_p();
    }
}


sub correction_feature {
    my($event) = @_;
    
    return($M - $nr_feats_on[$event]);
}


#
# Generalised Iterative Scaling Algorithm
# 
sub gis {
    my $k = 0;
    my $i;
    my $kl = 1e99;
    my $old_kl;
    my $diff;
    my $func;
    

    print "GIS:\n";
    print "it.\tD(p_ref||p)\tL(p)\t\tp(events)\n";
    printf "0\t%e\t%e\t%e\n", KL(\&p_ref,\&p), L(), p_events();
    for $i (0..$nr_features-1) {
	$pars[$i] = 1;
    }
    $correction_parameter = 1;
    curr_distr(1);
    do {
	# binary features
	for $i (0..$nr_features-1) {
	    $func = sub {
		my ($event) = @_;

		return($events[$event]->bit_test($i));
	    };
#	    print "\$E_ref[$i] : $E_ref[$i]\n";
#	    print "\$pars[\$i]: $pars[$i]\n";
	    $pars[$i] *= exp((log($E_ref[$i]) - log(E(\&p,$func))) / $M);
	}

	# correction feature parameter estimation
	$correction_parameter *= 
	    exp((log(E(\&p_ref,\&correction_feature)) - 
		 log(E(\&p,\&correction_feature))) / $M);

	$old_kl = $kl;
	curr_distr(1);
	$kl = KL(\&p_ref,\&p);
	$k++;
	printf "%u\t%e\t%e\t%e\n", $k, $kl, L(), p_events();
	if ($debug) {
	    print_distr();
	}
	$diff = $old_kl - $kl;
	if ($diff < 0) {
	    die "Something is very wrong; GIS is not converging!\n";
	}
    } until ($diff <= $KL_min || 
	     ($k > $KL_max_it));
    if ($debug) {
	curr_distr(1);
	check_p();
    }
}


#
# Field Induction Algorithm
#
# adds the best candidate!
sub fia {
    my($is) = @_;

    my ($i, $gain, $par);
    my $best_gain = 0;
    my $best_par = 0;
    my $best_i = 0;

    # find the best candidate
    for ($i = 0; $i < $nr_candidates; $i++) {
	# check if not already added
	if (!$is_added{$i}) {
	    ($gain,$par) = G($i, $is);
	    printf "%u\t gain: %e\talpha: %e\n", $i, $gain, $par;
	    if ($best_gain < $gain) {
		$best_gain = $gain;
		$best_i = $i;
		$best_par = $par;
	    }
	}
    }

    # extend the events with the best candidate
    print "Adding candidate $best_i\n";
    $nr_features++;
    for ($i = 0; $i < $nr_events; $i++) {
	$events[$i]->Interval_Substitute($candidates[$i],
					 $events[$i]->Size(),
					 0,
					 $best_i,
					 1);
    }
    $pars[$nr_features - 1] = $best_par;
    $is_added{$best_i} = 1;
    M();
    E_ref_i($nr_features-1);
    &$is();
}


##---------------------------------------------------------------------------##
##	Public routines
##---------------------------------------------------------------------------##

# read events, parameters, and candidates
sub init {
    my($e,$p,$c) = @_;
    
    if ($e) {
	read_events($e);
    }
    else {
	die "Event file is required\n";
    }
    if ($p) {
	read_parameters($p);
	$correction_parameter = 1;
    }
    else {
	random_parameters();
    }
    if ($c) {
	read_candidates($c);
    }

    # some computations
    compute_nr_feats_on();
    M();
    E_ref();
    curr_distr(0);

    # some checks
    check_initial_distr();
    check_features();
}


# create random event, parameter and candidate files
sub random_init {
    my($e,$p,$c) = @_;

    my($x,
       $f,
       %no_double_events,
       $vec,
       $freq);

    # randomise
    srand();

    # open the files
    if ($e) {
	open(EVENTS,">$e") ||
	    die "Could not open $e\n";
    }
    if ($c) {
	open(CANDS,">$c") ||
	    die "Could not open $c\n";
    }
    if ($p) {
	open(DISTR,">$p")  ||
	    die "Could not open $p\n";
    }

    $nr_features = ceil(rand($rand_max_features));
    print "Number of features: $nr_features\n";

    # greater than 1 AND smaller than (2^$nr_features - 1)
    $nr_events = 100;#1 + ceil(rand(1 << ($nr_features - 1)));
    print "Number of events: $nr_events\n";

    $nr_candidates = ceil(rand($rand_max_candidates));
    print "Number of candidates: $nr_candidates\n";

    # events
    $x = 0;
    while ($x < $nr_events) {
	$freq = 1 + floor(rand(100));
	$vec = Bit::Vector->new($nr_features);
	for ($f = 0; $f < $nr_features; $f++) {
	    if (rand() < 0.5) {
		$vec->Bit_On($f);
	    }
	}
	if (!$no_double_events{$vec->to_Bin}) {
	    $no_double_events{$vec->to_Bin()} = $x;
	    $p_ref[$x] = $freq;
	    $events[$x] = $vec->Clone();
	    $x++;
	}
	else {
	    $p_ref[$no_double_events{$vec->to_Bin()}] += $freq;
	}
    }
    for ($x = 0; $x < $nr_events; $x++) {
	print EVENTS $p_ref[$x]," ";
	print EVENTS $events[$x]->to_Bin(), "\n";
    }

    # candidates
    for ($x = 0; $x < $nr_events; $x++) {
	$candidates[$x] = Bit::Vector->new($nr_candidates);
	for ($f = 0; $f < $nr_candidates; $f++) {
	    if (rand() < 0.5) {
		$candidates[$x]->Bit_On($f);
	    }
	}
	print CANDS $candidates[$x]->to_Bin,"\n";
    }

    # throw away double administration
    undef %no_double_events;

    # parameters
    for ($f = 0; $f < $nr_features; $f++) {
	$pars[$f] = 1; #rand();
	print DISTR $pars[$f],"\n";
    }
    $correction_parameter = 1;
    
    # close the files
    close(EVENTS);
    close(DISTR);
    close(CANDS);    

    # some computations
    if ($normalise) {
	normalise_p_ref();
    }
    compute_nr_feats_on();
    M();
    E_ref();
    curr_distr(0);

    # some checks
    check_initial_distr();
    check_features();
}


# wrapper around fia()
# dies if no candidates are given
sub FI {
    my($nr_to_add, $is) = @_;

    my $i;

    print "FI: $nr_to_add\n";
    if ($nr_candidates == 0) {
	die "Don't have candidates\n";
    }
    if (!$nr_to_add) {
	$nr_to_add = 1;
    }
    if ($nr_to_add > $nr_candidates) {
	$nr_to_add = $nr_candidates;
    }
    if (!$is) {
	$is = \&IIS;
    }
    if ($nr_candidates > 0) {
	# fia() assumes scaling is done for current features
	curr_distr($is == \&GIS);
	&$is();
	for ($i = 0; $i < $nr_to_add; $i++) {
	    fia($is);
	}
	# scale the last time;
	# reference distribution of the last added feature
	&$is();
    }
}


# wrapper around iis()
sub IIS {
   if ($nr_events > 0) {
       iis();
   }
   else {
       die "Don't have events\n";
   }
}


# wrapper around gis()
sub GIS {
    my $i;

    if ($nr_events > 0) {
	gis();
    }
    else {
	die "Don't have events\n";
    }
}


# write new events, parameters, and candidates
# and clean up the mess
sub done {
   my($e,$p,$c,$l) = @_;

   if ($p) {
       write_parameters($p);
   }
   if ($e) {
       write_events($e);
   }
   if ($c) {
       write_candidates($c);
   }
   if ($l) {
       write_info($l);
   }
   undef @pars;
   undef @events;
   undef @candidates;
   undef $nr_candidates;
   undef $nr_events;
   undef $nr_features;
   undef @a;
   undef @p;
   undef @p_ref;
   undef @E_ref;
   undef %is_added;
   undef $correction_parameter;
}


# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is the stub of documentation for your module. You better edit it!

=head1 NAME

MaxEntropy - Perl module for Maximum Entropy Modeling

=head1 SYNOPSIS

  use Statistics::MaxEntropy;

  # debugging messages; default 0
  $Statistics::MaxEntropy::debug = 0;

  # normalise frequencies in events file; default 1
  $Statistics::MaxEntropy::normalise = 1;

  # maximum number of iterations for IIS; default 100
  $Statistics::MaxEntropy::NEWTON_max_it = 100;

  # minimal distance between new and old x for Newton's method; 
  # default 0.001
  $Statistics::MaxEntropy::NEWTON_min = 0.001;

  # maximum number of iterations for Newton's method; default 100
  $Statistics::MaxEntropy::KL_max_it = 100;

  # minimal distance between new and old x; default 0.001
  $Statistics::MaxEntropy::KL_min = 0.001;

  # configuration of Statistics::MaxEntropy::random_init
  # maximum number of candidates randomly generated
  $Statistics::MaxEntropy::rand_max_candidates = 20;

  # maximum number of features randomly generated
  $Statistics::MaxEntropy::rand_max_features = 20;

  # initialisation for GIS, IIS, or FI
  Statistics::MaxEntropy::init($events_file,
		       $parameters_file,
		       $candidates_file);

  # random initialisation for GIS, IIS, or FI
  Statistics::MaxEntropy::random_init($random_events_file,
                              $random_parameters_file,
                              $random_candidates_file);

  # Generalised Iterative Scaling
  Statistics::MaxEntropy::GIS();

  # Improved Iterative Scaling
  Statistics::MaxEntropy::IIS();

  # Feature Induction algorithm
  Statistics::MaxEntropy::FI($nr_candidates_to_add, $iterative_scaler);

  # writes new events, candidates, and parameters files
  Statistics::MaxEntropy::done($new_events_file,
                       $new_parameters_file,
        	       $new_candidates_file,
                       $information_file);


=head1 DESCRIPTION

This module is an implementation of the Generalised and Improved
Iterative Scaling (GIS, IIS) algorithms and the Feature Induction (FI)
algorithm as defined in (B<Darroch and Ratcliff 1972>, B<(Della Pietra
et al. 1997)>). The purpose of the scaling algorithms is to find the
maximum entropy distribution given a set of events and (optionally) an
initial distribution. Also a set of candidate features may be
specified; then the FI algorithm may be applied to find and add the
candidate feature(s) that give the largest `gain' in terms of Kullback
Leibler divergence when it is added to the current set of features.

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
    p(x) = \frac{1}{Z} e^{\sum_i \alpha_i f_i(x)}
\end{equation*}
where $Z$ is a normalisation factor given by
\begin{equation*}
    Z = \sum_x e^{\sum_i \alpha_i f_i(x)}
\end{equation*}

=end latex

=begin roff

where Z is a normalisation factor. The purpose of the IIS algorithm is
the find alpha_1..alpha_k such that D(p~||p), defined by

    D(p~||p) = 
       sum_x p~ . log(p~(x) / p(x)),

is minimal

=end roff

=begin latex


The purpose of the IIS algorithm is
the find $\alpha_1..\alpha_k$ such that $D(\tilde{p}||p)$, defined by


\begin{equation*}
    D(\tilde{p}||p) = 
       \sum_x \tilde{p} . \log (\frac{\tilde{p}(x)}{p(x)}),
\end{equation*}
is minimal.

=end latex



=head2 VARIABLES

=over 4

=item *

C<$Statistics::MaxEntropy::debug>

If set to "1", lots of debug information, and intermediate results will
be output. Default: "0".

=item *

C<$Statistics::MaxEntropy::normalise>

If set to "1", frequencies in the events file are normalised; this is
required if they not sum to "1"; default: "1".

=item *

C<$Statistics::MaxEntropy::NEWTON_max_it>

Sets the maximum number of iterations in Newton's method. Newton's
method is applied to find the new parameters \alpha_i of the features
f_i. Default: "100".

=item *

C<$Statistics::MaxEntropy::NEWTON_min>

Sets the minimum difference between x' and x in Newton's method; if
either the maximum number of iterations is reached or the difference
between x' and x is small enough, the iteration is stopped. Default:
"0.001".

=item *

C<$Statistics::MaxEntropy::KL_max_it>

Sets the maximum number of iterations applied in the IIS
algorithm. Default: "100".

=item *

C<$Statistics::MaxEntropy::KL_min>

Sets the minimum difference between KL divergences of two
distributions in the IIS algorithm; if either the maximum number of
iterations is reached or the difference between the divergences is
enough, the iteration is stopped. Default: "0.001".

=item *

C<$Statistics::MaxEntropy::rand_max_candidates>

Configures C<Statistics::MaxEntropy::random_init>; it sets the maximum number
candidate features that are randomly generated for each event. Values
"2" and higher are allowed.

=item *

C<$Statistics::MaxEntropy::rand_max_features>

Configures C<Statistics::MaxEntropy::random_init>; it sets the maximum number
of features that are randomly generated for each event. 

The minimum number of events that are randomly generates is "2". The
maximum number of events is bounded from above by 2^{number of
features randomly chosen}.

=back


=head2 FUNCTIONS

=over 4

=item *

C<Statistics::MaxEntropy::init($events_file,>
C<                     $parameters_file,>
C<                     $candidates_file);>

The events file is required. The candidate and initial parameter files
are optional. If the initial parameter file is specified as the empty
string then random parameters from [0,1] are chosen.

=item *

C<Statistics::MaxEntropy::random_init($random_events_file,>
C<                            $random_parameters_file,>
C<                            $random_candidates_file);>

All three files are required. Your program dies, if you leave one or
more unspecified. A caveat: 

=item *

C<Statistics::MaxEntropy::GIS();>

Calls the Generalised Iterative Scaling algorithm. See 
(B<Darroch and Ratcliff 1972>).

=item *

C<Statistics::MaxEntropy::IIS();>

Calls the Improved Iterative Scaling algorithm. If no events were
loaded or generated, your program exits. See (B<Della Pietra et al. 1997>).

=item *

C<Statistics::MaxEntropy::FI($nr_candidates_to_add, 
                   $iterative_scaler);>

Calls the Feature Induction algorithm. The parameter
C<$nr_candidates_to_add> is for the number of candidates it should
add. If this number is greater than the number of candidates, all
candidates are added. If it is not specified, it "1" is used.  The
parameter C<$iterative_scaler> can be used to provide one of scaling
algorithm C<GIS> or C<IIS>. The references to the sub's should be
used as follows:

    $iterative_scaler = \&GIS;
    $nr_candidates_to_add = 2;
    FI($nr_candidates_to_add, $iterative_scaler);

If it is not specified, the improved algorithm is used.

=item *

C<Statistics::MaxEntropy::done($new_events_file,>
C<                     $new_parameters_file,>
C<                     $new_candidates_file,>
C<                     $information_file);>

All parameters are optional. The candidate numbers that were added are
printed to the information file. The new events file contains the
events extended with the new features, the new candidates file
contains the set of candidates minus the candidates added. To the
distribution file the new parameters are printed. The information file
is used to print the indices of the candidates that were added to.
Before a new call of C<init> or C<random_init>, calling
C<Statistics::MaxEntropy::done> is a good thing to do.

=back


=head1 SYNTAX OF INPUT/OUTPUT FILES

Lines that start with a `#' and empty lines are ignored.

Below we give the syntax of in and output files.


=head2 EVENTS FILE (input/output)

Syntax of the event file (n features, and m events); the following
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
C<Bit::Vector> reads bit strings. Each fij is a bit and freqi a float
in the following schema:

    freq1 <white> f1n ... f13 f12 f11 <newline>
      .                     .
      .                     .
      .                     .
    freqi <white> fin ... fi3 fi2 fi1 <newline>
      .                     .
      .                     .
      .                     .
    freqm <white> fmn ... fm3 fm2 fm1

(m events, n features)


=head2 PARAMETERS FILE (input/output)

Syntax of the initial parameters file; one parameter per line:

    par1 <newline>
     .
     .
     .
    pari <newline>
     .
     .
     .
    parn

(n features)

The syntax of the output distribution is the same.


=head2 CANDIDATES FILE (input)

The syntax of the candidate feature file is more or less the same as
that for the events file:

=over 4

=item *

each line is an event (events specified in the same order as the
events file);

=item *

each column is a feature;

=item * 

constant feature functions are forbidden;

=item *

values are 0 or 1; 

=item *

no space between features;

=back

fij are bits:

    f1c ... f13 f12 f11 <newline>
	     .
	     .
             .
    fic ... fi3 fi2 fi1 <newline>
	     .
             .
             .
    fmc ... fm3 fm2 fm1

(m events, c candidate features)

=head2 INFO FILE (output)

The indices of the candidates added to the field, one per line. Bit
strings are always indexed right to left.



=head1 BUGS

=over 4

=item *

It's slow.

=item *

C<Statistics::MaxEntropy> communicates through files mainly. The reason for is
a practical one: for my own purposes I needed communication through
files; I'm using the module for large event bases (corpora), and I'm
not interested in (large) arrays that tell me what candidates have
been added, what parameters were found, or how the events file looks
like in an array. My (other) software components will read the files
again.

=back


=head1 ENVIRONMENT

No environment variables are used.


=head1 AUTHOR

=begin roff

Hugo WL ter Doest, terdoest@cs.utwente.nl

=end roff

=begin latex

Hugo WL ter Doest, \texttt{terdoest\symbol{'100}cs.utwente.nl}

=end latex


=head1 SEE ALSO

L<perl(1)>, L<Bit::Vector>, and L<POSIX>.


=head1 DIAGNOSTICS

The module dies with an appropriate message if

=over 4

=item *

it cannot open a specified events file;

=item *

if you specified a constant feature function (in the events file or
the candidates file); this may happen sometimes if you call
C<Statistics::MaxEntropy::random_init>;

=item *

if C<Statistics::MaxEntropy::FI> is called and no candidates are given;

=item *

if C<Statistics::MaxEntropy:IIS> is called and no events are given;

=item *

if the events file, candidates file, or the parameters file is not
consistent; causes (a.o.): insufficient features, or too many features
for some event; inconsistent candidate lines; insufficient, or to many
event lines in the candidates file.

=back


=head1 REFERENCES

=over 4

=item *

(Darroch and Ratcliff 1972) J. Darroch and D. Ratcliff, Generalised
Iterative Scaling for log-linear models, Ann. Math. Statist., 43,
1470-1480, 1972.

=item *

(Jaynes 1997) E.T. Jaynes, Probability theory: the logic of science,
1997, unpublished manuscript.

=item *

(Della Pietra, Della Pietra and Lafferty 1997) Stephen Della Pietra,
Vincent Della Pietra, and John Lafferty, Inducing features of random
fields, In: Transactions Pattern Analysis and Machine Intelligence,
19(4), April 1997.

=back


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
only under the terms of the GNU General Public License (version 2, or
later), which may be found in the distribution.

=cut
