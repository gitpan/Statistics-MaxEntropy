package Statistics::Candidates;

##---------------------------------------------------------------------------##
##  Author:
##      Hugo WL ter Doest       terdoest@cs.utwente.nl
##  Description: 
##      Object/methods for candidate features
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
	    $VECTOR_PACKAGE
	    $SPARSE);


##---------------------------------------------------------------------------##
##	Require libraries
##---------------------------------------------------------------------------##
use strict;
use diagnostics -verbose;
use Carp;
$SPARSE = 1;
if ($SPARSE) {
    $VECTOR_PACKAGE = "Statistics::SparseVector";
    use Statistics::SparseVector;
}
else {
    $VECTOR_PACKAGE = "Bit::Vector";
    use Bit::Vector;
}
require Exporter;
require AutoLoader;

@ISA = qw(Exporter AutoLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw(new
	     DESTROY
	     write
	     );


sub new {
    my($this, $arg) = @_;

    # for calling $self->new($someth):
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;
    if ($arg) {
	$self->read($arg);
    }
    return($self);
}


sub DESTROY {

}


# reads a candidates file
# dies if insufficient events or inconsistent lines
# syntax first line: <name> <tab> <name> <tab> ... <newline>
# syntax other lines: <bitvector>
sub read {
    my($self, $file) = @_;

    my($features,
       $sum,
       $event,
       $candidate_names);

    # prologue
    open(CANDS, $file) ||
	die "Could not open $file\n";
    print "Opened $file\n";

    # read candidate names, skip comments
    $candidate_names = "";
    do {
	$candidate_names = <CANDS>;
    } until ($candidate_names !~ /\#.*/);
    chomp $candidate_names;
    $self->{CANDIDATE_NAMES} = [split(/\t/,$candidate_names)];

    # read the candidate bitvectors
    $self->{NR_CANDIDATES} = 0;
    $self->{NR_CLASSES} = 0;
    while (<CANDS>) {
	if (!/\#.*/) {
	    chomp;
	    $features = $_;
	    if ($self->{NR_CLASSES} == 0) {
		$self->{NR_CANDIDATES} = length($features);
	    }
	    else {
		if ($self->{NR_CANDIDATES} != length($features)) {
		    croak "Candidate file corrupt ".
			"(line $self->{NR_CLASSES} has insufficient features)\n";
		}
	    }
	    $self->{CANDIDATES}[$self->{NR_CLASSES}++] = 
	      $VECTOR_PACKAGE->new_Bin($self->{NR_CANDIDATES}, $features);
	}
    }

    # epilogue
    close(CANDS);
    # check the candidates for constant functions
    $self->check();
    print "Read $self->{NR_CANDIDATES} candidates for $self->{NR_CLASSES} events; ";
    print "closed $file\n";
}


# check whether for all features f, \sum_x f(x) > 0, and
# \sum_x f(x) != nr_classes
sub check {
    my($self) = @_;

    my($x,
       $f, 
       $sum);

    for ($f = 0; $f < $self->{NR_CANDIDATES}; $f++) {
	$sum = 0;
	for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
	    $sum += $self->{CANDIDATES}[$x]->bit_test($f);
	}
	if (!$sum || ($sum == $self->{NR_CLASSES})) {
	    croak "Candidate ",$f+1, " is constant, remove it\n";
	}
    }
}


# writes remaining candidates to a file
# syntax: same as input candidates file
sub write {
    my($self, $file) = @_;

    my($x,
       $f);

    if (($self->{NR_CANDIDATES} > 0) && ($self->{NR_CLASSES})) {
	open(CANDIDATES,">$file") ||
	    die "Could not open $file\n";
	print "Opened $file\n";

	# write the list of candidate names that were not added
	for ($f = 0; $f < $self->{NR_CANDIDATES}; $f++) {
	    if (!$self->{ADDED}{$f}) {
		print CANDIDATES "$self->{CANDIDATE_NAMES}[$f]\t";
	    }
	}
	print CANDIDATES "\n";

	# write candidates that were not added
	for ($x = 0; $x < $self->{NR_CLASSES}; $x++) {
	    for ($f = 0; $f < $self->{NR_CANDIDATES}; $f++) {
		if (!$self->{ADDED}{$f}) {
		    print CANDIDATES $self->{CANDIDATES}[$x]->bit_test($f);
		}
	    }
	    print CANDIDATES "\n";
	}

	close CANDIDATES;
	print "Closed $file\n";
    }
}


sub clear {
    my($self) = @_;

    undef $self->{ADDED};
}

1;

__END__


# Below is the stub of documentation for your module. You better edit it!

=head1 NAME

Candidates - Perl5 module for manipulating candidate features (help module for C<Statistics::MaxEntropy>)

=head1 SYNOPSIS

  use Statistics::Candidates;

  # create a new candidates object and read candidate features
  $candidates = Statistics::Candidates->new($some_file);

  # checks for constant candidate features
  $candidates->check();

  # writes candidates that were not added to a file
  $candidates->write($some_other_file);

  # clear the administration about being added or not ...
  $candidates->clear();


=head1 DESCRIPTION

The C<Candidates> object is for storage, retrieval, and manipulation
of candidate features.

The reason for separating this code from C<Maxentropy.pm> is that a
set of candidate features should be considered a separate
object. Blessing them into C<MaxEntropy> would be unnatural.

Also this code is simpler and more stable than that in the
C<MaxEntropy> module.

This module requires C<Bit::SparseVector>.


=head1 METHODS

=over 4

=item C<new>

 $candidates = Statistics::Candidates->new($file);

=item C<check>

 $candidates->check();

=item C<write>

 $candidates->write($file);

=item C<clear>

 $candidates->clear();

=back


=head1 FILE SYNTAX

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

=item *

lines that start with C<#> are ignored.

=back

Below is a set of candidates for C<m> events, C<c> candidate features;
C<f_ij> are bits:

    name_c <tab> name_c-1 ... name_2 <tab> name_1 <newline>
    f_1c ... f_13 f_12 f_11 <newline>
	     .
	     .
             .
    f_ic ... f_i3 f_i2 f_i1 <newline>
	     .
             .
             .
    f_mc ... f_m3 f_m2 f_m1




=head1 SEE ALSO

L<Statistics::MaxEntropy>, L<Statistics::SparseVector>,
L<Bit::Vector>.


=head1 VERSION

Version 0.1


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

C<Statistics::Candidates> comes with ABSOLUTELY NO WARRANTY and may be copied
only under the terms of the GNU Library General Public License (version 2, or
later), which may be found in the distribution.

=cut
