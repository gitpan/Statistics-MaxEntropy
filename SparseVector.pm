package Statistics::SparseVector;

##---------------------------------------------------------------------------##
##  Author:
##      Hugo WL ter Doest       terdoest@cs.utwente.nl
##  Description: module for sparse bitvectors
##               method names equal that of Bit::Vector
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


use strict;
use vars qw($VERSION 
	    @ISA 
	    @EXPORT 
	    @EXPORT_OK);
use overload
'++'     =>    \&increment,
'='      =>    \&Clone,
'""'     =>    \&stringify;
require Exporter;
require AutoLoader;

@ISA = qw(Exporter AutoLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw(
	
);
$VERSION = '0.1';


# Preloaded methods go here.

# necessary for overloading ""
# in turn required by Data::Dumper that wants to stringify variables
sub stringify {
    my($self) = @_;

    return($self);
}


# create a new bitvector
# all bits are `off'
sub new {
    my($this, $n) = @_;

    # for calling $self->new($someth):
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;
    $self->{N} = $n;
    return($self);
}


# create a copy of a bitvector and return it
sub Clone {
    my($self) = @_;

    my($new);

    $new = Statistics::SparseVector->new($self->{N});
    for (keys %{$self->{VECTOR}}) {
	$new->{VECTOR}->{$_} = $self->{VECTOR}->{$_};
    }
    return($new);
}


# is called right before Perl destroys a bitvector
# this happens automatically if a vector has zero references
sub DESTROY {
    my($self) = @_;

}


# creates a new vector from a bitstring
sub new_Bin {
    my($this, $n, $bitvector) = @_;

    my($i, $self);

    $self = $this->new($n);
    $i = 0;
    while ($bitvector) {
	if (chop($bitvector) eq "1") {
	    $self->{VECTOR}->{$i} = 1;
	}
	$i++;
    }
    return($self);
}


# we assume a bit is `on' if its value is defined
sub bit_test {
    my($self, $i) = @_;

    return(defined($self->{VECTOR}->{$i}));
}


# turning off a bit is making its value undefined
sub Bit_Off {
    my($self, $i) = @_;

    undef $self->{VECTOR}->{$i};
}


# turns on a bit
sub Bit_On {
    my($self, $i) = @_;

    $self->{VECTOR}->{$i} = 1;
}


# flips a bit, i.e. makes it undef if it is defined,
# and makes it defined if it is undefined
sub bit_flip {
    my($self, $i) = @_;

    if (defined($self->{VECTOR}->{$i})) {
	undef $self->{VECTOR}->{$i};
    }
    else {
	$self->{VECTOR}->{$i} = 1;
    }
}


# increases the integer value of the bitvector by one
sub increment {
    my($self) = @_;

    my($carry, $i);

    $carry = 0;
    $i = 0;
    do {
	if ($self->{VECTOR}->{$i}) {
	    undef $self->{VECTOR}->{$i};
	    $carry = 1;
	}
	else {
	    $self->{VECTOR}->{$i} = 1;
	    $carry = 0;
	}
	$i++;
    } until (($carry == 0) || ($i == $self->{N}));
}


# fills the set
sub Fill {
    my($self) = @_;

    my $i;

    for ($i = 0; $i < $self->{N}; $i++) {
	$self->{VECTOR}->{$i} = 1;
    }
}


# clears all bits, empties the set
sub Empty {
    my($self) = @_;
    
    undef $self->{VECTOR};
}


# returns a bitstring
sub to_Bin {
    my($self) = @_;

    my($i, $s);

    $s = "";
    for ($i = 0; $i < $self->{N}; $i++) {
	if ($self->{VECTOR}->{$i}) {
	    $s = "1" . $s;
	}
	else {
	    $s = "0" . $s;
	}
    }
    return($s);
}


# returns a comma-separated list of features that are on
sub to_Enum {
    my($self) = @_;

    return(join(',', keys(%{$self->{VECTOR}})));
}


# expects a comma-separated list of numbers
sub new_Enum {
    my($this, $n, $s) = @_;

    my($self);

    $self = $this->new($n);
    for (split(/,/,$s)) {
	$self->{VECTOR}->{$_} = 1
    }
    return($self);
}


# returns the length of the vector
sub Size {
    my($self) = @_;

    return($self->{N});
}


# $len1 >= 0, $len2 >= 0, $vec2->Size() >= $off2 >= 0
# $vec1->Size() >= $off1 >= 0
# $vec1 may be undefined
# $len1 + $off1 < $vec->Size()
sub Interval_Substitute {
    my($vec2, $vec1, $off2, $len2, $off1, $len1) = @_;

    my($i,
       $oldvec2);

    # save  $vec2 in $oldvec2
    $oldvec2 = $vec2->Clone();
    # determine the new length for $vec2
    if ($off2 == $vec2->Size()) {
	# we are appending bits from source $vec1
	$vec2->{N} += $len1;
    }
    else {
	# we are inserting bits from $vec1
	$vec2->{N} += $len1 - $len2;
    }
    # target $vec2 changes only from $off2
    # copy the new bits from the source $vec1
    if (defined($vec1)) { # we have source bits
	for ($i = $off2; $i < $off2 + $len1; $i++) {
	    $vec2->{VECTOR}->{$i} = $vec1->{VECTOR}->{$i - $off2 + $off1};
	}
    }
    # append the rest of $oldvec2
    # index $i runs for $oldvec2, we correct for $vec2
    for ($i = $off2 + $len2; $i < $oldvec2->{N}; $i++) {
	$vec2->{VECTOR}->{$i-$len2+$len1} = $oldvec2->{VECTOR}->{$i};
    }
}


# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__

# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Statistics::SparseVector - Perl extension for manipulating sparse bitvectors

=head1 SYNOPSIS

 use Statistics::SparseVector;

 # methods that create new bitvectors
 $vec = Statistics::SparseVector->new($n);
 $vec2 = $vec1->Clone();
 $vec = Statistics::SparseVector->new_Enum($s);
 $vec = Statistics::SparseVector->new_Bin($s);

 # miscellaneous
 $vec2->Substitute_Vector($vec1, $of22, $len2, $off1, $len1);
 $vec->Size();
 $vec->to_Enum();
 $vec->to_Bin();
 $vec->Fill();
 $vec->Empty();
 $vec->increment();

 # manipulation on the bit level
 $vec->Bit_Off($i);
 $vec->Bit_On($i);
 $vec->bit_flip($i);
 $vec->bit_test($i);

 # overloaded operators
 # increment
 $vec++;
 # stringify
 "$vec"


=head1 DESCRIPTION

This module implements sparse bitvectors. Several methods for
manipulating bitvectors are implemented. On purpose the naming of the
methods is identical to the Bit::Vector package by Stephen Beyer; if
you find your vectors to be sparse (have little bits that are on) you
can easily switch to a less memory consuming representation.


=head2 Creation of bitvectors

=over 4

=item C<new>

=item C<Clone>

=item C<new_Enum>

=item C<new_Bin>

=back


=head2 Vector-wide manipulation of bits

=over 4

=item C<Substitute_Vector>

=item C<Fill>

=item C<Empty>

=item C<increment>

=item C<Bit_Off>

=item C<Bit_On>

=item C<bit_flip>

=item C<bit_test>

=back


=head2 Miscellany

=over 4

=item C<Size>

=item C<to_Enum>

=item C<to_Enum>

=item C<to_Bin>

=back


=head2 Overloaded operators

=over 4

=item C<++>

=item Double quotes

=back 


=head2 Remarks about the implementation

=over 4

=item *

Internally sparse vectors are represented by hashes.

=item *

Only a few methods from Bit::Vector are implemented. Maybe new ones
will follow in the future.

=item *

Method C<Substitute_Vector> is not thorougly debugged.

=back


=head1 VERSION

Version 0.1.


=head1 AUTHOR

=begin roff

Hugo WL ter Doest, terdoest@cs.utwente.nl

=end roff

=begin latex

Hugo WL ter Doest, \texttt{terdoest\symbol{'100}cs.utwente.nl}

=end latex


=head1 SEE ALSO

L<perl(1)>, L<Bit::Vector(3)>, L<Statistics::MaxEntropy(3)>, 
L<Statistics::ME.wrapper.pl(3)>, L<Statistics::Candidates(3)>.


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
