#! /Utils/bin/perl5
#! /usr/bin/perl

##---------------------------------------------------------------------------##
##  Author: Hugo WL ter Doest       terdoest@cs.utwente.nl
##  Description: Wrapper around Statistics::ME
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
use Getopt::Long;
use Statistics::ME;


##---------------------------------------------------------------------------##
##	Globals
##---------------------------------------------------------------------------##
#($main::PROG = $0) =~ s/.*\///;
#$main::VERSION = "0.4";




##---------------------------------------------------------------------------##
##	Routines
##---------------------------------------------------------------------------##

# parse command line and read some files
sub prologue {
    
    GetOptions("debug!" , \$Statistics::ME::debug,
	       "i_events=s" , \$main::events,
	       "i_parameters=s" , \$main::init_pars,
	       "i_candidates=s" , \$main::cands,
	       "o_events=s" , \$main::new_events,
	       "o_candidates=s" , \$main::new_cands,
	       "o_parameters=s" , \$main::new_pars,
	       "o_information=s" , \$main::information,
	       "KL_max_it=i", \$Statistics::ME::KL_max_it,
	       "NEWTON_max_it=i", \$Statistics::ME::NEWTON_max_it,
	       "KL_min=f", \$Statistics::ME::KL_min,
	       "NEWTON_min=f", \$Statistics::ME::NEWTON_min,
	       "normalise!", \$Statistics::ME::normalise,
	       "nr_to_add=i", \$main::nr_to_add
	       );

    Statistics::ME::init($main::events, 
			 $main::init_pars,
			 $main::cands);
}


sub run {
    if ($main::nr_to_add &&
	($main::cands)) {
      Statistics::ME::FI($main::nr_to_add);
    }
    else {
      Statistics::ME::IIS();
    }
}


sub epilogue {
  Statistics::ME::done($main::new_events,
		       $main::new_pars,
		       $main::new_cands,
		       $main::information);
}


##---------------------------------------------------------------------------##
				##------------##
				## Begin MAIN ##
				##------------##
prologue();
run();
epilogue();
				##------------##
				##  End MAIN  ##
				##------------##
##---------------------------------------------------------------------------##

