# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..2\n"; }
END {print "not ok 1\n" unless $loaded;}
use Statistics::MaxEntropy qw(GIS 
		      IIS 
		      FI 
		      init 
		      random_init 
		      done
		      $debug
		      $normalise
		      $NEWTON_max_it
		      $NEWTON_min
		      $KL_max_it
		      $KL_min
		      $rand_max_candidates
		      $rand_max_features);
use vars qw($TMP
	    $events_file
            $parameters_file
            $candidates_file
	    $new_events_file
            $new_parameters_file
            $new_candidates_file
            $info_file_1
	    $info_file_2
	    );
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):


# debugging messages; default 0
$debug = 0;
# normalise frequencies in events file; default 1
$normalise = 1;
# maximum number of iterations for IIS; default 100
$NEWTON_max_it = 100;
# minimal distance between new and old x for Newton's method; default 0.001
$NEWTON_min = 0.0001;
# maximum number of iterations for Newton's method; default 100
$KL_max_it = 100;
# minimal distance between new and old x; default 0.001
$KL_min = 0.001;

# configuration of Statistics::MaxEntropy::random_init
# maximum number of candidates randomly generated
$rand_max_candidates = 10;
# maximum number of features randomly generated
$rand_max_features = 20;

$TMP = "/tmp";
$events_file = "$TMP/events.txt";
$parameters_file = "$TMP/parameters.txt";
$candidates_file = "$TMP/candidates.txt";
$new_events_file = "$TMP/new.events.txt";
$new_parameters_file = "$TMP/new.parameters.txt";
$new_candidates_file = "$TMP/new.candidates.txt";
$info_file_1 = "$TMP/information.1.txt";
$info_file_2 = "$TMP/information.2.txt";

# test 1, with random events
random_init($events_file, $parameters_file, $candidates_file);
FI(2, \&IIS);
done($new_events_file, $new_parameters_file, $new_candidates_file,
     $info_file_1);
print "ok 1\n";


# test 2, events are in file
init($events_file, $parameters_file, $candidates_file);
FI(1, \&GIS);
done($new_events_file, $new_parameters_file, $new_candidates_file,
     $info_file_2);

print "ok 2\n";


