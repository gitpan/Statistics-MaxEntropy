#!/Utils/bin/perl5
#!/usr/bin/perl

print "1..7\n";
$i = 1;

use Statistics::MaxEntropy qw($debug
			      $NEWTON_max_it
			      $NEWTON_min
			      $KL_max_it
			      $KL_min
			      $SAMPLE_size);

use Statistics::Candidates;

use vars qw($scaling
	    $sampling
	    $i
	    $TMP
	    $events_file
            $parameters_file
            $candidates_file
	    $new_events_file
            $new_candidates_file
            $dump_file_1
            $dump_file_2);

# debugging messages; default 0
$debug = 0;
# maximum number of iterations for IIS; default 100
$NEWTON_max_it = 25;
# minimal distance between new and old x for Newton's method; default 0.001
$NEWTON_min = 0.0001;
# maximum number of iterations for Newton's method; default 100
$KL_max_it = 100;
# minimal distance between new and old x; default 0.0001
$KL_min = 0.00001;

$TMP = "/tmp";
$events_file = "data/events.txt";
$parameters_file = "data/parameters.txt";
$candidates_file = "data/candidates.txt";
$new_events_file = "$TMP/new.events.txt";
$new_candidates_file = "$TMP/new.candidates.txt";
$dump_file_1 = "$TMP/dump.1.txt";


# test the scalers for each of the sampling methods
$events=Statistics::MaxEntropy->new($events_file);
$SAMPLE_size = 100;
for $sampling ("enum", "corpus", "mc") {
    for $scaling ("iis", "gis") {
 	$events->clear();
 	$events->scale($sampling, $scaling);
 	print "ok $i\n";
	$i++;
    }
}

# dump the event space
$events->dump($dump_file_1);
print "ok $i\n";

__END__
