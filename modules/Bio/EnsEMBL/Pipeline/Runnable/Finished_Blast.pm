### Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast

package Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast;

use strict;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Symbol;
use Bio::EnsEMBL::Pipeline::Tools::BPlite;
use Bio::EnsEMBL::Pipeline::Config::Blast;
use BlastableVersion;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Blast);

$ENV{BLASTDB} = '/data/blastdb/Ensembl';

BEGIN {
    print "\nUSING " . __PACKAGE__ . "\n\n";
}

my %FASTA_HEADER;
my %BLAST_FLAVOUR;
my %REFILTER;
my %MAX_COVERAGE;
my %DISCARD_OVERLAPS;
my %SPLIT_GAPPED_ALIGNMENTS;
my %THRESHOLD_TYPE;
my %THRESHOLD;

foreach my $db (@$DB_CONFIG) {
    my ( $name,
	 $header,
	 $flavour,
	 $refilter,
	 $max_coverage,
	 $discard_overlaps,
	 $split_gapped_alignments,
	 $threshold_type,
	 $threshold ) = ( $db->{'name'}, 
			  $db->{'header'}, 
			  $db->{'flavour'}, 
			  $db->{'refilter'}, 
			  $db->{'max_coverage'}, 
			  $db->{'discard_overlaps'},
			  $db->{'split_gapped_alignments'},
			  $db->{'threshold_type'},
			  $db->{'threshold'});    
    if($db && $name){
      $FASTA_HEADER{$name} = $header; 
      $BLAST_FLAVOUR{$name} = $flavour; 
      $REFILTER{$name} = $refilter;
      $MAX_COVERAGE{$name} = $max_coverage;
      $DISCARD_OVERLAPS{$name} = $discard_overlaps;
      $SPLIT_GAPPED_ALIGNMENTS{$name} = $split_gapped_alignments;
      $THRESHOLD_TYPE{$name} = $threshold_type;
      $THRESHOLD{$name} = $threshold;
    }else{
      my($p, $f, $l) = caller;
      warn("Either DB: ".$db." or NAME: ".$name." isn't defined so can't work $f:$l\n");
    }
}
sub get_analysis {
    my ($self) = @_;

    my ($ana);
    unless ( $ana = $self->{'_analysis'} ) {
        my ($source) = $self->program =~ m{([^/]+)$}
          or $self->throw( "Can't parse last element from path: '" . $self->program . "'" );
        $ana = $self->{'_analysis'} = Bio::EnsEMBL::Analysis->new(
            -db              => $self->database,
            -db_version      => 1,                 # ARUUGA!!!
            -program         => $source,
            -program_version => 1,
            -gff_source      => $source,
            -gff_feature     => 'similarity',
            -logic_name      => 'blast',
        );
    }
    return $ana;
}
sub parse_results {
    my ( $self, $fh ) = @_;
    print "*** PARSING OUTPUT ***\n";

    my @parsers;

    if ( defined($fh) ) {
        @parsers = ( Bio::Tools::BPlite->new( -fh => $fh ) );
    }
    else {
        @parsers = $self->get_parsers;
    }
    my $db = $self->database;

    my $thresh_type  = defined($THRESHOLD_TYPE{$db}) ? $THRESHOLD_TYPE{$db} : $self->threshold_type;
    my $thresh       = defined($THRESHOLD{$db})      ? $THRESHOLD{$db}      : $self->threshold;

    print "SETTING [parse_results]: threshold type is <" . $thresh_type . ">\n";
    print "SETTING [parse_results]: threshold is <" . $thresh . ">\n";

    my $query_length = $self->query->length or $self->throw("Couldn't get query length");
    print "Query Length: <" . $query_length . ">\n";
    my $best_hits = {};

    
    my $re = $self->get_regex($db);
    if(!$re){
	$self->throw("no regex defined for ".$db);
    }
    print "SETTING [parse_results]: regex is </" . $re . "/>\n";
    foreach my $parser (@parsers) {
        while ( my $sbjct = $parser->nextSbjct ) {
	    my $fasta_header = $sbjct->name ;    
	    my ($name) = $fasta_header =~ /$re/;
	    unless ($name) {
		$self->throw("ERROR - Check BLAST config and/or blast headers. Couldn't get valid accession from '" . $fasta_header . "'");
	    }
#            warn "subject name = '$name'\n";
            my $hsps = [];
            while ( my $h = $sbjct->nextHSP ) {
                push ( @$hsps, $h );
            }
            # Is there an HSP in this subject above our threhold?
            my $above_thresh = 0;
            my ($best_value);
            if ( $thresh_type eq 'PID' ) {
                $best_value = 0;
                foreach my $h (@$hsps) {
                    my $percent_id = $h->percent;
                    if ( $percent_id >= $thresh ) {
                        $above_thresh = 1;
                    }
                    if ( $percent_id > $best_value ) {
                        $best_value = $percent_id;
                    }
                }
            }
            elsif ( $thresh_type eq 'PVALUE' ) {
                $best_value = 1e6;    # That's a pretty bad PVALUE!
                foreach my $hsp (@$hsps) {
                    my $p_val = $hsp->P;
                    if ( $p_val <= $thresh ) {
                        $above_thresh = 1;
                    }
                    if ( $p_val < $best_value ) {
                        $best_value = $p_val;
                    }
                }
            }
            else {
                $self->throw("Unknown threshold type '$thresh_type'");
            }
            next unless $above_thresh;

            my $top_score = 0;
            foreach my $hsp (@$hsps) {
                my $score = $hsp->score;
                if ( $score > $top_score ) {
                    $top_score = $score;
                }
            }
            my $best = $best_hits->{$best_value}{$top_score} ||= [];

            # Put as first element of hsps array
            unshift ( @$hsps, $name );

            push ( @$best, $hsps );
        }
    }
    $self->_apply_coverage_filter( $query_length, $best_hits );
}

sub _apply_coverage_filter {
    my ( $self, $query_length, $best_hits ) = @_;

    my %id_list;

    my $split_flag = defined($SPLIT_GAPPED_ALIGNMENTS{$self->database}) ? $SPLIT_GAPPED_ALIGNMENTS{$self->database} : 1;
    my $max_coverage = ( defined($MAX_COVERAGE{$self->database}) ? $MAX_COVERAGE{$self->database} : 0 ) % 256;
    my $discard_overlaps = $DISCARD_OVERLAPS{$self->database} || 0;
    
    print "SETTING [_apply_coverage_filter]: split_flag is <" . ( $split_flag ? "ON" : "OFF" ) . ">\n";
    print "SETTING [_apply_coverage_filter]: max_coverage is <" . $max_coverage . ">\n";
    $self->throw("Max coverage '$max_coverage' is beyond limit of method '255'") if $max_coverage > 255;
    print "SETTING [_apply_coverage_filter]: discard_overlaps is <" . ( $discard_overlaps ? "ON" : "OFF" ) . ">\n";

    # Make a string of nulls (zeroes) the length of the query
    my $coverage_map = "\0" x ( $query_length + 1 );

    my $ana = $self->get_analysis;

    # Loop through from best to worst according
    # to our threshold type.
    my (@bin_numbers);
    my $thresh_type = $self->threshold_type;
    if ( $thresh_type eq 'PID' ) {
	@bin_numbers = sort { $b <=> $a } keys %$best_hits;
    }
    elsif ( $thresh_type eq 'PVALUE' ) {
	@bin_numbers = sort { $a <=> $b } keys %$best_hits;
    }
    else {
	$self->throw("Unknown threshold type '$thresh_type'");
    }

    foreach my $bin_n (@bin_numbers) {
	       
	#print STDERR "\nLooking at bin: $thresh_type $bin_n\n";
	my $score_hits = $best_hits->{$bin_n};

	# Loop through from best to worst according
	# to score within threshold bin
	foreach my $top_score ( sort { $b <=> $a } keys %$score_hits ) {

	    #print STDERR "  Looking at hits scoring $top_score\n";
	    my $hits = $score_hits->{$top_score};

	    # Loop through each hit with this score
	    foreach my $hit (@$hits) {
                
		my $keep_hit = 0;
                
		my ( $name, @hsps ) = @$hit;

		# Don't keep multiple matches to the same sequence
		# at the same genomic location.
		@hsps = $self->_discard_worst_overlapping(@hsps) if $discard_overlaps;

		unless ($max_coverage == 0){
#		print STDERR "    Looking at $name ";
		    foreach my $hsp (@hsps) {
			my $q = $hsp->query;
			foreach my $i ( $q->start .. $q->end ) {
			    my $depth = unpack( 'C', substr( $coverage_map, $i, 1 ) );
			    if ( $depth < $max_coverage ) {
				$keep_hit = 1;
				$depth++;
				# Increment depth at this position in the map
				substr( $coverage_map, $i, 1 ) = pack( 'C', $depth );
			    }
			}
		    }
		}
		# Make FeaturePairs if we want to keep this hit
		if ($keep_hit || $max_coverage == 0) {
		    if ($split_flag) {
			# Split gapped HSPs into ungapped feature pairs
			foreach my $hsp (@hsps) {                            
			    $self->split_HSP( $hsp, $name, $ana );
			}
		    } else {
                        # Don't bother splitting gapped alignments
			my (@feat);
			foreach my $hsp (@hsps) {
			    my $q  = $hsp->query;
			    my $s  = $hsp->subject;
			    my $fp =
				$self->_makeFeaturePair( $q->start, $q->end, $q->strand, $s->start, $s->end, $s->strand,
							 $hsp->score, $hsp->percent, $hsp->P, $name, $ana, );
			    $self->growfplist($fp);
			}
		    }
		}
	    }
	}
    }
    print "*** FINISHED PARSING ***\n";
 }


sub _discard_worst_overlapping {
    my $self = shift;
    my @hsps = sort { $a->query->start <=> $b->query->start } @_;

    # Put all the hsps hits into overlapping bins
    my $first   = shift @hsps;
    my $start   = $first->start;
    my $end     = $first->end;
    my $current = [$first];        # The current bin
    my @bins    = ($current);
    while ( my $hsp = shift @hsps ) {
        my $q = $hsp->query;

        # Does this hsp overlap the current bin?
        if ( $q->start <= $end and $q->end >= $start ) {
            push ( @$current, $hsp );
            if ( $q->end > $end ) {
                $end = $q->end;
            }
        }
        else {
            $current = [$hsp];
            push ( @bins, $current );
            $start = $hsp->start;
            $end   = $hsp->end;
        }
    }

    foreach my $bin (@bins) {
        if ( @$bin == 1 ) {
            push ( @hsps, $bin->[0] );
        }
        else {

            # Remove the hsp with the lowest percent identity
            # (or lowest score or highest P value)
            my @bin_hsps = sort { $a->percent <=> $b->percent || $a->score <=> $b->score || $b->P <=> $a->P } @$bin;
            shift (@bin_hsps);

            if ( @bin_hsps == 1 ) {

                # We are left with 1 hsp, so add it.
                push ( @hsps, $bin_hsps[0] );
            }
            else {

                # Remaining hsps may not be overlapping
                push ( @hsps, $self->_discard_worst_overlapping(@bin_hsps) );
            }
        }
    }

    return @hsps;
}

sub fetch_databases {
    my ($self) = @_;
    
    my @databases;

    my $dbname = $self->database; 
    #print "fetching databases for ".$dbname."\n";
    $dbname =~ s/\s//g;
    foreach my $dbname(split(",", $dbname)){ # allows the use of a comma separated list in $self->database
	# prepend the environment variable $BLASTDB if
	# database name is not an absoloute path

	unless ($dbname =~ m!^/!) {
	    $dbname = $ENV{BLASTDB} . "/" . $dbname;
	}

	# If the expanded database name exists put this in
	# the database array.
	#
	# If it doesn't exist then see if $database-1,$database-2 exist
	# and put them in the database array
	if (-f $dbname) {
	    push(@databases,$dbname);
	    $self->db_version_searched($dbname);
	} else {
	    my $count = 1;

	    while (-f $dbname . "-$count") {
		push(@databases,$dbname . "-$count");
		$self->db_version_searched($dbname ."-".$count);
		$count++;
	    }
	    $! = undef; # to stop pollution as it will be "No such file or directory" after while loop above.
	}
    }
    if (scalar(@databases) == 0) {
	$self->throw("No databases exist for " . $dbname);
    }

    return @databases;

}
sub db_version_searched{
    my $self = shift;
    my $db = shift;
    unless($self->{'_db_version_searched'} && !$db){
    	my $ver = BlastableVersion->new();
    	my $dbv = $ver->get_version($db);
    	$self->{'_db_version_searched'} = $dbv;
    }
    return $self->{'_db_version_searched'};
}
1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast

=head1 DESCRIPTION

Unlike Bio::EnsEMBL::Pipeline::Runnable::Blast,
this module creates FeaturePairs from HSPs after
doing any depth filtering to save time and memory
when searching genomic sequences that generate
large numbers of blast matches.

=head2 usage of Bio::EnsEMBL::Pipeline::Config::Blast with this module


DB_CONFIG => [
        { 

	    name    => 'embl_vertrna,emnew_vertrna', # use test of field [analysis.db]

            type    => 'dna',

            header  => '\S+\s+(\S+)\s',

	    flavour => 'wu',

	    ungapped => 0,

	    discard_overlaps => 0,# boolean(run/not run) _discard_worst_overlapping

	    max_coverage => 254,  # set max coverage filter depth, max = 255

            split_gapped_alignments => 1, # boolean to split_gapped_alignments

            threshold => '1e4',

            threshold_type => 'PVALUE',
        },
        {
            name    => '',
            ...
        },
        ...
]

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

