
### Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast

package Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast;

use strict;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Blast);

sub split_gapped_alignments {
    my ( $self, $flag ) = @_;

    if ( defined $flag ) {
        $self->{'_split_gapped_alignments'} = $flag ? 1 : 0;
    }
    return $self->{'_split_gapped_alignments'} || 0;
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

    my @parsers;
    if ( defined($fh) ) {
        @parsers = ( Bio::Tools::BPlite->new( -fh => $fh ) );
    }
    else {
        @parsers = $self->get_parsers;
    }

    my $thresh_type  = $self->threshold_type;
    my $thresh       = $self->threshold;
    my $query_length = $self->clone->length or $self->throw("Couldn't get query length");
    warn "query length = $query_length\n";

    my $best_hits = {};
    foreach my $parser (@parsers) {
        while ( my $sbjct = $parser->nextSbjct ) {

            # Get the right identifier from the Subject name
            my $name = $sbjct->name;
                        
            if ( $name =~ /^[^\|]+\|([^\|]+)\|/ ) {
                $name = $1;
            }
            elsif ( $name =~ /^\S+\s+(\S+\.\d+)\b/ ) {
                $name = $1;
            }
            elsif ( $name =~ /^\S+\s+(\S+)/ ) {
                $name = $1;
            }
            elsif ( $name =~ /^>(\S+)$/ ) {
                $name = $1;
            }
            else {
                $self->throw("Can't parse name '$name'");
            }

            #warn "subject name = '$name'\n";

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
                    my $pc = $h->percent;
                    if ( $pc >= $thresh ) {
                        $above_thresh = 1;
                    }
                    if ( $pc > $best_value ) {
                        $best_value = $pc;
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

    my $split_flag = $self->split_gapped_alignments;

    ### Set max coveage -- should be a config parameter
    my $max_coverage = 10;
    $self->throw("Max coverage '$max_coverage' is beyond limit of method '255'") if $max_coverage > 255;

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
                @hsps = $self->_discard_worst_overlapping(@hsps);
                 
                    
                #print STDERR "    Looking at $name ";
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
                
                # Make FeaturePairs if we want to keep this hit
                if ($keep_hit) {
                    if ($split_flag) {

                        # Split gapped HSPs into ungapped feature pairs
                        foreach my $hsp (@hsps) {
                            
                            $self->split_HSP( $hsp, $name, $ana );
                        }
                    }
                    else {

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

1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

