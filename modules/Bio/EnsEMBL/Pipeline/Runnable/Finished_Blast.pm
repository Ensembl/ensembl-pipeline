
### Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast

package Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast;

use strict;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use vars '@ISA';
@ISA = qw{ Bio::EnsEMBL::Pipeline::Runnable::Blast };

sub get_analysis {
    my( $self ) = @_;
    
    my( $ana );
    unless ($ana = $self->{'_analysis'}) {
        my ($source) = $self->program =~ m{([^/]+)$}
            or $self->throw("Can't parse last element from path: '". $self->program ."'");
        $ana = $self->{'_analysis'} = Bio::EnsEMBL::Analysis->new(
            -db              => $self->database,
            -db_version      => 1,  # ARUUGA!!!
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
    my( $self, $fh ) = @_;

    my @parsers;
    if (defined($fh)) {
        @parsers = (Bio::Tools::BPlite->(-fh => $fh));
    } else {
        @parsers = $self->get_parsers;
    }

    my $thresh_type = $self->threshold_type;
    my $thresh      = $self->threshold;
    foreach my $parser (@parsers) {
        while (my $sbjct = $parser->nextSbjct) {
            my( @hsps );
	    while (my $hsp = $sbjct->nextHSP) {
	        push(@hsps, $hsp);
	    }
            
            # Is there an HSP in this subject above our threhold?
            my $above_thresh = 0;
            if ($thresh_type eq 'PID') {
                foreach my $h (@hsps) {
                    if ($h->percent >= $thresh) {
                        $above_thresh = 1;
                        last;
                    }
                }
            }
            elsif ($thresh_type eq 'PVALUE') {
                foreach my $h (@hsps) {
                    if ($h->P <= $thresh) {
                        $above_thresh = 1;
                        last;
                    }
                }
            }
            else {
                $self->throw("Unknown threshold type '$thresh_type'");
            }
            next unless $above_thresh;
            
            # Get the right identifier from the Subject name
            my $name = $sbjct->name;
            if ($name =~ /^[^\|]+\|([^\|]+)\|/) {
                $name = $1;
            }
            elsif ($name =~ /^(\S+)/) {
                $name = $1;
            }
            else {
                $self->throw("Can't parse name '$name'");
            }
            
            foreach my $hsp (@hsps) {
                my $q = $hsp->query;
                my $s = $hsp->subject;
                my $fp = $self->_makeFeaturePair(
                    $q->start, $q->end, $q->strand,
                    $s->start, $s->end, $s->strand,
                    $hsp->score, $hsp->percent, $hsp->P,
                    $name, $self->get_analysis,
                    );
                
                # growfplist is in RunnableI, but the output method in
                # Runnable::Blast accesses the array directly!  Bad!
                $self->growfplist($fp);
            }
        }
    }
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

