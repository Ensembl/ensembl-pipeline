
### Bio::EnsEMBL::Pipeline::Runnable::Finished_MiniEst2Genome

package Bio::EnsEMBL::Pipeline::Runnable::Finished_MiniEst2Genome;

use strict;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome;
# Inherit from MiniEst2Genome
use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use vars '@ISA';
@ISA = ('Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome');

sub run_blaste2g {
    my ( $self, $est, $features ) = @_;
    
    my $miniseq = $self->make_miniseq(@$features);

    my $hseq    = $self->get_Sequence($est) or $self->throw("Can't fetch sequence for id '$est'");

    my $eg = new Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome(
        -genomic => $miniseq->get_cDNA_sequence,
        -est     => $hseq,
    );
    $eg->run;

     foreach my $fp ( @{$eg->output} ) {

         my @converted = $miniseq->convert_FeaturePair($fp);

         if ( @converted > 1 ) {
             warn "feature converts into '" . scalar(@converted) . "' > 1 features - ignoring\n";
         } else {
            # convert_FeaturePair zaps strand and hseqname,0
            # so we put them back here.
            my $new = $converted[0];

#            $new->seqname( $fp->hseqname );  # 
#            $new->strand( $fp->strand );     #
#            $new->hseqname( $fp->hseqname ); #
#            $new->percent_id( $fp->percent_id ); #
#            $new->feature2->percent_id( $fp->percent_id ); #
# reverse logic here, change the $fp passed to equal the converted [miniseq]
# change seqname, hend, hstart, gsf_start, gsf_end
# hend & hstart appear to only be calculated from the length of the end & start
# producing the same length for the hit as the contig.  not good....
#	    $fp->hend(   $new->hend   );
#	    $fp->hstart( $new->hstart );
	    $fp->end(    $new->end   );
	    $fp->start(  $new->start );
#            $self->add_output(@converted);
#	    print STDERR "CIGAR_STRING : " . $fp->cigar_string . "\n";
#	    print STDERR $fp->gffstring . "\n";
            $self->add_output($fp);
        }
     }
}


sub add_output{
    my($self, @feat_pairs) = @_;
    push( @{ $self->{'_output'} } , @feat_pairs);
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_MiniEst2Genome

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

