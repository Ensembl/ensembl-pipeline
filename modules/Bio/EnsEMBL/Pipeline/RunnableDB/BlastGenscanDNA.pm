#!/usr/local/bin/perl -w

#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $genscan     = Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanDNA->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$genscan->fetch_input();
$genscan->run();
$genscan->output();
$genscan->write_output(); #writes to DB

=head1 DESCRIPTION

his object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides constructed from 
assembling genscan predicted features into contiguous DNA sequence. 
The resulting blast hits are written back as FeaturePairs.
The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanDNA;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::RunnableDB::Blast;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep);

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;
    #need to pass one peptide at a time
    $self->throw("Input must be fetched before run") unless ($self->genseq);
    foreach my $transcript ($self->transcripts)
    {
        print STDERR "Creating BioPrimarySeq ".$transcript->id."\n";
        my $gene_seq    = Bio::PrimarySeq->new(
                        -id         => $transcript->id,
                        -seq        => $transcript->dna_seq->seq,
                        -moltype    => 'protein' );
                                        
        my %parameters = $self->formatted_parameters();
        my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(%parameters);
        $runnable->clone($gene_seq);
        $runnable->threshold($self->threshold());
        $runnable->run();
        $self->runnable($runnable);                                        
    }

    my (@output);
    my @runnable = $self->runnable();
    #align the output to the contig
    foreach my $run (@runnable)
    {
        foreach my $transcript ($self->transcripts)
        {
            if ($run->clone->id eq $transcript->id)
            {
                print STDERR "MATCHED: ".$run->clone->id." with ".$transcript->id."\n";
                $self->align_hit_to_contig($run, $transcript);
                last;
            }
        }
    }
}

sub align_hit_to_contig {
    my ($self, $run, $trans) = @_;
    my @exon_aligns;
    
    #map each exon onto a single dna sequence
    foreach my $exon ($trans->each_Exon)
    {
        my %ex_align;
        my $mrna    = $trans->dna_seq->seq 
            or $self->throw("Unable to read sequence from Transcript ".$trans->name."\n");
        my $ex_seq  = $exon->seq->seq
            or $self->throw("Unable to read sequence from Exon ".$exon->name."\n");
        
        $ex_align {'name'}      = $exon->id;
        $ex_align {'gen_start'} = $exon->start;
        $ex_align {'gen_end'}   = $exon->end;
        $ex_align {'strand'}    = $exon->strand;
        $ex_align {'phase'}     = $exon->phase;
        $ex_align {'end_phase'} = $exon->end_phase;
        $ex_align {'cdna_start'}= index($mrna, $ex_seq) +1;
        $ex_align {'cdna_end'}  = ($ex_align {'cdna_start'} + length($ex_seq)) -1;   
        push (@exon_aligns, \%ex_align);
    }
    
    #map each feature onto 1 or more exons
    foreach my $fp ($run->output)
    { 
        if ($fp->strand != 1 && $fp->strand != -1)
        {
            warn("Strand not set to 1 or -1 for Feature "
                   .$fp->seqname."(".$fp->strand.")\n");
            next;
        }
        
        my @aligned_exons;
        foreach my $ex (@exon_aligns)
        {
            if ($ex->{'strand'} != 1 && $ex->{'strand'} != -1)
            {
                warn("Strand not set to 1 or -1 for Exon "
                       .$ex->{'name'}."(".$ex->{'strand'}.")\n");
                next;
            }            
            
            if ($ex->{'cdna_end'} >= $fp->start)
            {
                push (@aligned_exons, $ex);
                last if ($ex->{'cdna_end'} >= $fp->end);
            }
        }
        $self->create_dna_featurepairs($fp, @aligned_exons);
    }
}

sub create_dna_featurepairs {
my  ($self, $fp, @aligned_exons) = @_;

    #split features across exons if necessary
    foreach my $ex (@aligned_exons)
    {
        #my ($start, $end, $hstart, $hend);
        my ($start_phase, $end_phase, $start, $end, $hstart, $hend);
    
        if ($$ex{'cdna_start'} < $fp->start)
        {
            $start  =   $ex->{'gen_start'} + ($fp->start - $ex->{'cdna_start'});
            $hstart =   $fp->hstart; 
            $start_phase = 0; 
        }
        else
        {
            $start  =   $ex->{'gen_start'};
            $hstart =   ($ex->{'cdna_start'} - $fp->start) + $fp->hstart;
            $start_phase = $ex->{'start_phase'};
        }
        
        if ($$ex{'cdna_end'} > $fp->end)
        {
            $end    =   $ex->{'gen_start'} + ($fp->end - $ex->{'cdna_start'});
            $hend   =   $fp->hend;
            $end_phase = 0;
        }
        else
        {
            $end    =   $ex->{'gen_end'};
            $hend   =   $fp->hstart + ($ex->{'cdna_end'} - $fp->start); 
            $end_phase = $ex->{'end_phase'};
        }
        
        #print STDERR "NAME: ".$fp->hseqname
        #             ."\tEx: ".($ex->{'cdna_end'} - $ex->{'cdna_start'})
        #             ."\tF1: ".($fp->end - $fp->start +1)
        #             ."\tf2: ".($fp->hend - $fp->hstart +1)
        #             ."\th1: ".($end - $start + 1)
        #             ."\th2: ".($hend - $hstart +1)."\n";
        
        #print STDERR "Ex ".$ex->{'cdna_start'}." - ".$ex->{'cdna_end'}
        #             ." Feat ".$fp->start." - ".$fp->end
        #             ." Feat2 ".$fp->hstart." - ".$fp->hend
        #             ." H1 $start - $end H2 $hstart - $hend"
        #             ." Ex ".$ex->{'gen_start'}." - ".$ex->{'gen_end'}."\n";
                     
        my $contig_feat = Bio::EnsEMBL::SeqFeature->new (
                                -seqname    =>  $ex->{'name'},
                                -start      =>  $start, 
                                -end        =>  $end,
                                -strand     =>  $fp->strand,
                                -score      =>  $fp->score,
                                -p_value    =>  $fp->p_value,
                                -percent_id =>  $fp->percent_id,
                                -analysis   =>  $fp->analysis,
                                -primary_tag=>  $fp->primary_tag,
                                -source_tag =>  $fp->source_tag, 
                                -phase      =>  $start_phase,  
                                -end_phase  =>  $end_phase );
                                
        my $hit_feat = Bio::EnsEMBL::SeqFeature->new (
                                -seqname    =>  $fp->hseqname,
                                -start      =>  $hstart, 
                                -end        =>  $hend,
                                -strand     =>  $fp->hstrand,
                                -score      =>  $fp->score,
                                -p_value    =>  $fp->p_value,
                                -percent_id =>  $fp->percent_id,
                                -analysis   =>  $fp->analysis,
                                -primary_tag=>  $fp->primary_tag,
                                -source_tag =>  $fp->source_tag, 
                                -phase      =>  $start_phase,  
                                -end_phase  =>  $end_phase );
        
        my $featurepair = Bio::EnsEMBL::FeaturePair->new (
                                -feature1   => $contig_feat,
                                -feature2   => $hit_feat );
        
        $self->featurepairs($featurepair);
    }
}

