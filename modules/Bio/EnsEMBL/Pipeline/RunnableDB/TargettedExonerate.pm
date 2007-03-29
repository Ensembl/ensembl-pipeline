package  Bio::EnsEMBL::Pipeline::RunnableDB::TargettedExonerate;

use strict;

use vars qw(@ISA);
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Sequences qw (
							     GB_PROTEIN_INDEX
							     GB_PROTEIN_SEQFETCHER
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							     GB_GW_DBNAME
                                                             GB_GW_DBPORT                                          
							    );

use Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGenewise;
use Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes;



@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGenewise);



sub fetch_input{
  my ($self) = @_;
  my $runnable = Bio::EnsEMBL::Analysis::RunnableDB::Exonerate2Genes->new
    (
     -input_id => $self->input_id,
     -db => $self->db,
     -analysis => $self->analysis,
    );
  $self->runnable($runnable);
}


sub run{
  my ($self) = @_;
  eval{
    $self->runnable->fetch_input();
    $self->runnable->run();
  };
  if($@){
    throw("Failed in run of ".$self->runnable." $@");
  }
  $self->validate_genes();
}


sub validate_genes{
  my ($self) = @_;

  my @genes;

  my @seqfetchers;
  push (@seqfetchers, $self->seqfetcher);
  
  #print "VALIDATING GENES\n";

  GENE:foreach my $gene(@{$self->runnable->output}){
      foreach my $transcript(@{$gene->get_all_Transcripts}){

        #print "TRANSCRIPT: ",$transcript,"\tStart: ",$transcript->start,"\tEnd: ",$transcript->end,"\n";
 
        my $valid_transcripts = 
          Bio::EnsEMBL::Pipeline::Tools::GeneUtils->validate_Transcript
              ($transcript,
               $self->query,
               $GB_TARGETTED_MULTI_EXON_COVERAGE,
               $GB_TARGETTED_SINGLE_EXON_COVERAGE,
               $GB_TARGETTED_MAX_INTRON,
               $GB_TARGETTED_MIN_SPLIT_COVERAGE,
               \@seqfetchers
              );

        next GENE unless defined $valid_transcripts;
        foreach my $checked_transcript (@$valid_transcripts){

          # add a start codon if appropriate
          $checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_start_codon($checked_transcript);
          
          # add a stop codon if appropriate
          $checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon($checked_transcript);
          
          # attach analysis to supporting features
          foreach my $exon(@{$checked_transcript->get_all_Exons}){
            foreach my $sf(@{$exon->get_all_supporting_features}){
              $sf->analysis($self->analysis);
            }
          }
          my $gene   = new Bio::EnsEMBL::Gene;
          $gene->type($self->analysis->logic_name);
          $gene->analysis($self->analysis);
          $gene->add_Transcript($checked_transcript);
          push(@genes,$gene);
        }
      }
    }
  $self->output(@genes);
}
