#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CombinedGeneBuild.pm
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CombinedGeneBuild

=head1 SYNOPSIS

=head1 DESCRIPTION

Runs the various stages of the gene build process on a clone

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_CombinedGeneBuild;

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::Contig_TargettedGeneE2G;
use Bio::EnsEMBL::Pipeline::RunnableDB::Contig_BlastMiniGenewise;
use Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Gene_Builder;
use Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::PmatchFeature;
use Bio::EnsEMBL::Pipeline::SeqFetcher::getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
# config file; parameters searched for here if not passed in as @args
require "Bio/EnsEMBL/Pipeline/GB_conf.pl";

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
			   -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Pipeline::Analysis (optional) 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  # dbobj, input_id, seqfetcher, and analysis objects are all set in
  # in superclass constructor (RunnableDB.pm)

  return $self;
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: checks in GB_conf::seqfetch_conf for $indexname; if it exists, 
           returns a Bio::EnsEMBL::Pipeline::SeqFetcher::getseqs, otherwise 
           returns a Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    : $indexname - string


=cut

sub make_seqfetcher {
  my ( $self, $indexname, $conf_hash ) = @_;
  # need to treat contents of $conf_hash as name of hash
  my $index = $::{$conf_hash}{$indexname};
  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::getseqs(
								  '-db' => \@db,
								 );
  }
  else{
    # default to Pfetch
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  }

  return $seqfetcher;
}

=head2 fetch_input

 Title   : fetch_input
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_input{
  my ($self,@args) = @_;

  my $cloneid     = $self->input_id;
  my $clone       = $self->dbobj->get_Clone($cloneid);
  foreach my $contig  ($clone->get_all_Contigs()){
    $self->make_targetted_runnables($contig);
    $self->make_similarity_runnable($contig);
  }
#  $self->make_riken_runnable;
  $self->make_genebuild_runnable;
}

=head2 make_targetted_runnables

 Title   : make_targetted_runnables
 Usage   :
 Function: find all proteins that pmatched to input contig and set up Contig_TargettedGeneE2G jobs for them. 
 Example :
 Returns : 
 Args    :


=cut

sub make_targetted_runnables {
  my ($self, $contig) = @_;

  # set up seqfetchers
  my $protein_fetcher = $self->make_seqfetcher("protein_index", "targetted_conf");
  my $cdna_fetcher    = $self->make_seqfetcher("cdna_index", "targetted_conf");

  # we need to find all the proteins that pmatch into this region
  # take a note of those that fall across the ends of the vc? and do what, precisely?
  # extend the VC? that will completely screw up the final genebuild. Hmmm.
  # do it, track it & see how many are affected.

  my $pmfa = new Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor( $self->dbobj );
  
  # pmatch_feature table was designed to use chrname, start & end but should 
  # work just as well with contigname, start & end
  my $contigname = $contig->id;
  my $start   = 1;
  my $end     = $contig->length - 1;

  foreach my $feat($pmfa->get_PmatchFeatures_by_chr_start_end($contigname, $start, $end)){
    my $input = $feat->chr_name   . ":" . 
                $feat->start      . "," . 
                $feat->end        . ":" . 
                $feat->protein_id . ":" . 
                $feat->cdna_id;

    print STDERR "TGE input: $input\n";

    my $tgr = new Bio::EnsEMBL::Pipeline::RunnableDB::Contig_TargettedGeneE2G(
				       -dbobj           => $self->dbobj,
				       -input_id        => $input,
				       -seqfetcher      => $protein_fetcher,
				       -cdna_seqfetcher => $cdna_fetcher,
#				       -analysis => $analysis,
									     );
    $self->targetted_runnable($tgr);
  }

}

=head2 targetted_runnable

 Title   : targetted_runnable
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub targetted_runnable{
  my ($self,$arg) = @_;

  if (!defined($self->{'_targetted_runnables'})) {
      $self->{'_targetted_runnables'} = [];
  }
  
  if (defined($arg)) {
      if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	  push(@{$self->{'_targetted_runnables'}},$arg);
      } else {
	  $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
      }
  }
  
  return @{$self->{'_targetted_runnables'}};  
}

=head2 make_similarity_runnable

 Title   : make_similarity_runnable
 Usage   :
 Function: makes one Contig_BlastMiniGenewise runnable per contig
 Example :
 Returns : 
 Args    :


=cut

sub make_similarity_runnable {
  my ($self, $contig) = @_;
 
#  analysis setting up? 
#  my $analysis = $self->dbobj->get_Analysis_Adaptor->??;

  my $seqfetcher = $self->make_seqfetcher("protein_index", "seqfetch_conf");
  
  my $sim = new Bio::EnsEMBL::Pipeline::RunnableDB::Contig_BlastMiniGenewise(
									  -dbobj      => $self->dbobj,
									  -input_id   => $contig->id,
									  -seqfetcher => $seqfetcher,
#									  -analysis   => $analysis,
									 );
  $self->similarity_runnable($sim);
}

=head2 similarity_runnable

 Title   : similarity_runnable
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub similarity_runnable{
  my ($self,$arg) = @_;

  if (!defined($self->{'_similarity_runnables'})) {
      $self->{'_similarity_runnables'} = [];
  }
  
  if (defined($arg)) {
      if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	  push(@{$self->{'_similarity_runnables'}},$arg);
      } else {
	  $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
      }
  }
  
  return @{$self->{'_similarity_runnables'}};  

}

=head2 make_riken_runnable

 Title   : make_riken_runnable
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_riken_runnable {
  my ($self) = @_;
 
#  analysis setting up? 
#  my $analysis = $self->dbobj->get_Analysis_Adaptor->??;
  my $seqfetcher = $self->make_seqfetcher("riken_index", "riken_conf");

  if($seqfetcher->isa("Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch")){
    $self->warn("no riken_index defined in GB_conf::riken_conf - cannot run Riken_BlastMiniGenewise\n");
    return 0;
  }



  my $rik = new Bio::EnsEMBL::Pipeline::RunnableDB::Riken_BlastMiniGenewise(
									    -dbobj      => $self->dbobj,
									    -input_id   => $self->input_id,
									    -seqfetcher => $seqfetcher,
#									    -analysis   => $analysis,
									   );
  $self->riken_runnable($rik);
}

=head2 riken_runnable

 Title   : riken_runnable
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub riken_runnable {
  my ($self, $arg) = @_;

  if (!defined($self->{'_riken_runnables'})) {
      $self->{'_riken_runnables'} = [];
  }
  
  if (defined($arg)) {
      if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	  push(@{$self->{'_riken_runnables'}},$arg);
      } else {
	  $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
      }
  }
  
  return @{$self->{'_riken_runnables'}};  

}

=head2 make_genebuild_runnable

 Title   : make_genebuild_runnable
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_genebuild_runnable{
  my ($self) = @_;
 
#  analysis setting up? 
#  my $analysis = $self->dbobj->get_Analysis_Adaptor->??;
  # assumes that fpc contigs have been generated using pipeline/scripts/make_clone_gp
  my $fpcctg = "ctg." . $self->input_id;

  my $gbr = new Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Gene_Builder(
								  -dbobj    => $self->dbobj,
								  -input_id => $fpcctg,
#								  -analysis => $analysis,
								 );
  $self->genebuild_runnable($gbr);
}


=head2 genebuild_runnable

 Title   : genebuild_runnable
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub genebuild_runnable {
  my ($self,$arg) = @_;

  if (!defined($self->{'_genebuild_runnables'})) {
      $self->{'_genebuild_runnables'} = [];
  }
  
  if (defined($arg)) {
      if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	  push(@{$self->{'_genebuild_runnables'}},$arg);
      } else {
	  $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
      }
  }
  
  return @{$self->{'_genebuild_runnables'}};  
}


=head2 run

 Title   : run
 Usage   :
 Function: runs the various runnables controlling the whole gene build in the right order
 Example :
 Returns : 
 Args    :


=cut

sub run {
  my ($self) = @_;
#$self->throw("bailing before run\n");

print STDERR "***Running targetted build***\n";
  foreach my $tge($self->targetted_runnable){
    $tge->fetch_input;
    $tge->run;
   $tge->write_output;
  }

print STDERR "***Running similarity build***\n";
  foreach my $sgw($self->similarity_runnable){
    $sgw->fetch_input;
    $sgw->run;
    $sgw->write_output;
  }

#print STDERR "***Running riken build***\n";
#  foreach my $rgw($self->riken_runnable){
#    $rgw->fetch_input;
#    $rgw->run;
#    $rgw->write_output;
#  }

print STDERR "***Running final build***\n";
  foreach my $gb($self->genebuild_runnable){
    $gb->fetch_input;
    $gb->run;
    $gb->write_output;
  }

}

=head2 output

 Title   : output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub output{
  my ($self) = @_;
  # nothing to output ...
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   
    Args    :   

=cut


sub write_output {
  my ($self) = @_;
  # data has already been written out, so no need to do anything here.
}

sub convert_output {
  my ($self) = @_;
  # nothing to do
}

1;
