#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneE2G.pm
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneE2G

=head1 SYNOPSIS

=head1 DESCRIPTION

Runs all the TargettedGeneE2G jobs needed for a chunk of genomic sequence

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneE2G;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneE2G;
use Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::PmatchFeature;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_GOLDEN_PATH
					 GB_TARGETTED_PROTEIN_INDEX
					 GB_TARGETTED_CDNA_INDEX
					);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
			   -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneE2G object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  # dbobj, input_id, seqfetcher, and analysis objects are all set in
  # in superclass constructor (RunnableDB.pm)

  # golden path
  my $path = $GB_GOLDEN_PATH;
  $path = 'UCSC' unless (defined $path && $path ne '');
  $self->dbobj->assembly_type($path);

  return $self;
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: if $index exists, 
           returns a Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs, otherwise 
           returns a Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    : $indexname - string


=cut

sub make_seqfetcher {
  my ( $self, $index ) = @_;

  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs(
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

  $self->make_targetted_runnables;
}

=head2 make_targetted_runnables

 Title   : make_targetted_runnables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub make_targetted_runnables {
  my ($self) = @_;

  # set up seqfetchers
  my $protein_fetcher = $self->make_seqfetcher($GB_TARGETTED_PROTEIN_INDEX);
  my $cdna_fetcher    = $self->make_seqfetcher($GB_TARGETTED_CDNA_INDEX);

  # we need to find all the proteins that pmatch into this region
  # take a note of those that fall across the ends of the vc? and do what, precisely?
  # extend the VC? that will completely screw up the final genebuild. Hmmm.
  # do it, track it & see how many are affected.

  my $pmfa = new Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor( $self->dbobj );
  my $input_id = $self->input_id;
  my $msg = "input_id $input_id has invalid format - expecting chr_name.start-end";
  $self->throw($msg) unless $input_id =~ /(\w+)\.(\d+)-(\d+)/;
  
  my $chr_name = $1;
  my $start   = $2;
  my $end     = $3;

print STDERR "fetching features for $chr_name $start $end\n";

  foreach my $feat($pmfa->get_PmatchFeatures_by_chr_start_end($chr_name, $start, $end)){
    my $input = $feat->chr_name   . ":" . 
                $feat->start      . "," . 
                $feat->end        . ":" . 
                $feat->protein_id . ":" . 
                $feat->cdna_id;

    print STDERR "TGE input: $input\n";

    my $tgr = new Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneE2G(
								       -dbobj         => $self->dbobj,
								       -input_id      => $input,
								       -seqfetcher    => $protein_fetcher,
								       -cdna_seqfetcher => $cdna_fetcher,
#								       -analysis => $analysis,
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

print STDERR "***Running targetted build***\n";
  foreach my $tge($self->targetted_runnable){
    $tge->fetch_input;
    $tge->run;
    $tge->write_output;
  }

$self->{'_targetted_runnables'} = [];

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
