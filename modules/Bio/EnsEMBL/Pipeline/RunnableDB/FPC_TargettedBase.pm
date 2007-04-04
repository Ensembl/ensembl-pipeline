
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedBase.pm
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedBase

=head1 SYNOPSIS

=head1 DESCRIPTION

Runs all the TargettedGenewise jobs needed for a chunk of genomic sequence

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedBase;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw( throw ) ; 
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Sequences
  qw (
      GB_PROTEIN_SEQFETCHER
     );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts
  qw (
      GB_KILL_LIST
     );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
			   -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis);

    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneWise object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Args    :   -db:      A Bio::EnsEMBL::DBSQL::DBAdaptor,
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 

=cut

sub new {
  my ($class, @args) = @_;
  #print STDERR "args @args\n";
  my $self = $class->SUPER::new(@args);

  # db, input_id, seqfetcher, and analysis objects are all set in
  # in superclass constructor (RunnableDB.pm)

  return $self;
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: if $index exists, 
           returns a Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs, otherwise throws
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    : $indexname - string


=cut

sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;

  my $seqfetcher;
  
  (my $class = $seqfetcher_class) =~ s/::/\//g;

  throw ("Configuration-error !! There's no entry for GB_PROTEIN_SEQFETCHER in Sequences.pm\n") if (length($class)==0)  ;

  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );
    
    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );
  }
  else{
    $self->throw("can't make seqfetcher\n");
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
  my $count = 0;
  #print STDERR "***Running targetted build***\n\n";
 TARGETTED:
  foreach my $tge($self->targetted_runnable){

    eval{
      $tge->fetch_input;
    };

    if($@){
      $self->warn("problems fetching input for Targetted run (TargettedGenewise or TargettedExonerate): [$@]\n");
      next TARGETTED;
    }

    eval{
      $tge->run;
    };

    if($@){
      $self->warn("problems running Targetted run (TargettedGenewise or TargettedExonerate): [$@]\n");
      next TARGETTED;
    }

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
  my @output;
  foreach my $tge($self->targetted_runnable){
    push(@output, $tge->output);
  }
  return @output;
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
 TGE:
  foreach my $tge($self->targetted_runnable){
    eval{
      my @genes = $tge->write_output;
    };
    
    if($@){
      $self->warn("problems writing output for Targetted: [$@]\n");
      next TGE;
    }
  }
  # data has already been written out, so no need to do anything here.
}

sub convert_output {
  my ($self) = @_;
  # nothing to do
}


sub output_db {
    my( $self, $output_db ) = @_;

    if ($output_db) 
    {
	$output_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
	    || $self->throw("Input [$output_db] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
	$self->{_output_db} = $output_db;
    }
    return $self->{_output_db};
}


=head2 fill_kill_list

 Title   : fill_kill_list
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub fill_kill_list {
  my ($self) = @_;
  my %kill_list;
  open (KILL_LIST, "< $GB_KILL_LIST") or die "can't open $GB_KILL_LIST";
  while (<KILL_LIST>) {

    chomp;
    my @list = split;
    next unless scalar(@list); 	# blank or empty line
    $kill_list{$list[0]} = 1;
  }

  close KILL_LIST or die "error closing $GB_KILL_LIST\n";

  return \%kill_list;
}


1;
