#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneWise.pm
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneWise

=head1 SYNOPSIS

=head1 DESCRIPTION

Runs all the TargettedGeneWise jobs needed for a chunk of genomic sequence

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_TargettedGeneWise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::PmatchFeature;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBNAME
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Sequences qw (
							     GB_PROTEIN_INDEX
							     GB_PROTEIN_SEQFETCHER
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General   qw (
							     GB_INPUTID_REGEX
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts    qw (
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
  my $protein_fetcher = $self->make_seqfetcher($GB_PROTEIN_INDEX, $GB_PROTEIN_SEQFETCHER);

  # we need to find all the proteins that pmatch into this region
  # take a note of those that fall across the ends of the vc? and do what, precisely?
  # extend the VC? that will completely screw up the final genebuild. Hmmm.
  # do it, track it & see how many are affected.
  #input_id cb25.fpc4118.1-298757 has invalid format - expecting chr_name.start-end
  my $pipeline_db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(-dbname => $self->db->dbname,
								-user => $self->db->username,
								-pass => $self->db->password,
								-host => $self->db->host);

  my $pmfa = $pipeline_db->get_PmatchFeatureAdaptor();
  #print STDERR "have ".$pmfa." adaptor\n";
  my $input_id = $self->input_id;
  my $msg = "input_id $input_id has invalid format - expecting chr_name.start-end";
  $self->throw($msg) unless $input_id =~ /$GB_INPUTID_REGEX/;
  
  my $chr_name = $1;
  my $start   = $2;
  my $end     = $3;

  #print STDERR "fetching features for $chr_name $start $end\n";
  my @features = $pmfa->get_PmatchFeatures_by_chr_start_end($chr_name, $start, $end);
  #print STDERR "have ".@features." to run with\n";
  my $genewise_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						       '-host'   => $GB_GW_DBHOST,
						       '-user'   => $GB_GW_DBUSER,
						       '-pass'   => $GB_GW_DBPASS,
						       '-dbname' => $GB_GW_DBNAME,
						      );
  
  
  $genewise_db->dnadb($self->db);
  $self->output_db($genewise_db);
  my %kill_list = %{$self->fill_kill_list};

  foreach my $feat($pmfa->get_PmatchFeatures_by_chr_start_end($chr_name, $start, $end)){
    
    #reject any proteins that are in the kill list
    if(defined $kill_list{$feat->protein_id}){
      print STDERR "skipping " . $feat->protein_id . "\n";
      next;
    }


    my $input = $feat->chr_name   . ":" . 
                $feat->start      . "," . 
                $feat->end        . ":" . 
                $feat->protein_id . ":" ;

    #print STDERR "TGW input: $input\n";

    my $tgr = new Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise(
									-db => $self->db,
									-input_id => $input,
									-seqfetcher => $protein_fetcher,
									-analysis => $self->analysis,
									-output_db => $self->output_db,
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
  my $count = 0;
  #print STDERR "***Running targetted build***\n";
 TGE:   
  foreach my $tge($self->targetted_runnable){
    
    eval{
      $tge->fetch_input;
    };
    
    if($@){
      $self->warn("problems fetching input for TargettedGeneWise: [$@]\n");
      next TGE;
    }

    eval{
      $tge->run;
    };

    if($@){
      $self->warn("problems running TargettedGeneWise: [$@]\n");
      next TGE;
    }
    
    eval{
     my @genes = $tge->write_output;
     $count += @genes;
    };

    if($@){
      $self->warn("problems writing output for TargettedGeneWise: [$@]\n");
      next TGE;
    }
  }
  
  print "there were ".$count." genes produce\n";
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
