
=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::SlamDB

=head1 SYNOPSIS

  get a Bio::EnsEMBL::Pipeline::RunnableDB::SlamDB object:

  $obj = new Bio::EnsEMBL::Pipeline::RunnableDB::SlamDB (
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis 
                                                       );

  $slamdb->fetch_input();
  $slamdb->run();
  $slamdb->output();
  $slamdb->write_output();

=head1 DESCRIPTION

 This object wraps Bio::EnsEMBL::Pipeline::Runnable::Slam (and uses also Avid and ApproxAlign)
 to add functionality to read and write to databases.
 A Bio::EnsEMBL::Pipeline::DBSQL::Obj is required for databse access.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut






=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

  package Bio::EnsEMBL::Pipeline::RunnableDB::SlamDB;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::DB::RandomAccessI;

use Bio::EnsEMBL::Pipeline::Runnable::Avid;
use Bio::EnsEMBL::Pipeline::Runnable::Slam;
use Bio::EnsEMBL::Pipeline::Tools::ApproxAlign;


# vars from Slamconf-File
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Slamconf qw (
                                                            SLAM_COMP_DB_USER
                                                            SLAM_COMP_DB_PASS
                                                            SLAM_COMP_DB_NAME
                                                            SLAM_COMP_DB_HOST
                                                            SLAM_COMP_DB_PORT
                                                           );



use vars qw(@ISA);


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;

  my $self = {};
  bless $self, $class;

  my ($db, $input_id, $seqfetcher, $analysis) = $self->_rearrange([qw(
                                                                      DB
                                                                      INPUT_ID
                                                                      SEQFETCHER
                                                                      ANALYSIS )],
                                                                  @args);

  &throw("No database handle input for first organsim\n") unless defined($db);
  &throw("No analysis object input") unless defined($analysis);

  $self->analysis($analysis);
  $self->store_regions($input_id);

  $self->db($db);          #super db()
  $self->db_comp;

  return $self;
}
##   one db can be used as "base-database" like homo sapiens,
##   the other database can read out of the config-file.

##    But how to passs the diffrent slice-ids to the runnable ?
##    ==> Here we can code the two input-id's in one string and split the string.

##   nice example GeneBuilder !!!
##
##  ./modules/Bio/EnsEMBL/Pipeline/Config/GeneBuild/GeneBuilder.pm.example
##   use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder q
##   building an id like CHRNAME_start-end




sub fetch_input {
  my $self = shift;

  $self->fetch_slice1( $self->db, $self->base_region);
  $self->fetch_slice2( $self->db_comp, $self->comp_region);
}


sub fetch_slice1{
  my ($self,$db,$region) = @_;

  my $sa = $db->get_SliceAdaptor();
  my @items = split/_/, $region; # region-format CHR_START-END
  my $chr = $items[0];
  my ($start, $end)  = split/-/, $items[1];
  my $slice = $sa->fetch_by_region('chromosome' , $chr, $start, $end) ;
  $self->{_slice1} = $slice;

  return $self->{_slice1};
}

sub fetch_slice2{
  my ($self,$db,$region) = @_;

  my $sa = $db->get_SliceAdaptor();
  my @items = split/_/, $region; # region-format CHR_START-END
  my $chr = $items[0];
  my ($start, $end)  = split/-/, $items[1];
  my $slice = $sa->fetch_by_region('chromosome' , $chr, $start, $end) ;
  $self->{_slice2} = $slice;

  return $self->{_slice2};
}


=head2 db_comp

    Title   :   db_comp
    Usage   :   $self->db_comp($obj);
    Function:   Gets or sets the value of db_comp
    Returns :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor compliant object
                (which extends Bio::EnsEMBL::DBSQL::DBAdaptor)
    Args    :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor compliant object

=cut

sub db_comp {
  my( $self) = shift;

  my  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                                -user   => $SLAM_COMP_DB_USER,
                                                -dbname => $SLAM_COMP_DB_NAME,
                                                -host   => $SLAM_COMP_DB_HOST,
                                                -pass   => $SLAM_COMP_DB_PASS,
                                                -port   => $SLAM_COMP_DB_PORT,
                                                -driver => 'mysql'
                                               );
  $self->{'_db_comp'} = $db;
  return $self->{'_db_comp'};
}





sub store_regions {
  my ($self,$input_id) = @_;

  if (defined $input_id) {
    my @input = split /---/,$input_id;
    $self->base_region($input[0]);
    $self->comp_region($input[1]);
  }
}


sub base_region {
  my ($self,$input) = @_;

  if (defined $input) {
    $self->{_base_region} = $input;
  }
  return $self->{_base_region};
}

sub comp_region {
  my ($self,$input) = @_;

  if (defined $input) {
    $self->{_comp_region} = $input;
  }
  return $self->{_comp_region};
}

1;


