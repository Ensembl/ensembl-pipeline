
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


# vars read from Slamconf-File
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Slamconf qw (
                                                            SLAM_ORG1_NAME
                                                            SLAM_ORG2_NAME
                                                            SLAM_BIN
                                                            SLAM_PARS_DIR
                                                            SLAM_MINLENGTH
                                                            SLAM_COMP_DB_USER
                                                            SLAM_COMP_DB_PASS
                                                            SLAM_COMP_DB_NAME
                                                            SLAM_COMP_DB_HOST
                                                            SLAM_COMP_DB_PORT
                                                            SLAM_ORG2_RESULT_DB_USER
                                                            SLAM_ORG2_RESULT_DB_PASS
                                                            SLAM_ORG2_RESULT_DB_NAME
                                                            SLAM_ORG2_RESULT_DB_HOST
                                                            SLAM_ORG2_RESULT_DB_PORT
                                                           );



use vars qw(@ISA);


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

############################################################

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


  $self->db($db);               #super db()
  $self->db_org2;

  return $self;
}

###============================== RUN ==============================

sub run{
  my ($self) = shift;

  my $avid =  new Bio::EnsEMBL::Pipeline::Runnable::Avid (
                                                          -slice1      => $self->fetch_slice1,
                                                          -slice2      => $self->fetch_slice2,
                                                          -slam_output => 'true',
                                                         );

  $avid->run;

  # modify approximate alignment
  my $ApproxAlign = new Bio::EnsEMBL::Pipeline::Tools::ApproxAlign(
                                                                   -aat =>  $avid -> parsed_binary_filename, # name of parsed binary
                                                                   -seqY => $avid -> fasta_filename1, # first filename
                                                                   -seqZ => $avid -> fasta_filename2 # second filename
                                                                  );

  ####=================================== START ApproxAlign ===================================#
  #### OPTIONS for approximate Aligement - perhaps we can put this in a configuration-file of SlamDB or ApproxAlign ?

  my $acc = 20;                 # acc-slam-default: 20
  my $don = 10;                 # don-slam-default: 10
  my $stp = 5;                  # stp-slam-default: 5
  my $fatten = undef;           # fatten-default  : undef
  my $fwdOnly = 0;              # fwdOnly         : 0
  my $weight1 = $ApproxAlign->weight();

  ### EXON BOUNDARIES
  my $fwdStp1 = ['TAA', 0, 'zOnly', $stp];
  my $fwdStp2 = ['TAG', 0, 'zOnly', $stp];
  my $fwdStp3 = ['TGA', 0, 'zOnly', $stp];
  my $revStp1 = ['TTA', 2, 'zOnly', $stp];
  my $revStp2 = ['CTA', 2, 'zOnly', $stp];
  my $revStp3 = ['TCA', 2, 'zOnly', $stp];
  #my $fwdDon = ['GT',  2, 'points',  $don];
  #my $revDon = ['AC', -1, 'points', -$don, -($don+1) ],
  #my $fwdAcc = ['AG', -1, 'points', -$acc, -($acc+1) ];
  #my $revAcc = ['CT',  2, 'points',  $acc];
  my $fwdDon = ['GT',  2, 'line',   $don];
  my $revDon = ['AC', -1, 'line', -($don + 1)]; #10-->-11 WATCH
  my $fwdAcc = ['AG', -1, 'line', -($acc + 1)];
  my $revAcc = ['CT',  2, 'line',   $acc];
  my $exonBounds;

  if ($fwdOnly) {               # fwdOnly normally not used
    $exonBounds = [
                   $fwdStp1, $fwdStp2, $fwdStp3,
                   $fwdDon, $fwdAcc
                  ];
  } else {
    $exonBounds = [
                   $fwdStp1, $fwdStp2, $fwdStp3,
                   $revStp1, $revStp2, $revStp3,
                   $fwdDon, $revDon,
                   $fwdAcc, $revAcc
                  ];
  }
  # not used because slam-default doesn't fatten
  push(@{$exonBounds}, ['fatten', undef, undef, $fatten]) if(defined($fatten));

  $ApproxAlign->expand($exonBounds);
  my $weight2 = $ApproxAlign->weight();
  $ApproxAlign->makeConsistent();
  my $weight3 = $ApproxAlign->weight();

  printf STDERR ("Raw aat weight:          %6d\n",$weight1) ;
  printf STDERR ("Intermediate aat weight: %6d\n",$weight2) ;
  printf STDERR ("Final aat weight:        %6d\n",$weight3) ;

  my $aat_filename = Bio::EnsEMBL::Pipeline::RunnableI->get_tmp_file("/tmp","approxAlignOutput","aat"); #CHANGE THIS!!!!!!!!!!!!1

  if ($ApproxAlign->isConsistent) {
    $ApproxAlign->write("$aat_filename"); #can perhaps return the filename instead of getting one
    print "data written to $aat_filename\n"  ;
  } else {
    die "Error: final aat is not consistent (shouldn't have happened).\n"
  }
  ##=================================== END ===================================#


  print "aproxalignfile -$aat_filename-\n"  ;

  my $slamobj = new Bio::EnsEMBL::Pipeline::Runnable::Slam (
                                                            -slice1        => $self->fetch_slice1,
                                                            -slice2        => $self->fetch_slice2,
                                                            -fasta1        => $avid->fasta_filename1,
                                                            -fasta2        => $avid->fasta_filename2,
                                                            -approx_align  => $aat_filename,
                                                            -org1          => $SLAM_ORG1_NAME,
                                                            -org2          => $SLAM_ORG2_NAME,
                                                            -slam_bin      => $SLAM_BIN,
                                                            -slam_pars_dir => $SLAM_PARS_DIR,
                                                            -minlength     => $SLAM_MINLENGTH,
                                                            -debug         => 0,
                                                            -verbose       => 0

                                                           );

  # running slam and parsing results

  $slamobj->run;

  $self->predicted_transcripts_org1 ( $slamobj ->predtrans1);
  $self->predicted_transcripts_org2 ( $slamobj ->predtrans2);
}


# attach analysis and slice and write the output to the database



sub write_output {
  my ($self) = @_;

  #writing output for first organism
  print "in writeoutput\n";
  $self->write_output_new($self->db, $self->fetch_slice1, $self->predicted_transcripts_org1 );
  $self->write_output_new($self->db_org2, $self->fetch_slice2, $self->predicted_transcripts_org2 );
}



# caller write_output_new ( $db, $analysis,  )
sub write_output_new {
  my ($self,$db,$slice,$pt) = @_;
  print "in write_output_new\n";
  my $pred_adp = $db->get_PredictionTranscriptAdaptor;
  my @pred_trans = @{$pt};


#  if ($#pred_trans>0) {
    # we have some predictionTranscript, so let's put'em in the db
#    print STDERR "Write output have ".@pred_trans." features\n";
    my $analysis = $self->analysis;

#    my @manipulate;
#    foreach my $trans (@pred_trans) {
#      print "new feature $trans\n";
#      $trans->analysis($analysis);
#      push @manipulate,$trans;
#    }
#    foreach my $tr (@manipulate) {

    foreach my $tr (@pred_trans) {
      $tr->analysis($analysis);
      foreach my $exon (@{$tr->get_all_Exons}) {

        $exon->slice($slice);
#        my $analysis = $self->analysis;
#        $exon->analysis($analysis);
      }
    }

    $pred_adp->store(@pred_trans);
#  } else {
#    print "no output\n";
#  }
}





## attach analysis and slice and write the output to the database
#sub write_output_old {
#  my ($self) = @_;

#  my $db       = $self->db();
#  my $pred_adp = $self->db->get_PredictionTranscriptAdaptor;


#  my @features = @{$self->predicted_transcripts_org1};
#  # convert to diffrent slice !
#  my $slice = $self->fetch_slice1;



#  print STDERR "Write output have ".@features." features\n";

#  # foreach transcript
#  foreach my $f (@features) {
#    print "new feature $f\n";

#    # attaching analysis to Bio::EnsEMBL::Feature-object
#    $f->analysis($self->analysis);

#    foreach my $exon (@{$f->get_all_Exons}) {
#      print "new exon\n";
#      $exon->slice($slice);
#    }
#    $pred_adp->store(@features);
#  }
#}





sub predicted_transcripts_org1 {
  my ($self,$arrayref) = @_;

  if (defined $arrayref) {
    $self->{_store_pred_trans_org1} = $arrayref;
  }
  return $self->{_store_pred_trans_org1};
}


sub predicted_transcripts_org2 {
  my ($self,$arrayref) = @_;

  if (defined $arrayref) {
    $self->{_store_pred_trans_org2} = $arrayref;
  }
  return $self->{_store_pred_trans_org2};
}




















############################################################

sub fetch_input {
  my $self = shift;

  $self->fetch_slice1( $self->db, $self->base_region);
  $self->fetch_slice2( $self->db_org2, $self->comp_region);
}


############################################################

sub fetch_slice1{
  my ($self,$db,$region) = @_;

  if (defined $db && defined $region) {
    my $sa = $db->get_SliceAdaptor();
    my @items = split/_/, $region; # region-format CHR_START-END
    my $chr = $items[0];
    my ($start, $end)  = split/-/, $items[1];
    my $slice = $sa->fetch_by_region('chromosome' , $chr, $start, $end) ;
    $self->{_slice1} = $slice;
  }
  return $self->{_slice1};
}



############################################################

sub fetch_slice2{
  my ($self,$db,$region) = @_;
  if (defined $db && defined $region) {
    my $sa = $db->get_SliceAdaptor();
    my @items = split/_/, $region; # region-format CHR_START-END
    my $chr = $items[0];
    my ($start, $end)  = split/-/, $items[1];
    my $slice = $sa->fetch_by_region('chromosome' , $chr, $start, $end) ;
    $self->{_slice2} = $slice;
  }
  
  return $self->{_slice2};
}



############################################################

=head2 db_org2

    Title   :   db_org2
    Usage   :   $self->db_org2($obj);
    Function:   Gets or sets the value of db_org2
    Returns :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor org2liant object
                (which extends Bio::EnsEMBL::DBSQL::DBAdaptor)
    Args    :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor compliant object

=cut

sub db_org2 {
  my( $self) = shift;

  # data of db for writing results of second organism analysis (out of Conf/Genebuild/Slamconf.pm)

  my  $db_result_org2 = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                                            -user   => $SLAM_ORG2_RESULT_DB_USER,
                                                            -dbname => $SLAM_ORG2_RESULT_DB_NAME,
                                                            -host   => $SLAM_ORG2_RESULT_DB_HOST,
                                                            -pass   => $SLAM_ORG2_RESULT_DB_PASS,
                                                            -port   => $SLAM_ORG2_RESULT_DB_PORT,
                                                            -driver => 'mysql'
                                                           );

  # attaching dna-db for data retreival

  my  $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                                   -user   => $SLAM_COMP_DB_USER,
                                                   -dbname => $SLAM_COMP_DB_NAME,
                                                   -host   => $SLAM_COMP_DB_HOST,
                                                   -pass   => $SLAM_COMP_DB_PASS,
                                                   -port   => $SLAM_COMP_DB_PORT,
                                                   -driver => 'mysql'
                                                  );



  $db_result_org2 -> dnadb($dnadb);
  $self->{'_db_org2'} = $db_result_org2;

  return $self->{'_db_org2'};
}




############################################################

sub store_regions {
  my ($self,$input_id) = @_;

  if (defined $input_id) {
    my @input = split /---/,$input_id;
    $self->base_region($input[0]);
    $self->comp_region($input[1]);
  }
}

############################################################

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


