
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
                                                            SLAM_MAXLENGTH
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
  $self->regions($input_id);
  $self->db($db);               #super db()
  $self->db_org2;
  return $self;
}


# gets a reference to an array of slices to run avid on
# returns reference to an avid-object
sub avid{
  my ($self,$slices) = @_;

  my $avid =  new Bio::EnsEMBL::Pipeline::Runnable::Avid (
                                                          -slice1      => ${$slices}[0],
                                                          -slice2      => ${$slices}[1],
                                                         );
  $avid->run;
  return $avid;
}


# gets name of parsed binary and fastanames
# returns name of written modified approximate alignment
sub approx_align{
  my ($parsed_bin,$fasta1,$fasta2);

  my $approx_obj = new Bio::EnsEMBL::Pipeline::Tools::ApproxAlign(
                                                                   -aat =>  $parsed_bin,   # /path/to/parsedbinaryfile
                                                                   -seqY => $fasta1,       # /path/to/firstfasta.fasta
                                                                   -seqZ => $fasta2        # /path/to/secondfasta.fasta
                                                                  );
  $approx_obj->expand($approx_obj->exonbounds);
  $approx_obj->makeConsistent();

  if ($approx_obj->isConsistent) {
    my $aatfile = $approx_obj->write(); #can perhaps return the filename instead of getting one
    print "data written to $aatfile\n"  ;
  } else {
    die "Error: final aat is not consistent (shouldn't have happened).\n"
  }
  return $approx_obj;
}


sub run{
  my ($self) = shift;
  my @subslices;

  ####  if seqlength > maxlength we have to split seq ###

  if ( (${$self->slices}[0]->length  || ${$self->slices}[0]->length ) > $SLAM_MAXLENGTH) {

    # run avid on original slices
    my $avid = &avid(${$self->slices}[0],${$self->slices}[1] );

    # run ApproxAlign on org slices
    my $approx = &approx_align($avid->parsed_binary_filename,$avid->fasta_filename1,$avid->fasta_filename2);

    # cut the first seq according to the positions of the repeats
    # and than try to find equal positions in the second seq by
    # using the approximate Alignment-aatfile

    my @cuts = @{ $self->calculate_cutting ($approx->aatfile) } ;

    # we got the cutting positions in @cuts, we build & store the subslices in @subslices
    # @cuts is an array of arrays [ [start1,end1,start2,end2],[start1,end1,start2,end2] ]
    for my $subseqs (@cuts) {
      # store the subslices in 2nd array of arrays
      my ($start1, $end1) = @{$subseqs}[0,1];
      my ($start2, $end2) = @{$subseqs}[2,3];

      my $subslice1 = ${$self->slices}[0] -> sub_Slice( $start1, $end1 );
      my $subslice2 = ${$self->slices}[1] -> sub_Slice( $start2, $end2 );
      push @subslices, [$subslice1,$subslice2];
    }
  }else{
    # we don't need any cutting
    push @subslices, $self->slices;
  }

  for my $slices (@subslices) {

    my $avid = &avid(${$slices}[0],${$slices}[1] );

    # run ApproxAlign on org slices
    my $approx = &approx_align($avid->parsed_binary_filename,$avid->fasta_filename1,$avid->fasta_filename2);


    # make new slam-run with subslice
    my $slamobj = new Bio::EnsEMBL::Pipeline::Runnable::Slam (
                                                            -slice1      => ${$slices}[0],
                                                            -slice2      => ${$slices}[1],
                                                            -fasta1        => $avid->fasta_filename1,
                                                            -fasta2        => $avid->fasta_filename2,
                                                            -approx_align  => $approx->aafile,
                                                            -org1          => $SLAM_ORG1_NAME,
                                                            -org2          => $SLAM_ORG2_NAME,
                                                            -slam_bin      => $SLAM_BIN,
                                                            -slam_pars_dir => $SLAM_PARS_DIR,
                                                            -minlength     => $SLAM_MINLENGTH,
                                                            -debug         => 0,
                                                            -verbose       => 0
                                                           );
  # run slam, parse results
    $slamobj->run;

  # set ref to arrays with predicted transcripts for both organisms
    $self->predtrans_both_org ( $slamobj ->predtrans );

  # get back the parsed gff->hang the gff's together
  # write the results to database
  # AND WHAT ABOUT THE COORDINATES ???
  }


  # POSSIBILITES :
  # compare repeats of first seq in db with repeats after RM-run
  # transfer db-repeats in RM-outfile for first seq 
  # or 
  # transfer RM-outfile-repeats in array which is used by this script (format START - END)
  #
  # cut the sequence in diffrent parts --- but do we have to to it ? what are the needs for it ?
  # Which programs are working with it ?
  # Which programs need to read the files, which get slices ?
  # Why is the RM done ? (logic in slam.pl)
  # AND WHAT ABOUT THE COORDINATES ???



}  # end run




  ### make cuts according to position of repeats in first sequence

    # output-format of aat-file (one row for each base, slam needs the same input)
    #   base lowerBound upperBound
    #   0        0          881
    #   1        2          882
    #   2        843        883
    #   3        843        884
    #   4        849        885
    #   5        850        886
    #   6        851        887
    #   7        852        888
    #   8        853        889
    #   9        856        898


sub calculate_cutting{
  my ($self,$ApproxAlign) = @_;

  # getting attributes of object
  my @slices  = @{$self->slices};
  my $targetcut = $SLAM_MAXLENGTH;
  my @all_repeats = @{$self->get_repeat_features};

  my @cuts1 = (1);
  my $len1 = $slices[0]->length;

  while ($targetcut < $len1) {

    my $cut = undef;
    while (1) {
      # $all_repeats[0]->[0] = start of repeat
      # $all_repeats[0]->[1] = end of repeat

      if ((@all_repeats==0) || ($all_repeats[0]->[0] > $targetcut)) {

        # No repeats or startpos of first repeat is bigger than targetcut
        # so there are no repeats before target-cuttingposition, so cut
        last;

      }elsif ($all_repeats[0]->[1] >= $targetcut) {
        # end of repeat is "bigger" than targetcut
        # repeat spans target (cool), see example1 above
        $cut = $targetcut;
        last;

      }else {
        # Store end of repeat as the best-yet value, then move on.
        $cut = $all_repeats[0]->[1];
        shift(@all_repeats);
      }
    }

    if ( (!defined($cut)) || (($targetcut-$cut+1) > ($SLAM_MAXLENGTH/2)) ) {
      # If no repeats before targetcut or cut too far away, then cut at target anyway
      $cut = $targetcut;
    }
    push(@cuts1,$cut);
    $targetcut = $cut + $SLAM_MAXLENGTH;
  } # while(1)

  # last cut is length of seq
  # now we got cutting-positions for the first sequence

  push(@cuts1,$len1) if($cuts1[$#cuts1] < $len1);

  ################################################################################
  #
  # Make array of matching cuts in other sequence using the approximate alignement
  #
  # what is the lower bound for the first cut ? What generally is a lowerBound ?

  # get first cut for second sequence
  my @cuts2 = (1+$ApproxAlign->lowerBound($cuts1[0]-1));


  for (my $i=1; $i < (scalar(@cuts1)-1); $i++) {
    push(@cuts2,sprintf("%d",($ApproxAlign->lowerBound($cuts1[$i]-1) + $ApproxAlign->upperBound($cuts1[$i]-1))/2.0));
  }
  push(@cuts2,(1+$ApproxAlign->upperBound($cuts1[$#cuts1]-1)));



  $cuts1[0] = $cuts1[0]-1;
  $cuts2[0] = $cuts2[0]-1;

  # now the splits

  my @splits = ();                   # nr of cuts
  for(my $i=0, my $cutCount=0; $i < (scalar(@cuts1)-1); $i++) {
    if($cuts2[$i]+1 > $cuts2[$i+1]) {
      # skip if we have an insertion in the base seqeunce.
      next;
    } else {
      $cutCount++;                                                      # cuts in the first seq   # cuts in the second seq
##      push(@splits,[sprintf("%s.%d","mycutfile_dir",$cutCount),$cuts1[$i]+1,$cuts1[$i+1],$cuts2[$i]+1,$cuts2[$i+1]]);
      push(@splits,[ $cuts1[$i]+1, $cuts1[$i+1], $cuts2[$i]+1, $cuts2[$i+1] ] );
      # Format of cutfile:
      # human_contig.fasta_mice_contig.fasta.cut.1      1       100500  1       93943
      # human_contig.fasta_mice_contig.fasta.cut.2      100501  119071  93944   100090
    }
  }
  return(\@splits);
}






############################################################

sub fetch_input {
  my $self = shift;
  $self->slices([$self->db,$self->db_org2]);
}


# fetching slices for each org out of specified db

sub slices{
  my ($self ,$db) = @_;

  if ($db) {
    my @slices;
    my @coords = @{$self->regions};
    for (my $i=0;$i<=1;$i++) {
      my $sa = ${$db}[$i]->get_SliceAdaptor();
      my ($chr,$start,$end) = splice(@coords, 0,3);
      my $slice = $sa->fetch_by_region('chromosome' , $chr, $start, $end) ;
      push @slices, $slice;
    }
    $self->{_slices}=\@slices;
  }
  return $self->{_slices};
}

# splits input-id and sets the diffrent regions
sub regions {
  my ($self,$input_id) = @_;

  if (defined $input_id) {
    my @input = split /---/,$input_id; # format chr2-start1-end1---chr2-start2-end2
    my ($chr1, $start1, $end1) = split/-/, $input[0];
    my ($chr2, $start2, $end2) = split/-/, $input[1];
    $self->{_regions}=[$chr1, $start1, $end1, $chr2, $start2, $end2];
  }
  return $self->{_regions};
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


=head2 predtrans_both_org

  Title    : predtrans
  Usage    : $obj->predtrans
  Function : Sets/gets the Predicted transcripts for the first organism
  Returns  : Ref. to an Array of Arrayrefs. to Bio::EnsEMBL::PredictionTranscript
  Args     : References to two arrays

=cut

sub predtrans_both_org{
  my ($self,$ref_predtrans) = @_;

  if ($ref_predtrans) {
    $self->{_ref_predtrans} = $ref_predtrans;
  }
  return $self->{_ref_predtrans};
}


sub write_output {
  my ($self) = @_;

  #writing output for both organisms to diffrent databases
  $self->write_dbresults ( $self->db,${$self->slices}[0], ${$self->predtrans_both_org}[0] );
  $self->write_dbresults ( $self->db_org2, ${$self->slices}[1], ${$self->predtrans_both_org}[1] );
}



# writing the results to the given database (mouse/human/rat)
# gets a database to write to, a slice and a reference to an array of predicted transcripts
# looks up for the analysis

sub write_dbresults {
  my ($self,$db,$slice,$pt) = @_;

  my $pred_adp = $db->get_PredictionTranscriptAdaptor;
  my @pred_trans = @{$pt};
  my $analysis = $self->analysis;
  foreach my $tr (@pred_trans) {
    $tr->analysis($analysis);
    foreach my $exon (@{$tr->get_all_Exons}) {
      $exon->slice($slice);
    }
  }
  $pred_adp->store(@pred_trans);
}





#                           REPEATS AND SPLITTING
# get the start/end-positions of repeats in the first sequence (LTRs, LINEs and SINEs)
# out of the db and look for a good place to cut. Good positions are if the expexted
# cuttingposition (which is the max. length of seqslice to compare) lies in a region 
# with repeats (ex1)
# or
#
# ex1:  ..._repeat_repeat_repeat_CUTTINGPOSITION_repeat_repeat_repeat_...
#
# ex2:  ..._

sub get_repeat_features {
  my $self = shift;

  my @slices = @{$self->slices};
  my (@all_rpt);

  # repeat-types to look for
  my @repeats =('LTRs','Type I Transposons/LINE','Type I Transposons/SINE');

  my @all;
  # get all repeats of the given types (above)
  for my $arpt(@repeats){
    for my $rpt ( @{$slices[0]->get_all_RepeatFeatures(undef,"$arpt")}) {
      my $rpt_start = $rpt->start;
      my $rpt_end = $rpt->end;
      push (@all_rpt , [ $rpt_start , $rpt_end ] );
      push (@all , [ $rpt_start , $rpt_end , $rpt->display_id] );
    }
  }
# routine of comparing repeats of db with repeats from repeatmasker-run
#  for (@all){
#    my @r = @{$_};
#    for(@r){
#      print "\t $_";
#    }
#    print "\n";
#  }
  # sort startpos of rpt ascending order
  @all_rpt = sort { $a->[0] <=> $b->[0] } @all_rpt;
  # check the sorting
#  for my $r (@all_rpt){
#    my @r=@{$r};
#    for(@r){
#      print "\t $_";
#    }
#    print "\n";
#  }
  return \@all_rpt;
}
