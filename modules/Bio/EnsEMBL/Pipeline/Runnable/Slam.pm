#
# Written by Jan-Hinnerk Vogel
#
# jhv [at] sanger.ac.uk
#
# Copyright GRL/EBI 2004
#
# You may distribute this module under the same terms as perl itself


=pod

=head1 NAME

 Bio::EnsEMBL::Pipeline::Runnable::Slam

=head1 SYNOPSIS

 To construct a new Bio::EnsEMBL::Pipeline::Runnable::Slam object use

 $obj = new Bio::EnsEMBL::Pipeline::Runnable::Slam (
                                                      -slice1        => $slice1,
                                                      -slice2        => $slice2,
                                                      -fasta1        => 'seq1.fasta',
                                                      -fasta2        => 'seq2.fasta',
                                                      -approx_align  => 'myaatfile.aat',
                                                      -org1          => 'H.sapiens',
                                                      -org2          => 'M.musculus',
                                                      -slam_bin      => '/nfs/acari/jhv/bin/slam'
                                                      -slam_pars_dir => '',
                                                      -minlength     => 250,
                                                      -debug         => 0,
                                                      -verbose       => 1
                                                      );

or

 $obj = Bio::EnsEMBL::Pipeline::Runnable::Slam->new(
                                                    -slice1        => $slice1,
                                                    -slice2        => $slice2,
                                                    -fasta1        => 'seq1.fasta',
                                                    -fasta2        => 'seq2.fasta',
                                                    -approx_align  => 'myaatfile.aat',
                                                    );



 $obj -> run;


=head1 DESCRIPTION

Slam takes two Fasta-Sequences and a modified approximate alignment-file
(see Bio::EnsEMBL::Pipeline::Tools::ApproxAlign). Presently, Slam-runs can
be performed on 3 diffrent organisms: H.sapiens, M.musculus and R.norvegicus.


=head1 OPTIONS

=over

=item  B<-slice1>        first slice (Bio::SeqIO-object)

=item  B<-slice2>        second slice (Bio::SeqIO-object)

=item  B<-fasta1>        /path/to/fastafile1

=item  B<-fasta2>        /path/to/fastafile2

=item  B<-approx_align>  path+filename of aat-file

=item  B<-org1>          species of first organism (def: H.sapiens)

=item  B<-org1>          species of second organism (def: M.musculus)

=item  B<-slam_bin>      path to slam-binary (def: /nfs/acari/jhv/bin/)

=item  B<-slam_pars_dir> path to slam-libary-directory 

=item  B<-minlength>     minium sequence length (def: 250 bp)

=item  B<-debug>         debug-option for slam-binary (def: 0)

=item  B<-verbose>       verbose-option for slam-run (def: 0)



=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 METHODS

The rest of the documentation details each of the object methods.
Internal (private) methods are usually preceded with a "_".

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Slam;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::PredictionExon;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;



# Inherits from Interface Bio::EnsEMBL::Pipeline::RunnableI
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
  my ($class,@args) = @_;
  my $self={};
  bless $self,$class;

  my (
      $slice1,                  # refrence to first slice
      $slice2,                  # refrence to second slice
      $fasta1,                  # filename 1st fasta-seq
      $fasta2,                  # filename 2nd fasta-seq
      $approx_align,            # filename modified approximate alignement (.aat-File)
      $org1,                    # name 1st organism (H.sapiens M.musculus R.norvegicus)
      $org2,                    # name 2nd organism (H.sapiens M.musculus R.norvegicus)
      $slam_bin,                # location slam_binary
      $slam_pars_dir,           # location of Pars-dir
      $minlength,               # min seq length
      $workdir,                 # opt. work-directory
      $debug,                   # opt. debug-option for slam
      $verbose                  # opt. verbose-option for runnable
     ) = $self->_rearrange([qw(
                               SLICE1
                               SLICE2
                               FASTA1
                               FASTA2
                               APPROX_ALIGN
                               ORG1
                               ORG2
                               SLAM_BIN
                               SLAM_PARS_DIR
                               MINLENGTH
                               WORKDIR
                               DEBUG
                               VERBOSE
                              )
                           ], @args
                          );

  $self->slice1($slice1);
  $self->slice2($slice2);


  $self->fasta1($fasta1);
  $self->fasta2($fasta2);

  $self->approx_align($approx_align);

  $self->org1($org1);
  $self->org2($org2);

  $self->slam_bin($slam_bin);
  $self->slam_pars_dir($slam_pars_dir);

  $self->minlength($minlength);

  $self->workdir($workdir);
  $self->debug($debug);
  $self->verbose($verbose);

  $self->printvars if ($self->verbose);
  return $self;
}



sub DESTROY {
  my $self = shift;
  #  $self->deletefiles;
}

=pod

=head2 _parse_results (String, String)

  Title    : run
  Usage    : $obj->_parse_ results
  Function : parses the gff-resultfile created in the slam-run
  Returns  : none
  Args     : 1 String which specify the file-basename of the fasta-file

=cut

# example of gff-file written by slam
# (see http://www.sanger.ac.uk/Software/formats/GFF/)

# <seqname> <source> <feature>   <start>  <end> <score> <strand> <frame> [attributes]   |    [comments]        |
#=========================part1========================================================   ---------part2------ | -----part3--------
#  $seq     $src     $feat        $start  $end  $score  $strand  $frame   $gid   $ginr  |   $tr_id      $tr_nr |   $ex     $extype

#  seq1      SLAM     CDS         18350   18443    .        +       1     gene_id "001" ; transcript_id "001.1"; exontype "internal"
#  seq2      SLAM     CDS         18525   18594    .        +       0     gene_id "001" ; transcript_id "001.1"; exontype "internal"
#  seq2      SLAM     CDS         18731   18747    .        +       2     gene_id "001" ; transcript_id "001.1"; exontype "terminal"
#  seq2      SLAM     stop_codon  18745   18747    .        +       .     gene_id "001" ; transcript_id "001.1";
#  seq2      SLAM     SLAM_CNS    18749   18816    .        .       .     cns_id "005"  ; identity "85.1"
#  seq2      SLAM     SLAM_CNS    20900   20986    .        .       .     cns_id "006"  ; identity "69.8"
#  seq2      SLAM     stop_codon  21269   21271    .        -       .     gene_id "002" ; transcript_id "002.1";
#  seq2      SLAM     CDS         21269   21444    .        -       2     gene_id "002" ; transcript_id "002.1"; exontype "terminal"
#  seq2      SLAM     CDS         21521   21656    .        -       0     gene_id "002" ; transcript_id "002.1"; exontype "internal"
#  seq2      SLAM     CDS         21938   22578    .        -       2     gene_id "002" ; transcript_id "002.1"; exontype "internal"
#  seq2      SLAM     CDS         22944   23925    .        -       0     gene_id "002" ; transcript_id "002.1"; exontype "initial"
#  seq2      SLAM     start_codon 23923   23925    .        -       .     gene_id "002" ; transcript_id "002.1";
#  seq2      SLAM     SLAM_CNS    26181   26242    .        .       .     cns_id "007"  ; identity "68.9"


# interesting lines
# <seqname> <source> <feature>   <start>  <end> <score> <strand> <frame> [attributes]   [comments]
#  seq1      SLAM     CDS         18350   18443    .        +       1     gene_id "001" ; transcript_id "001.1"; exontype "internal"
#  seq2      SLAM     CDS         18525   18594    .        +       0     gene_id "001" ; transcript_id "001.1"; exontype "internal"
#  seq2      SLAM     CDS         18731   18747    .        +       2     gene_id "001" ; transcript_id "001.1"; exontype "terminal"
#  seq2      SLAM     stop_codon  18745   18747    .        +       .     gene_id "001" ; transcript_id "001.1";
#  seq2      SLAM     stop_codon  21269   21271    .        -       .     gene_id "002" ; transcript_id "002.1";
#  seq2      SLAM     CDS         21269   21444    .        -       2     gene_id "002" ; transcript_id "002.1"; exontype "terminal"
#  seq2      SLAM     CDS         21521   21656    .        -       0     gene_id "002" ; transcript_id "002.1"; exontype "internal"
#  seq2      SLAM     CDS         21938   22578    .        -       2     gene_id "002" ; transcript_id "002.1"; exontype "internal"
#  seq2      SLAM     CDS         22944   23925    .        -       0     gene_id "002" ; transcript_id "002.1"; exontype "initial"
#  seq2      SLAM     start_codon 23923   23925    .        -       .     gene_id "002" ; transcript_id "002.1";


sub _parse_results {

  # parsing the results and adding predicted exons (predictionExon-obj) to transcript-objects
  my ($self,$slice,$basefile) = @_;

  my $path = $self->workdir;
  my $gff  = $self->workdir."/".$basefile.".gff";

  # what to do if resultfile contains no data ?
  $gff  = "/tmp/seq1.38047.gff";

  # we will walk through the file
  # if a line has a gene_id we will store the gene-id in a hash-keyhash
  # we also store the line in another hash
  # later, we'll process the two hashes
  my (%genes,%predex);

  open(IN,"$gff") || die "could not read file $gff";

  while (<IN>) {
    chomp;
    if (/gene_id/) {
      my @line = split /;/; #only the attributes 0-9 are stable
      my @attributes = /\s+/,$line[0];
      my $key = "$attributes[8] $attributes[9]";
      $genes{$key}=();  # key: -->gene_id "001"<-- FIND THE ERROR !
      $predex{$_}=\@attributes;
    }
  } # EOF
  close(IN);

## now let's process the predicted genes !!!
#  my ($gene,$predict);
#  foreach(keys %genes) {
#    $gene = $_;
#    print "\n-$gene-\n";
#    foreach(keys %predex){
#      if (/$gene/){
#        # reference to an array !!!
#        my @temp = @{$predex{$_}};
#        print "$temp[0]\t$temp[1]\n";
#      }
#    }
#  }
##    # processing the whole line

#      my ($seq,$src,$feat,$start,$end,$score,$strand,$frame,$gid_act,$ginr) = split/\s+/, $line[0];
#      my ($tr_id,$tr_nr) = split/\s+/, $line[1];
#      my ($et, $extype) = split/\s+/, $line[2];

#    if ($ln=~m/gene_id/) {      # we have coding sequence and we have a transcript

#      if ($geneid_old ne $geneid_act) {
#      $strand = 1 if ($strand eq "+");
#      $strand = -1 if ($strand eq "-");
#      $strand = 0 if ($strand eq ".");



#      my $ex = new Bio::EnsEMBL::PredictionExon(
#                                                -START     => $start,
#                                                -END       => $end,
#                                                -STRAND    => $strand, #valid values (Feature.pm: 1,-1,0)
#                                                -SLICE     => $slice,
#                                                -DBID      => 100, # EDIT!
#                                                -P_VALUE   => 0,
#                                                -SCORE     => 0
#                                               );
#    }
#    $geneid_old = $geneid_act;

}



=pod

=head2 run

  Title    : run
  Usage    : $obj->run
  Function : runs the slam-algorithm on the two supplied fasta-sequences. The run needs the modified approximate alignemnt 
             (see Bio::EnsEMBL::Pipeline::Tools::ApproxAlign ).
  Returns  : none
  Args     : none

=cut


sub run {
  my ($self) = @_;

  my $gcdir  = $self->_getgcdir;
  my $fasta1 = $self->fasta1;
  my $fasta2 = $self->fasta2;

  my $command =  $self->slam_bin .
    " -a ".$self->approx_align .
      " -p ".$gcdir . " ".$fasta1 . " " . $fasta2 .
        " -org1 ".$self->org1 .
          " -org2 ".$self->org2;

  $command .= " -v " if $self->verbose;
  $command .= " -debug " if $self->debug;

  print "slam-command; $command\n" if $self->verbose;

  open( EXO, "$command |" ) || $self->throw("Error running Slam $!");
  close(EXO);


  $fasta1=~s/(.+)\.(fasta|fa)/$1/; # get rid of suffix (.fasta or .fa)
  $fasta2=~s/(.+)\.(fasta|fa)/$1/; # get rid of suffix (.fasta or .fa)

  my $wdir = $self->workdir;

  $fasta1=~s/$wdir\///;         # get rid of workdir-prefix (/tmp/)
  $fasta2=~s/$wdir\///;         # get rid of workdir-prefix 



  $self->file($fasta1."_".$fasta2.".cns");

  $self->_parse_results($self->slice1,$fasta1); #  only file-basename
  $self->_parse_results($self->slice2,$fasta2); #  only file-basename
  return 1
}

sub files_to_delete {
  my ($self,$file) = @_;

  $file=~s/(.+)\.(fasta|fa)/$1/; # get rid of suffix (.fasta or .fa)
  my @suffix = qw (.gff .rna .pep );
  foreach (@suffix) {
    $self->file($file.$_);
  }
  return;
}




sub _getgcdir{
  my $self = shift;

  my $seq1 = $self->fasta1;
  my $seq2 = $self->fasta2;

  my $org1 = $self->org1;
  my $org2 = $self->org2;

  my $gcdirs = {
                'H.sapiens_M.musculus' => [
                                           [ 0,  43],
                                           [43,  51],
                                           [51,  57],
                                           [57, 100]
                                          ],
               },

                 my $pairName = sprintf("%s_%s",$org1,$org2);
  my $gcdir;

  if (exists($gcdirs->{$pairName})) {
    print "We have paramter bins defined for this organism pair.\n" if $self->verbose;

    my($seqStream,$seqObj);

    open(SEQ1,$seq1) || die "Can't open $seq1 for read: $!.\n";
    $seqStream = Bio::SeqIO->new(-fh => \*SEQ1, -format => 'Fasta');
    $seqObj = $seqStream->next_seq();
    my $gc1 = &_gccontent($seqObj);
    close SEQ1;

    open(SEQ2,$seq2) || die "Can't open $seq2 for read: $!.\n";
    $seqStream = Bio::SeqIO->new(-fh => \*SEQ2, -format => 'Fasta');
    $seqObj = $seqStream->next_seq();
    my $gc2 = &_gccontent($seqObj);
    close SEQ2;

    my $len1 = &_seqlen($seq1); ## HERE
    my $len2 = &_seqlen($seq2); ## HERE
    my $gc = (($gc1 * $len1) + ($gc2 * $len2)) / ($len1 + $len2);

    $gcdir = undef;
    my $nBins = scalar(@{$gcdirs->{$pairName}});
    for (my $i=0; $i<$nBins; $i++) {
      if (($gc > $gcdirs->{$pairName}[$i][0]) && ($gc <= $gcdirs->{$pairName}[$i][1])) {
        $gcdir = sprintf("bin%d",$i+1);
        last;
      }
    }
    die "Didn't find a GC dir for GC content $gc.\n" if(!defined($gcdir));
  } else {
    # No binning for this organism pair.
    $gcdir = "bin1";
  }
  my $pardir = sprintf("%s/%s_%s/%s",$self->slam_pars_dir,$org1,$org2,,$gcdir);
  return($pardir);
}


sub _gccontent{
  my($seq) = @_;

  my $tempSeq = $seq->seq();
  my $n = length($tempSeq);
  my(@nGC) = ($tempSeq =~ /[GC]/gi);

  return(100*scalar(@nGC)/$n);
}

sub _seqlen {
  my($seqFile) = @_;

  open(SEQFILE,$seqFile) || die "Can't open $seqFile for read\n";
  my $seqStr = Bio::SeqIO->new(-fh => \*SEQFILE, -format => 'Fasta' );
  my $seq = $seqStr->next_seq();
  my $len = length($seq->seq());
  close SEQFILE;

  return $len;
}




sub fasta1 {
  my $self = shift;

  $self->{_fasta1} = shift   if (@_);
  return $self->{_fasta1}
}


sub fasta2 {
  my $self = shift;

  $self->{_fasta2} = shift   if (@_);
  return $self->{_fasta2}
}


sub approx_align {
  my $self = shift;

  $self->{_approx_align} = shift   if (@_);
  return $self->{_approx_align}
}


sub org1{
  my ($self,$org1) = @_;

  $self->{_org1} = 'H.sapiens' if (!defined $org1 && !defined $self->{_org1});
  $self->{_org1} = $org1 if (defined $org1);
  return $self->{_org1};
}

sub org2{
  my ($self,$org2) = @_;

  $self->{_org2} = 'M.musculus' if (!defined $org2 && !defined $self->{_org2});
  $self->{_org2} = $org2 if (defined $org2);
  return $self->{_org2};
}


sub slam_bin {
  my ($self, $location) = @_;

  if (defined $location) {
    $self->{_slam_bin} = $location;
    $self->throw("Slam not found at $location: $!\n") unless (-e $location);
  }
  $self->{_slam_bin} = '/nfs/acari/jhv/bin/slam' if (!defined $location && !defined $self->{_slam_bin});

  return $self->{_slam_bin};
}


sub slam_pars_dir {
  my ($self, $pars_dir) = @_;

  if (defined $pars_dir) {
    $self->{_slam_pars_dir} = $pars_dir;
    $self->throw("Slam not found at $pars_dir: $!\n") unless (-e $pars_dir);
  }
  $self->{_slam_pars_dir} = '/nfs/acari/jhv/lib/slam_pars_dir' if (!defined $pars_dir && !defined $self->{_slam_pars_dir});
  return $self->{_slam_pars_dir};
}



sub minlength{
  my ($self,$minlength) =@_;

  $self->{_minlength} = 250 if (!defined $minlength && !defined $self->{_minlength});
  $self->{_minlength} = $minlength if (defined $minlength);
  return $self->{_minlength};
}

sub slice1{
  my ($self,$slice1) = @_;

  $self->{_slice1} = $slice1 if (defined $slice1);
  $self->throw("Slam needs a first slice to work on!\n") if (! defined $self->{_slice1} && !defined $slice1);
  return $self->{_slice1};
}

sub slice2{
  my ($self,$slice2) = @_;

  $self->{_slice2} = $slice2 if (defined $slice2);
  $self->throw("Slam needs a first slice to work on!\n") if (! defined $self->{_slice2} && !defined $slice2);
  return $self->{_slice2};
}




############################################################

sub debug{
  my ($self,$debug) = @_;

  $self->{_debug} = '0' if (!defined $debug && !defined $self->{_debug});
  $self->{_debug} = $debug if (defined $debug);
  return $self->{_debug};
}


sub verbose{
  my ($self,$verbose) = @_;

  $self->{_verbose} = '0' if (!defined $verbose && !defined $self->{_verbose});
  $self->{_verbose} = $verbose if (defined $verbose);
  return $self->{_verbose};
}




sub printvars {
  my $self = shift;
  my %tmp = %{$self};

  print "Key:\t\t\tValue:\n";
  print "------------------------------\n";
  foreach (keys %tmp) {
    print "-$_-\t\t\t-$tmp{$_}-\n";
  }
  print "------------------------------\n";
}

1;

