#
# Written by Jan-Hinnerk Vogel
#
# jhv [at] sanger.ac.uk
#
# Copyright GRL/EBI 2004 a
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
                                                      -max_memory    => 1572864,
                                                      -minlength     => 250,
                                                      -debug         => 0,
                                                      -verbose       => 0
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

=item  B<-max_memory>    maximum memory size of slam-process (default 1.5 Gb)

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
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::PredictionExon;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use diagnostics;




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
      $max_memory,              # max. memory the slam-process is allowed to use, default 1.5 Gb
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
                               MAX_MEMORY
                               MINLENGTH
                               WORKDIR
                               DEBUG
                               VERBOSE
                              )
                           ], @args
                          );


  # setting defaults

  $self->max_memory_size(1572864);


  $self->slices($slice1,$slice2);
  $self->fasta($fasta1,$fasta2);

  $self->approx_align($approx_align);

  $self->slam_bin($slam_bin);
  $self->slam_pars_dir($slam_pars_dir);
  $self->max_memory_size($max_memory);
  $self->minlength($minlength);

  $self->workdir($workdir);
  $self->debug($debug);

  $self->verbose($verbose);
  $self->verbose("1");

  $self->printvars if ($self->verbose);

  return $self;
}



sub DESTROY {
  my $self = shift;
     $self->deletefiles;
}

=pod

=head2 _parse_results ()

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


sub parse_results {
  my $self = shift;

  # parsing results of first organism
  # returns ref to array of predicted transcripts
  # could be empty arrays or arrays full of predicted transcripts (paradise :-)

  my $arrayref1 = $self->_parser( ${$self->slices}[0], ${$self->fasta}[0] );
  my $arrayref2 = $self->_parser( ${$self->slices}[1], ${$self->fasta}[1] );

  #     [HPT HPT HPT] [HM HM HM] or  [ [][] ]
  $self->predtrans( $arrayref1,$arrayref2 );
}


# parsing the results and adding predicted exons to container of PredictionTranscripts
# only the coding sequences are stored (CDS).
# returns an array-reference to a PredictionTranscript with stored PredictionExons

sub _parser {
  my ($self,$slice,$gff) = @_;
  my (%transcripts);
  print "gff: $gff\n" if $self->verbose;
  $gff=~s/(\.fasta)/\.gff/;     # subst. .fasta-suffix with .gff-suffix

  ## data for error-message if slam-run fails due to out-of-memory-error
  my $e_start = $slice->start;
  my $e_end = $slice->end;
  my $e_chr = $slice->seq_region_name;


  # processing the written gff-file
    open(IN,"$gff") || $self->throw("Slam.pm: OUT OF MEMORY-ERROR !\n XXXX\nCould not read $gff: $0 \n");
    while (<IN>) {
      chomp;
      my $aline = $_;


#### NEW REDESIGN WITH START/STOPCODONS

      if ($aline =~m/gene_id/) {
        # line contains "CDS"

        if ( $aline !~/start_codon/  &&  $aline !~/stop_codon/) {
         # if the line is no start-or stopcodon
          my @line = split /;/;   # split line in 3 parts, bcs only the attributes 0-9 are stable
          my @attributes = split /\s+/,$line[0];
          my @transc_id =  split /\s+/,$line[1];

          # building a unique key for each transcript
          my $key = "$attributes[8] $attributes[9]"; # key: -->gene_id "001"<--

          # we concatenate the first and second part of the line
          push (@attributes, @transc_id);
          my $attr_ref = \@attributes;        

          # we store the ref of line with attributes of the exon using the unique transcript-key as key
          # an HASH of ARRAYS of ARRAYS
          push (@{ $transcripts{$key}}, $attr_ref );

	  # {gene_id "001"}=[ "chr1 SLAM CDS 19053348 19053557 . - 1 gene_id "M4H1U1D4-28.002" transcript_id bla",
	  #                   "chr1 SLAM CDS 19053348 19053557 . - 1 gene_id "M4H1U1D4-28.002" transcript_id bla"....
	  #   	             ....
	  #		    ]

        }
	}









##### OLD PART ######

#      if ($aline =~m/gene_id/) {
#        #if the line contains a CDS
#        if ( $aline !~/start_codon/  &&  $aline !~/stop_codon/) {
#          # if the line is no start-or stopcodon
#
#          my @line = split /;/;   # split line in 3 parts, bcs only the attributes 0-9 are stable
#          my @attributes = split /\s+/,$line[0];
#          my @transc_id =  split /\s+/,$line[1];
#
#          # building a unique key for each transcript
#          my $key = "$attributes[8] $attributes[9]"; # key: -->gene_id "001"<--
#
#          # we concatenate the first and second part of the line
#          push (@attributes, @transc_id);
#          my $attr_ref = \@attributes;
#
#          # we store the ref of line with attributes of the exon using the unique transcript-key as key
#          # an HASH of ARRAYS of ARRAYS
#          push (@{ $transcripts{$key}}, $attr_ref );
#        }
#      }


    }
    close(IN);

  # array of length 0
  my @all_predicted_transcripts=();
  my $pred_trans;

  # for every predicted transcript
  for my $transkey (sort(keys %transcripts)) {

    $pred_trans = new Bio::EnsEMBL::PredictionTranscript();

    # get all the exons of the predicted transcript
    my @exon_refs = @{$transcripts{$transkey}};

    # process every predicted exon which belongs to the transcript
    for my $exons (@exon_refs) {
      my  @attributes = @{$exons}; # @attributes contains all items of a gff-fileline (with CDS)

      # extracting startbp, endbp strand and phase
      my ($start,$end,$strand,$phase) = @attributes[3,4,6,7];

      $strand = 1 if ($strand eq "+");
      $strand = -1 if ($strand eq "-");
      $strand = 0 if ($strand eq ".");

      my $pred_exon = new Bio::EnsEMBL::PredictionExon(
                                                       -START     => $start,
                                                       -END       => $end,
                                                       -STRAND    => $strand, # valid values (Feature.pm: 1,-1,0)
                                                       -SLICE     => $slice,
                                                       -P_VALUE   => 0, #try undef here
                                                       -PHASE     => $phase , # same as frame ?
                                                       -SCORE     => 0 #try undef here
                                                      );
      # add predicted exon to the transcript-container (Bio::EnsEMBL::Feature)
      $pred_trans->add_Exon( $pred_exon );
    }
    push @all_predicted_transcripts, $pred_trans; # [HPT HPT HPT]
  }
  return \@all_predicted_transcripts; # retrun \[HPT] or \[] emtpy array
}



=pod

=head2 run

  Title    : run
  Usage    : $obj->run
  Function : runs the slam-algorithm on the two supplied fasta-sequences and reads the modified approximate alignemnt-file.
             (see Bio::EnsEMBL::Pipeline::Tools::ApproxAlign ). The results (gff-files) will be parsed.
  Returns  : none
  Args     : none

=cut


sub run {
  my ($self) = @_;

  my $gcdir  = $self->_getgcdir;
  my $fasta1 = ${$self->fasta}[0];
  my $fasta2 = ${$self->fasta}[1];
  my $slice1 = ${$self->slices}[0];
  my $slice2 = ${$self->slices}[1];

  my $t = ${$self->org}[1];
  my $command =  $self->slam_bin .
    " -a ".$self->approx_align .
        " -p ".$gcdir . " ".$fasta1 . " " . $fasta2 .
          " -org1 ".${$self->org}[0] .
            " -org2 ".${$self->org}[1];
#  $command .= " -v " if $self->verbose;
#  $command .= " -debug " if $self->debug;

  print "Slam.pm: Slam-command:\n $command\n" if $self->verbose;

  # kick of fasta suffix for gff-filenames and cns-tempfile
  (my $base1 = $fasta1) =~s/(.+)\.(fasta|fa)/$1/; # get rid of suffix (.fasta or .fa)
  (my $base2 = $fasta2) =~s/(.+)\.(fasta|fa)/$1/; # get rid of suffix (.fasta or .fa)

  my $gff1 = $base1.".gff";
  my $gff2 = $base2.".gff";

  my $slength1 = &_seqlen($fasta1);
  my $slength2 = &_seqlen($fasta2);

  print "Slam.pm: Length of Sequence1: $slength1\t Sequenc2: $slength2\tSlam-run follows\n" if $self->verbose;

  if (($slength1<$self->minlength) || ($slength2 < $self->minlength)) {

    $self->warn("Sequencelength of cutted sequence is smaller than min. Offset  (Slamconf.pm). The Sequence will be skipped (no Slam-run)");
    # if one of the sequences is too short to process, just make empty gff files.
    open(GFF, ">$gff1") ||  die "$0: $gff1: $!\n";
    close (GFF);
    open(GFF, ">$gff2") ||  die "$0: $gff2: $!\n";
    close (GFF);
  } else {
    print "Slam.pm: try to evalute / run Slam\n" if $self->verbose;
    my $pstat;
    my $max_memory = $self->max_memory_size;

    print "Slam.pm: max-memory-size: $max_memory\n" if $self->verbose;

    eval {  $pstat = system ("ulimit -v  $max_memory  -c 0 ; $command")    };

#    unless ($pstat ne "0") {
#      print "WARNING\t SLAM RUN NEEDED LOTS OF MEM\t BROKE UP\n" ;
#      print "Slam.pm: Slam run seems to be buggy\n";
#    }
  }
    # add cns-file to list of tempfiles
    my $wdir = $self->workdir;
    $base1=~s/$wdir\///;       # get rid of workdir-prefix (/tmp/)
    $base2=~s/$wdir\///;       # get rid of workdir-prefix 
    $self->file($base1."_".$base2.".cns");

  # parse results if resultfiles are written
    if ( (-e $gff1) &&( -e $gff2 )) {
      $self->parse_results;  # parse existing gff's
  }else{
    # error-message, evtl re-analysis
    my $e1_start = $slice1->start;    my $e1_end = $slice1->end;    my $e1_chr = $slice1->seq_region_name;
    my $e2_start = $slice2->start;    my $e2_end = $slice2->end;    my $e2_chr = $slice2->seq_region_name;

    print "\t\t\t\t------------------\n\t\t\t\tOUT OF MEMORY:\n";
    print "\tFailed to run Slam on region $e1_chr-$e1_start-$e1_end---$e2_chr-$e2_start-$e2_end\n\t\t\t\t------------------\n";
    # store empty array because the run failed
    my $aref = [];
    $self->predtrans($aref,$aref);
  }

  $self->files_to_delete(${$self->fasta}[0]);
  $self->files_to_delete(${$self->fasta}[1]);
  return 1;
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

  my $seq1 = ${$self->fasta}[0];
  my $seq2 = ${$self->fasta}[1];

  my ($org1,$org2) = @{$self->org};


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

    my $len1 = &_seqlen($seq1); ## CHANGE HERE
    my $len2 = &_seqlen($seq2); ## CHANGE HERE
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

  open(SEQFILE,$seqFile) || die "Slam.pm: Can't open $seqFile for read\n";
  my $seqStr = Bio::SeqIO->new(-fh => \*SEQFILE, -format => 'Fasta' );
  my $seq = $seqStr->next_seq();
  my $len = length($seq->seq());
  close SEQFILE;
  return $len;
}


=pod

=head2 fasta ( filename )

  Title    : fasta
  Usage    : $obj->fasta
  Function : sets/gets the path and name of the fasta-files
  Returns  : Array-reference
  Args     : two strings containing /path/to/fata-files

=cut


sub fasta {
  my ($self,$fasta1,$fasta2) = @_;

  $self->{_fasta} = [$fasta1,$fasta2] if ($fasta1 && $fasta2);
  return $self->{_fasta};
}


=head2 predtrans (ref1,ref2)

  Title    : predtrans
  Usage    : $obj->predtrans
  Function : Sets/gets the Predicted transcripts for both organisms
  Returns  : Ref. to an Array of Arrayrefs. to Bio::EnsEMBL::PredictionTranscript
  Args     : References to two arrays

=cut

sub predtrans{
  my ($self,$ref_predtrans1,$ref_predtrans2) = @_;

  if ($ref_predtrans1 && $ref_predtrans2) {
    # Storing the refernces to the array of predicted transcripts in an array of predicted transcripts
    $self->{_ref_predtrans} = [$ref_predtrans1,$ref_predtrans2];
  }
  return $self->{_ref_predtrans};
}


=head2 max_memory_size (int)

  Title    : max_memory_size
  Usage    : $obj->max_memory_size( value )
  Function : Sets/gets the maximum memory size which the slam process is allowed to allocate
  Returns  : integer
  Args     : integer

=cut


sub max_memory_size {
  my ($self,$mem) = @_;

  $self->{_max_memory} = $mem if defined $mem;
  return $self->{_max_memory}
}


sub approx_align {
  my $self = shift;

  $self->{_approx_align} = shift   if (@_);
  return $self->{_approx_align}
}



sub org{
  my ($self,$org1,$org2) = @_;

  if (!defined $self->{_org}) {
    $org1 = 'H.sapiens'  if (!defined $org1);
    $org2 = 'M.musculus' if (!defined $org2);
    $self->{_org} = [$org1,$org2];
  }
  return $self->{_org};
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



sub slices{
  my ($self,$slice1,$slice2) = @_;

  if (defined $slice2 && defined $slice1) {
    $self->{_slices} = [$slice1,$slice2];
  }

  return $self->{_slices};
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

