#
# Written by Jan-Hinnerk Vogel
#
# jhv [at] sanger.ac.uk
#
# Copyright GRL/EBI 2004
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=pod

=head1 NAME

 Bio::EnsEMBL::Pipeline::Runnable::Avid

=head1 SYNOPSIS

 To construct a new Bio::EnsEMBL::Pipeline::Runnable::Avid-Object use

 $obj = Bio::EnsEMBL::Pipeline::Runnable::Avid->new(
                                                      -slice1        => $slice1,
                                                      -slice2        => $slice2,
                                                      -avid          => '/path/to/avid',
                                                      -avid_options  => '-obin' or '-opsl',
                                                      -workdir       => 'path/to/workdir',
                                                      -slam_output   => 'true' or 'false'
                                                     );

or

 $obj = Bio::EnsEMBL::Pipeline::Runnable::Avid->new(
                                                      -slice1        => $slice1,
                                                      -slice2        => $slice2,
                                                     );


 $obj -> run;

 $result = $obj -> results;
 $fasta1 = $obj -> fasta_filename1;
 $fasta2 = $obj -> fasta_filename2;

=head1 DESCRIPTION

Avid takes two Bio::Seq (or Bio::PrimarySeq) objects and aligns them
against each other using their FASTA-sequences and their repeatmasked
sequences. If only two slices are passed to Avid, the default options 
are choosen (means: workdir = '/tmp/', avid_option = '-obin', -slam-output='true'
and avid = '/nfs/acari/jhv/bin/avid').


=head1 OPTIONS

=over

=item  B<-slice1>         reference to the first an Bio::EnsEMBL::Slice - object

=item  B<-slice2>         reference to the second an Bio::EnsEMBL::Slice - object

=item  B<-avid>           optional path to the used avid-binary

=item  B<-avid_options>   format of avid-output (-obin || -opsl, default:-obin)

=item  B<-workdir>        optional working-directory (default: /tmp/)

=item  B<-slam_output>    parse binary file and write aat-file (true or false)





=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 METHODS

The rest of the documentation details each of the object methods.
Internal (private) methods are usually preceded with a "_".

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Avid;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);





sub new {

  my ($class,@args) = @_;
  my $self = {};
  bless $self,$class;

  my (
      $slice1,
      $slice2,
      $avid,
      $avid_options,
      $workdir,
      $slam_output ) = $self->_rearrange([qw(
                                             SLICE1
                                             SLICE2
                                             AVID
                                             AVID_OPTIONS
                                             WORKDIR
                                             SLAM_OUTPUT
                                            )
                                         ], @args);
  $self->{_verbose} = 0 ;

  # location of avid-binary avid-2.1b0
  $self->_avid_binary($avid);

  # storing slices
  $self->_slice1 ( $slice1 );
  $self->_slice2 ( $slice2 );

  # output-option for avid (def: binary output -obin)
  $self->_avid_options( $avid_options );

  # set names of fasta-files
  $self->_filename1;
  $self->_filename2;

  # should the binary-file be parsed for a slam-run?
  $self->_slam_output_opt($slam_output);

  $self->workdir($workdir);
  $self->checkdir;

  return $self;
}


sub DESTROY {
  my $self = shift;
  #  $self->deletefiles;
}



=pod

=head2 run

  Title    : run
  Usage    : $obj->run
  Function : aligns the masked and unmasked sequences of the two given Bio::PrimarySeqI-Objects.
  Returns  : nothing
  Args     : nothing

=cut

sub run {
  my ($self) = @_;

  # Write masked & unmasked query sequences (4 files)
  $self->_write_sequences;

  # fasta filenames which avid needs (there MUST also be 2 .masked-files !)
  my $fa_first = $self->workdir."/".$self->_filename1.".fasta";
  my $fa_secnd = $self->workdir."/".$self->_filename2.".fasta";

  my $command = $self->_avid_binary  . " " .$self->_avid_options;
  $command .=  " " . $fa_first    . " ";
  $command .=  " " . $fa_secnd         ;

  print STDERR "avid-command : $command\n"  if $self->_verbose;

  open( AVID, "$command |" ) || $self->throw("Error running avid $!");
  close(AVID);

  # parse binary output for slam-run
  $self->_parse_binary if ($self->_slam_output_opt);

  # register files written by avid (mout,psl,...)
  $self->_avid_files;

  return 1
}


############################################################

sub _avid_files {
  my ($self) = shift;

  my $workdir   = $self->workdir;
  my $filename1 = $self->_filename1;
  my $filename2 = $self->_filename2;

  # option -opsl writes out: .minfo .info .out .psl
  # option -obin writes out: .minfo .mout

  my $basename = $workdir."/".$filename1."_".$filename2;

  $self->file   ( $basename . ".minfo"); #adding to list of files to delete

  if ($self->_avid_options =~ /opsl/) {
    $self->results( $basename  . ".psl");
    $self->file( $basename . ".psl");
    $self->file( $basename . ".out");
    $self->file( $basename . ".info");
  } else {
    $self->results ( $basename . ".mout" );
    $self->file ( $basename . ".mout" );
  }
  return;
}



sub _slice1{
  my ($self,$slice1) = @_;

  if (defined $slice1) {
    unless ($slice1->isa("Bio::PrimarySeqI")) {
      $self->throw("Submitted Sequence is not a Bio::PrimarySeqI");
    }
    $self->{_slice1} = $slice1;
  }

  return $self->{_slice1};
}


sub _slice2{
  my ($self,$slice2) = @_;

  if (defined $slice2) {
    unless ($slice2->isa("Bio::PrimarySeqI")) {
      $self->throw("Submitted Sequence is not a Bio::PrimarySeqI");
    }
    $self->{_slice2} = $slice2;
  }

  return $self->{_slice2};
}


sub _avid_binary {
  my ($self, $location) = @_;

  if (defined $location) {
    print "Location of Avid-File is $location\n" if ($self->_verbose);
    $self->throw("Avid not found at $location: $!\n") unless (-e $location);
    $self->{_avid} = $location ;
  }
  $self->{_avid}='/nfs/acari/jhv/bin/avid' if ( !defined $location  && !defined $self->{_avid});

  return $self->{_avid};
}


sub _avid_options {
  my ($self,$opt) =  @_;

  $self->{_avid_options} = $opt if (defined $opt);
  $self->{_avid_options} = '-obin' if ( ( !defined $opt) && ( !defined  $self->{_avid_options}));

  return $self->{_avid_options};
}


sub _slam_output_opt {
  my ($self,$out) =  @_;

  if (defined $out) {
    if ($out eq "1" || $out=~m/(true|t)/i) {
      $self->{_slam_output} = '1';
    } else {
      $self->{_slam_output} = '0' if ( ( !defined $out) && ( !defined  $self->{_slam_output}));
    }
  }
  return $self->{_slam_output};
}



############################################################

sub _write_sequences {
  my $self = shift;

  $self->_write_unmasked_sequence($self->_slice1,$self->_filename1);
  $self->_write_unmasked_sequence($self->_slice2,$self->_filename2);

  $self->_write_masked_sequence($self->_slice1,$self->_filename1);
  $self->_write_masked_sequence($self->_slice2,$self->_filename2);
  return;
}


# gets a slice-object and writes unmaked sequence

sub _write_unmasked_sequence{
  my ($self,$slice,$nam) = @_;

  my $file = $self->workdir."/".$nam.".fasta";

  # writing unmasked sequence
  my  $seqobj = Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');

  $seqobj->write_seq($slice);
  print "Writing sequence $file\n" if ($self->_verbose);

  $self->file($file);
}


# writing masked fasta-sequence
sub _write_masked_sequence{
  my ($self,$slice,$nam) = @_;

  my $file = $self->workdir."/".$nam.".fasta.masked";

  my $rm_slice = $slice->get_repeatmasked_seq();
  my  $seqobj =  Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
  $seqobj -> write_seq($rm_slice);

  print "Writing sequence $file\n" if ($self->_verbose);

  $self->file($file);
  return;
}


sub _filename1{
  my ($self) = shift;

  if (! defined $self->{_filename1} ) {
    my $tmp = $self->get_tmp_file("","seq1","");
    $tmp =~s/^\/(.+)\.$/$1/;    #snipping away first slash and last dot
    $self->{_filename1} = $tmp;
  }
  return $self->{_filename1};
}

sub _filename2{
  my ($self) = shift;;

  if (! defined $self->{_filename2} ) {
    my $tmp = $self->get_tmp_file("","seq2","");
    $tmp =~s/^\/(.+)\.$/$1/;    #snipping away first slash and last dot
    $self->{_filename2} = $tmp;
  }
  return $self->{_filename2};
}


=pod

=head2 fasta_filename1

  Title    : fasta_filename1
  Usage    : $obj->results
  Function : returns the path+name of first fasta-sequence
  Returns  : String
  Args     : nothing

=cut

sub fasta_filename1{
  my ($self) = shift;

  $self->{fasta_filename1}=$self->workdir."/".$self->_filename1.".fasta";
  return $self->{fasta_filename1};
}


=pod

=head2 fasta_filename2

  Title    : fasta_filename2
  Usage    : $obj->results
  Function : returns the path+name of second fasta-sequence
  Returns  : String
  Args     : nothing

=cut

sub fasta_filename2{
  my ($self) = shift;

  $self->{fasta_filename2}=$self->workdir."/".$self->_filename2.".fasta";
  return $self->{fasta_filename2};
}



sub _parse_binary {
  my ($self) =  @_;

  $self->checkdir;
  my $wdir = $self->workdir;

  # constructing name of written binary .mout-file seq1.2234_seq2.8952.mout


  my $filename1 = $self->_filename1;
  my $filename2 = $self->_filename2;


  my $basefile = $self->workdir."/".$filename1."_".$filename2;
  my $binfile = $basefile.".mout";
  my $outfile = $basefile."_parsed_binary.txt";

  $self->parsed_binary_filename($outfile);

  ####### START PARSING BINARY OUTPUT OF AVID-ALIGNEMENT #######

  my $lInd = -1;
  my $rInd = -1;
  my $byte;
  my @x;

  open(FH,"<$binfile") || $self->throw("Could not open Avid's binary output-file: $binfile\n");
  while (read(FH,$byte,1)) {
    my $intByte = unpack("C",$byte);
    my($lHalf,$rHalf) = ( ($intByte >> 4), ($intByte & 0xF) );
    $lInd++ if($lHalf>0);
    $rInd++ if($rHalf>0);
    $x[$lInd] = $rInd if($lHalf > 0);
  }
  close(FH);

  my $lLen = $lInd+1;
  my $rLen = $rInd+1;

  my $halfWin = 1;
  my($i,$min,$max);
  my $lastMax = -1;
  open(OUT,">$outfile") || $self->throw("Could not write to output-file $outfile\n");
  for ($i=0; $i < $lLen-1; $i++) {
    $min = $x[$i] - $halfWin;
    $min = 0 if($min < 0);
    $min = $lastMax+1 if($min > ($lastMax+1));
    $max = $x[$i] + $halfWin;
    $max = $rLen-1 if($max >= $rLen);
    print OUT  join("\t",$i,$min,$max) . "\n";
    $lastMax = $max;
  }
  $min = $x[$i] - $halfWin;
  $min = 0 if($min < 0);
  $min = $lastMax+1 if($min > ($lastMax+1));
  $max = $rLen-1;

  print OUT join("\t",$i,$min,$max) . "\n";
  close(OUT);

  ####### STOP PARSING BINARY OUTPUT OF AVID-ALIGNEMENT #######

  $self->file($outfile);

  return $outfile;
}


sub _verbose {
  my $self = shift;

  $self->{_verbose} = shift if @_ ;
  return $self->{_verbose}
}

=pod

=head2 parsed_binary_filename

  Title    : parsed_binary_filename
  Usage    : $obj->parsed_binary_filename
  Function : returns the path+name of the parsed binary output (if option -slam_output => true)
  Returns  : String
  Args     : String

=cut

sub parsed_binary_filename{
  my ($self,$fn) = @_;
  $self->{_parsed_binary_filename}=$fn if (defined $fn && !defined $self->{_parsed_binary_filename});
  return $self->{_parsed_binary_filename};
}

=pod

=head2 results

  Title    : results
  Usage    : $obj->results
  Function : returns the path+name of the written alignement-file (binary .mout-file or .psf-file)
  Returns  : String
  Args     : nothing

=cut

1;

