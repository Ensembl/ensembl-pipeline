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
                                                      -avid          => $avid,
                                                      -avid_options  => $avid_options,
                                                      -workdir       => $workdir,
                                                      -file_prefix   => \@file_prefix,
                                                      -slam_output   => $slam_output
                                                     );

or

 $obj = Bio::EnsEMBL::Pipeline::Runnable::Avid->new(
                                                      -slice1        => $slice1,
                                                      -slice2        => $slice2,
                                                     );


 $obj -> run;

 $result = $obj -> results;

=head1 DESCRIPTION

Avid takes two Bio::Seq (or Bio::PrimarySeq) objects and aligns them
against each other using their FASTA-sequences and their repeatmasked
sequences. The setting of a working directory and a prefix for the written
.fasta and .fasta.masked files is optional.

=head1 OPTIONS

=over

=item  B<-slice1>         reference to the first an Bio::EnsEMBL::Slice - object

=item  B<-slice2>         reference to the second an Bio::EnsEMBL::Slice - object

=item  B<-avid>           optional path to the used avid-binary (default:)

=item  B<-avid_options>   options which are passed to avid

=item  B<-workdir>        optional working-directory (default: /tmp/)

=item  B<-file_prefix>    reference to an array of strings which contain the

                          prefix for written temp-files 

                          (default: ("seq1","seq2")


$slam_output  = an optional string to convert the binary output to a text-file



$obj->run;
my $resultfile = $obj->results;



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
      $file_prefix,
      $slam_output ) = $self->_rearrange([qw(
                                             SLICE1
                                             SLICE2
                                             AVID
                                             AVID_OPTIONS
                                             WORKDIR
                                             FILE_PREFIX
                                             SLAM_OUTPUT
                                            )
                                         ], @args);
  $self->{_verbose} = 1 ;

  ########### STORE NAME OF WRITTEN FASTA'S #############
  $self->{_fasta_filenames}=[];


  ########### TESTING TYPES OF PASSED OBJECTS ###########
  $self->query_sequences($slice1,$slice2);


  ########### LOCATION OF BINARY default av9d-2.1b0 #####
  $self->avid_binary($avid);


  ########### PASSED AVID_OPTIONS #######################
  $self->avid_options( $avid_options );


  ########## FILE_PREFIX ################################
  $self->file_prefix ($file_prefix);


  ########## WRITING SLAM OUTPUT ########################
  defined $slam_output ? $self->slam_output_opt($slam_output) : $self->slam_output_opt("False");
  print "slam putput ooption is set to:".$self->slam_output_opt."\n";

  ########## SETTING WORKDIR ############################
  $self->workdir($workdir);


  $self->checkdir;
  return $self;
}


sub DESTROY {
  my $self = shift;
  $self->deletefiles;
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

  # Write masked and unmasked query sequences to 4 files and store their names
  $self->write_sequences;


  # names of fasta files on which avid runs on (there must also be a .masked-file)

  my $fa_first = $self->workdir.@{$self->file_basename}[0];
  my $fa_secnd = $self->workdir.@{$self->file_basename}[1];

  $fa_first .="fasta";
  $fa_secnd .="fasta";

  my $command = $self->avid_binary  . " " .$self->avid_options;
  $command .=  " " . $fa_first    . " ";
  $command .=  " " . $fa_secnd         ;

  print STDERR "Avid command : $command\n"  if $self->verbose;

  open( AVID, "$command |" ) || $self->throw("Error running avid $!");
  close(AVID);

  # write output for slam 
  $self->parse_binary if ($self->slam_output_opt);

  # pass written minfo mout.. files for deletion to $self->file
#  $self->files_to_delete($fa_first,$fa_secnd); #313 

  $self->printvars;

  return 1
}


############################################################

sub files_to_delete {
  my ($self,$ff,$sf) = @_;

  my $workdir        =  $self->workdir;

  # substitution of workdir-prefix and .fasta-suffix to get filenames written by avid
  $ff             =~s/$workdir\/(.+)\.fasta/$1/;
  $sf             =~s/$workdir\/(.+)\.fasta/$1/;

  my $basename = $ff."_".$sf;


  # option -opsl writes out: .minfo .info .out .psl
  # option -obin writes out: .minfo .mout

  $self->file   ( $basename . ".minfo"); #adding to list of files to delete

  if ($self->avid_options =~ /obin/) {
    $self->results ( $basename . ".mout" );
    $self->file ( $basename . ".mout" );
    $self->file ( $basename . ".minfo" );
  }

  if ($self->avid_options =~ /opsl/) {
    $self->results( $basename  . ".psl"  );
    $self->file   ( $basename . ".info" );
    $self->file   ( $basename . ".out"  );
  }
}



sub query_sequences {
  my ($self,$slice1,$slice2) = @_;

  if (defined $slice1 && defined $slice2 ) {
    unless ($slice1->isa("Bio::PrimarySeqI") || $slice2->isa("Bio::PrimarySeqI")) {
      $self->throw("Submitted Sequence is not a Bio::PrimarySeqI");
    }
    my  @seqs;
    push @seqs, $slice1, $slice2;
    $self->{_query_sequences} = \@seqs;
  }
  return $self->{_query_sequences};
}

############################################################


sub avid_binary {
  my ($self, $location) = @_;

  if (defined $location) {
    print "Location of Avid-File is $location\n" if ($self->verbose);
    $self->throw("Avid not found at $location: $!\n") unless (-e $location);
    $self->{_avid} = $location ;
  }
  $self->{_avid}='/nfs/acari/jhv/bin/avid' if ( !defined $location  && !defined $self->{_avid});

  return $self->{_avid};
}

############################################################

sub avid_options {
  my ($self,$opt) =  @_;

  $self->{_avid_options} = $opt if (defined $opt);
  $self->{_avid_options} = '-obin' if ( ( !defined $opt) && ( !defined  $self->{_avid_options}));

  return $self->{_avid_options};
}

############################################################

sub slam_output_opt {
  my ($self,$out) =  @_;

  $self->{_slam_output} = "1" if (defined $out);
  $self->{_slam_output} = 'FALSE' if ( ( !defined $out) && ( !defined  $self->{_slam_output}));

  return $self->{_slam_output};
}

############################################################
# writing out 4 files for avid-input : 2 fasta-files and 2 repeatmasked fastas

sub write_sequences {
  my $self = shift;

  my $count = 0;
  foreach my $slice ( @{$self->query_sequences} ) {

    my $tmp_filename=$self->get_tmp_file("",@{$self->file_prefix}[$count],"");

    # push basefilename in an array of basenames
    $self->file_basename($tmp_filename);

    my $file = $self->workdir.$tmp_filename;

    # store names of fasta-files
    $self->fasta_filenames($file);

    # add .fasta
    $file .= "fasta";

    # writing unmasked sequence
    my  $seqobj = Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
    $seqobj->write_seq($slice);
    print "Writing sequence $file\n" if ($self->verbose);


    # inherited method to store filenames to unlink in @array
#    $self->file($file); #313

    # getting filename of (first) masked sequence 
    $file .= ".masked";

    # writing masked fasta-sequence
    my $rm_slice = $slice->get_repeatmasked_seq();
    $seqobj =  Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
    $seqobj -> write_seq($rm_slice);

    print "Writing sequence $file\n" if ($self->verbose);

    # unlink masked fasta
#    $self->file($file); #313
    $count++;
  }
}

############################################################

sub parse_binary {
  my ($self) =  @_;

  $self->checkdir;
  my $wdir = $self->workdir;

  #getting name of binary slam outputfile
  my $ff = @{$self->file_basename}[0];
  $ff =~s/^\/(.+)\.$/$1/;                #snipping away first slash and last dot

  my $sf = @{$self->file_basename}[1];
  $sf =~s/^\/(.+)\.$/$1/;                #snipping away last dot

  my $binfile = $self->workdir."/".$ff."_".$sf.".mout";
  my $outfile = $self->workdir."/".$ff."_".$sf."_slam_input.txt";

  $self->parsed_binary_output_filename($outfile);

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
 $self->file($outfile);   #uncomment to delete the aat-outputfile for slam txt outfiles for slam !!!
  return $outfile;
}


sub verbose {
  my $self = shift;

  $self->{_verbose} = shift if @_ ;
  return $self->{_verbose}
}


sub file_prefix{
  my ($self,$prefix) = @_;

  $self->{_file_prefix} = \@$prefix   if (defined $prefix && ref($prefix) eq "ARRAY");
  $self->{_file_prefix} = ["seq1","seq2"]   if (!defined $prefix && !defined $self->{_file_prefix});

  return $self->{_file_prefix};
}

sub parsed_binary_output_filename{
  my ($self,$fn) = @_;
  $self->{_parsed_binary_output_filename}=$fn if (defined $fn && !defined $self->{_parsed_binary_output_filename});
  return $self->{_parsed_binary_output_filename};
}

sub file_basename {
  my ($self,$base) = @_;

  push @{$self->{_file_basename}},$base if (!defined $self->{file_basename} && defined $base);
  $self->{_file_basename} = [$base] if (defined $base && !defined $self->{_file_basename});
  return $self->{_file_basename}; # /tmp/seq1.88796.
}


sub fasta_filenames {
  my ($self,$name) = @_;

  # if a filename is supplied, add it to the array
  push @{$self->{_fasta_filenames}},$name."fasta" if ( defined $name);
  return $self->{_fasta_filenames};
}

sub printvars {
  my $self = shift;
  my %tmp = %{$self};

  print "Key:\t\t\tValue:\n";
  print "------------------------------\n";
  foreach (keys %tmp ) {
    print "-$_-\t\t\t-$tmp{$_}-\n";
  }
}

1;

