#
# Written by Jan-Hinnerk Vogel
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
@sequences   = an array Bio::Seq objects,
$avid        = an optional string which specifies the location of the used avid binary
$options     = an optional string with  options (default: binary output (.mout-file))
$workdir     = an optional string containing the location of the working-directory (def: /tmp:)
$file_prefix = an optional array of strings used as file_prefixes for the temporarily written output-files

my $obj = Bio::EnsEMBL::Pipeline::Runnable::Avid->new(
                                                      -query_seqs    => \@sequences,
                                                      -avid          => $avid,
                                                      -options       => $options,
                                                      -workdir       => $workdir,
                                                      -file_prefix   => \@file_prefix
                                                     );

or

  my $obj = Bio::EnsEMBL::Pipeline::Runnable::Avid->new(
                                                        -query_seqs    => \@sequences,
                                                       );


$obj->run;
my $resultfile = $obj->results;


=head1 DESCRIPTION

Avid takes two Bio::Seq (or Bio::PrimarySeq) objects and aligns them
against each other using their FASTA-sequences and their repeatmasked
sequences. The setting of a working directory and a prefix for the written 
.fasta and .fasta.masked files is optional.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

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
  my $self = {};                #construct empty hash
  bless $self,$class;

  my (
      $query_seqs,
      $avid,
      $options,
      $workdir,
      $file_prefix ) = $self->_rearrange([qw(
                                             QUERY_SEQUENCES
                                             AVID
                                             OPTIONS
                                             WORKDIR
                                             FILE_PREFIX
                                            )
                                         ], @args);
  $self->{_verbose} = 1 ;


  ########### TESTING TYPES OF PASSED OBJECTS ###########
  $self->query_sequences($query_seqs);


  ########### LOCATION OF BINARY default av9d-2.1b0 #####
  $self->avid_binary($avid);


  ########### PASSED OPTIONS ############################
  $self->options( $options );


  ########## FILE_PREFIX ####################
  $self->file_prefix ($file_prefix);


  ########## SETTING WORKDIR ####################
  defined $workdir ? $self->workdir($workdir) : $self->workdir('/tmp/');

  return $self;
}




sub DESTROY {
  my $self = shift;
  $self->deletefiles;
}


################ RUN METHOD ###########################################

=head2 run

Usage    : $obj->run
Function : runs the avid alignment algorithm against the sequences of the given set of Bio::PrimarySeqI
objects and their according masked sequences and puts the results into the file $resultfile =$obj->results.

=cut

sub run {
  my ($self) = @_;

  # Write masked and unmasked query sequences to 4 files and store their names
  $self->write_sequences;


  # names of fasta files on which avid runs on (there must also be a .masked-file)
  my $fa_first =  @{$self->run_avid_on_files}[0];
  my $fa_secnd =  @{$self->run_avid_on_files}[1];


  my $command  = $self->avid_binary  . " " .$self->options;
  $command .=  " " . $fa_first    . " ";
  $command .=  " " . $fa_secnd         ;

  print STDERR "Avid command : $command\n"  if $self->verbose;

  open( AVID, "$command |" ) || $self->throw("Error running avid $!");
  close(AVID);


  # add the files written by avid to the @file-array for deletion
  $self->files_to_delete($fa_first,$fa_secnd);

  return 1
}


############################################################

sub files_to_delete {
  my ($self,$ff,$sf) =  @_;
  my @temp;

  my $workdir        =  $self->workdir;

  # substitution of workdir-prefix and .fasta-suffix to get filenames
  $ff             =~s/$workdir(.+)\.fasta/$1/;
  $sf             =~s/$workdir(.+)\.fasta/$1/;


  # option -opsl writes out: .info .out .psl
  # option -obin writes out: .minfo .mout

  if ($self->options =~ /obin/) {
    #    @temp     = qw ( .minfo );                            # add .mout
    $self->results( $ff . "_" . $sf . ".mout" );
  }
  if ($self->options =~ /opsl/) {
    #    @temp     = qw ( .minfo .info .out );                 # add .psl
    $self->results( $ff . "_" . $sf . ".psl" );
  }



  # adding the filenames to delete to @file-Array
  for (@temp) {
    $self->file( $ff."_". $sf . $_ );
  }
}



sub query_sequences {
  my ($self,$seq_ref) = @_;

  if ($seq_ref) {
    if (ref($seq_ref) ne "ARRAY") {
      $self->throw("[$seq_ref] is not a reference to an Array\n");
    } else {
      foreach (@{$seq_ref}) {
        $self->throw("[$_] is not a Bio::PrimarySeqI") unless $_->isa("Bio::PrimarySeqI");
      }
    }
    my @seqs                 = @{$seq_ref};
    $self->{_query_sequences} = \@seqs;
  }
  return @{$self->{_query_sequences}};
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

sub options {
  my ($self,$opt) =  @_;

  $self->{_options} = $opt if (defined $opt);
  $self->{_options} = '-obin' if ( ( !defined $opt) && ( !defined  $self->{_options}));
  return $self->{_options};
}


############################################################
# writing out 4 files for avid-input : 2 fasta-files and 2 repeatmasked fastas

sub write_sequences {
  my $self = shift;

  my $count = 0;
  foreach my $slice ( $self->query_sequences ) {

    # getting filename of (first) unmasked sequence
    my $filename = $self->get_tmp_file($self->workdir , @{$self->file_prefix}[$count] , "fasta");

    # inherited method to store filenames to unlink in @array
    $self->file($filename);

    # storing the fasta-filenames for avid-run (a little redundant...)
    $self->run_avid_on_files($filename);

    # writing unmasked fasta-sequence
    my  $seqobj = Bio::SeqIO->new(-file => ">$filename" , '-format' => 'Fasta');
    $seqobj->write_seq($slice);
    print "Writing sequence $filename\n" if ($self->verbose);

    # getting filename of (first) masked sequence 
    $filename .= ".masked";

    # inherited method to store filenames to unlink in @array
    $self->file($filename);

    # writing masked fasta-sequence
    my $rm_slice = $slice->get_repeatmasked_seq();
    $seqobj =  Bio::SeqIO->new(-file => ">$filename" , '-format' => 'Fasta');
    $seqobj -> write_seq($rm_slice);
    print "Writing sequence $filename\n" if ($self->verbose);
    $count++;
  }
}



########### STORING FILENAMES FOR AVID RUN ###########

sub run_avid_on_files {
  my ($self,$file) = @_;

  if (!defined  $self->{_run_avid_on_files} ) {
    $self->{_run_avid_on_files} = [];
  }
  if (defined $file) {
    push  @{$self->{_run_avid_on_files}} , $file ;
  }
  return $self->{_run_avid_on_files} ;
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


1;

