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
@sequences = a list of Bio::Seq objects,
$avid      = a location for the binary,
$options   = a string with options (default is to produce binary output)

  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Avid->new(
								 -query_seqs    => \@sequences,
                                                                 -avid          => $avid,
								 -options       => $options,
								);

 $runnable->run;
 my @results = $runnable->output;


=head1 DESCRIPTION

Avid takes two Bio::Seq (or Bio::PrimarySeq) objects and aligns them
against each other using their FASTA-sequences and their repeatmasked
sequences.  The resulting output file is parsed to produce a set of features.


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


#################### NEW ####################

sub new {

  my ($class,@args) = @_;
  my $self={};  #construct empty hash
  bless $self,$class;

  my (
      $query_seqs,
      $avid,
      $options) = $self->_rearrange([qw(
					QUERY_SEQUENCES
					AVID
					OPTIONS
				       )
				    ], @args);
  $self->{_output} = [];
  $self->{_verbose} = 1 ;




  ########### TESTING TYPES OF PASSED OBJECTS ###########

  if(ref($query_seqs) ne "ARRAY") {
    $self->throw("[$query_seqs] is not a reference to an Array\n"); 
  }else{
    foreach(@{$query_seqs}) {
	    $self->throw("[$_] is not a Bio::PrimarySeqI") unless $_->isa("Bio::PrimarySeqI");
    }
  }
  $self->query_sequences($query_seqs);



  ########### LOCATION OF USED BINARY ###################
  # if not given, we use the  default:avid-2.1b0

  if ($avid) {
      $self->avid_binary($avid);
  } else {
      $self->avid_binary('/nfs/acari/jhv/bin/avid');
  }



  ########### PASSED OPTIONS ############################

  print "options $options\n";
  if ($options) {
    $self->options( $options );
  }else{
    $self->options("-obin");
  }



  $self->filenames($$);
  $self->avid_files($$);


  $self->printvars;
  return $self;
}




sub DESTROY {
  my $self = shift;
  unlink $self->filenames;
}



############################################################
#
# Analysis methods
#
############################################################

=head2 run

Usage   :   $obj->run($workdir, $args)
Function:   Runs avid script and puts the results into the file $self->results
            It calls $self->parse_results, and results are stored in $self->output
=cut

sub run {
  my ($self) = @_;

  # Set name of results file containing PID
  $self->results(     $self->workdir . "/results.$$"     );

  # Write query sequences into files
   $self->write_sequences;

  my ($file1, $file2) = ($self->filenames)[0,2];

  my $command   = $self->avid_binary  . " " .$self->options;
     $command .=  " " . $file1 . " " . $file2;


  $self->avid_files(333);

  print STDERR "Avid command : $command\n"  if $self->_verbose;


  open( EXO, "$command |" ) || $self->throw("Error running avid $!");
#  $self->parse_results(\*EXO);
  close(EXO);
  return 1
}



############################################################
# get/set methods

sub query_sequences {
  my ($self,$seq_ref)= @_;
  if ($seq_ref) {
    my @seqs           = @{$seq_ref};
    $self->{_query_sequences} = \@seqs;
  }
  return @{$self->{_query_sequences}};
}

############################################################



sub avid_binary {
  my ($self, $location) = @_;
  if ($location) {                #if two values are submitted check if the submit. binary exists
    $self->throw("Slam not found at $location: $!\n") unless (-e $location);
    $self->{_slam} = $location ;
  }
  return $self->{_slam};
}

############################################################

sub options {
  my ($self,$opt) = @_; 
    $self->{_options}= $opt if $opt;
  return $self->{_options};
}

############################################################

sub output {
  my ($self, @output) = @_;

  if (@output) {
    unless( $self->{_output} ){
      $self->{_output} = [];
    }
    push( @{$self->{_output}}, @output );
  }
  return @{$self->{_output}};
}

############################################################
# writing out 4 files for avid-input : 2 fasta-files and 2 repeatmasked fastas 


sub write_sequences {
  my $self   = shift;

   my @filenames   = $self->filenames;

  foreach my $slice ( $self->query_sequences ){
    my $file       = shift  @filenames;
    my  $seqobj    =  Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
        $seqobj    -> write_seq($slice);
        print "Writing out sequence $file\n" if ($self->_verbose);


       $file       =  shift  @filenames;
    my $rm_slice   = $slice->get_repeatmasked_seq();
       $seqobj     = Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
       $seqobj     -> write_seq($rm_slice);
        print "Writing out sequence $file\n" if ($self->_verbose);
  }
}



########### BUILDING FILENAMES WITH PATH ###########

sub filenames {
  my ($self,$pid) = @_;

  my @filenames;

  if ( $pid ) {
    for(my $i=1;$i<=2;$i++) {
      push @filenames, $self->workdir."/seq$i.$pid.fasta"; 
      push @filenames, $self->workdir."/seq$i.$pid.fasta.masked";
    }
    $self->{_filenames}=\@filenames;
  }
  return @{$self->{_filenames}} ;
}


sub avid_files {
  my ($self,$pid) = @_;
 ;
}

sub _verbose {
  my $self = shift;
    $self->{_verbose} = shift if @_ ;
    return $self->{_verbose}
}


sub seq_names {
  my ($self, $fn ) = @_ ;
      my @names = @{$fn};
  return; 
}

sub printvars {
  my $self = shift;
  my %obj = %{$self};
  print map { "$_ => $obj{$_}\n" }  keys %obj ;
}




1;

