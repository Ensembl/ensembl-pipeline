#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneComparison::cDNA_Comparison

=head1 SYNOPSIS

  my $obj = Bio::EnsEMBL::Pipeline::GeneComparison::cDNA_Comparison->new( 					
									 '-transcripts'    => $transcsripts,
									 '-cdnas'          => $cdnas,
									 '-exonerate'      => $exonerate,
									 '-exonerate_args' => $args,
									);
 
$obj->run;

my @features = $obj->output;

=head1 DESCRIPTION

Runs exonerate between two set of cDNAs. Tipically one set corresponds to the cDNAs of some set of
predicted transcripts and the other set corresponds to good quality full-length cdnas to compare to.
The reasons to compare them can be:
We want to know how many extra genes are in the full-length cdna set.
We want to run some quality checks on the produced transcripts.

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneComparison::cDNA_Comparison;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::PrimarySeqI;
use Bio::SeqIO;
#use Bio::DB::RandomAccessI;
use Bio::EnsEMBL::Pipeline::ESTConf qw(
				       EST_EXONERATE_ARGS
				      );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

# arguments
#  -transcripts: all the transcripts (in fasta file?) (Bio::PrimarySeqI objects?)
#  -cdnas:      a bunch of cdnas    (in fasta file?) (Either path to file containing est seqs or reference to array of Bio::Seq?)
#  -exonerate:  path to exonerate executable (optional)
#  -exonerate_args: arguments to be passed to exonerate. memory, word size etc

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my( $transcripts, $cdnas, $exonerate, $exargs ) = $self->_rearrange([qw(TRANSCRIPTS
									  CDNAS
									  EXONERATE
									  EXONERATE_ARGS)],
								      @args);
  
  print STDERR "transcripts: $transcripts cdnas: $cdnas exonerate: $exonerate exargs: $exargs\n"; 
  print STDERR "input_id: ".$self->input_id."\n";

  $transcripts = $self->input_id unless defined($transcripts);
  $self->throw("No transcript sequences in the input")           
    unless defined($transcripts);
  $self->transcripts($transcripts) if defined($transcripts);
  
  $self->throw("No cdnas specified") 
    unless defined($cdnas);
  $self->cdnas($cdnas) if defined($cdnas);
  
  $self->exonerate($exonerate) if defined $exonerate;
  
  $self->exonerate_args($exargs) if defined($exargs);
  
  unless ( $self->exonerate_args ){
    $self->exonerate_args($EST_EXONERATE_ARGS);
  }
  
  return $self; 
}

=head2 transcripts

Function:   Get/set method for transcripts

=cut
  
sub transcripts {
  my( $self, $value ) = @_;    
  if ($value) {
    $self->{'_transcripts'} = $value;
  }
  return $self->{'_transcripts'};
}

=head2 cdnas

    Function:   Get/set method for cdnas
    Args    :   name of a file containing cdna seq(s), OR reference to an array of Bio::Seq???

=cut

sub cdnas {
    my( $self, $cdnas ) = @_;   
    if ($cdnas) { 
      if (ref($cdnas) eq 'ARRAY') {

	# I'm not at all sure this is right
	my $time = time; chomp($time);
	my $cdnafile = "/tmp/cdnafile_.$$.$time.fn";
	$self->cdnafilename($cdnafile);

	foreach my $cdna(@$cdnas) {
	  $cdna->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
	}

	$self->{'_cdna_sequences'} = $cdnas;
      }
      else {
	# it's a filename - check the file exists
	$self->throw("[$cdnas] : file does not exist\n") unless -e $cdnas;
	$self->cdnafilename($cdnas);
	$self->{'_cdna_sequences'} = $cdnas;
    }
  }
  
  #NB ref to an array of Bio::Seq
  return $self->{'_cdna_sequences'};

  }

=head2 exonerate

 Title   : exonerate
 Usage   : $obj->exonerate($exonerate)
 Function: get/set for exonerate
 Returns : path to exonerate
 Args    : exonerate (optional)


=cut

sub exonerate {
   my ($self, $exonerate) = @_;

   if (!defined $self->{'_exonerate'}){
      $self->{'_exonerate'} = "";
   }

   if( defined $exonerate ) {
      $self->{'_exonerate'} = $exonerate;
    }
    return $self->{'_exonerate'};

}


=head2 exonerate_args

  Title   :   exonerate_args
  Usage   :   $self->exonerate_args($args)
  Function:   Get/set method for arguments to exonerate
  Returns :   string
  Args    :   string (optional)

=cut

sub exonerate_args {
    my( $self, $exargs ) = @_;    
    if(!defined $self->{'_exonerate_args'}){
      $self->{'_exonerate_args'} = "";
    }
    if ($exargs) {
      $self->{'_exonerate_args'} = $exargs;
    }
    return $self->{'_exonerate_args'};
  }


sub cdnafilename {
  my ($self, $filename) = @_;
  $self->{'_cdnafilename'} = $filename if ($filename);
  return $self->{'_cdnafilename'};
}


sub fetch_input {
my ($self) = @_;
return;
}




=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs exonerate vs input ests, prints out the results, filters them and returns the filtered results
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;
    # filter ESTs using exonerate
    my $exonerate_res = $self->run_exonerate(); # ref to an array
}

=head2 run_exonerate

  Title   : run_exonerate
  Usage   : $self->run_exonerate()
  Function: Runs exonerate vs input cdnas
  Returns : array of Bio::EnsEMBL::FeaturePair
  Args    : 

=cut

sub run_exonerate {
  my ($self) = @_;
  my @res;

  my $cdnaseqs  = $self->cdnas;
  my $cdnafile = $self->cdnafilename;

  # do we need to write out the est sequences?
  if(ref($cdnaseqs) eq 'ARRAY'){
    eval{
      if (-e $cdnafile) { 
	$self->throw("already using $cdnafile\n"); 
      }
      my $cdnaOutput = Bio::SeqIO->new(-file => ">$cdnafile" , '-format' => 'Fasta')
	or $self->throw("Can't create new Bio::SeqIO from $cdnafile '$' : $!");
      
      foreach my $cdna_seq(@$cdnaseqs) {
	$cdnaOutput->write_seq($cdna_seq);
      }
    };
    
    if($@){
      $self->warn("couldn't run exonerate - problem writing cdnafile\n");
      return;
    }
    
  }
  
  my $exr = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
							     '-exonerate' => $self->exonerate,
							     '-genomic'   => $self->transcripts,
							     '-est'       => $self->cdnafilename,
							     '-args'      => $self->exonerate_args,
							     '-print'     => 1,
							   );
  
  my $ungapped = 1;
  $exr->run($ungapped);

  my $res = $exr->output; # ref to an array

  # clean up temp files
  if(ref($cdnaseqs) eq 'ARRAY'){
    unlink $cdnafile;
  }

  return $res;
  
}
    
=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
  my ($self,$feat) = @_;
  
  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }
  
  if(defined $feat){
    push(@{$self->{'_output'}},@{$feat});
  }
  
  return $self->{'_output'}; #ref to an array
}

1;


