#
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

Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs->new('-genomic'    => $genseq,
								   '-ests'       => $ests,
								   '-exonerate'  => $exonerate,
								   '-exonerate_args' => $args,
								  )

    $obj->run;

    my @features = $obj->output;


=head1 DESCRIPTION
Runs exonerate between a contig sequence and an estfile. 
Exonerate output is printed to STDOUT by Exonerate - this is controlled by a flag to the Exonerate constructor.

Runnable:ExonerateESTs sits between RunnableDB::ExonerateESTs and Runnable::Exonerate

we currently call Runnable::Exonerate::run with a flag for ungapped; this will probably change as Exonerate matures.


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::PrimarySeqI;
use Bio::SeqIO;
use Bio::DB::RandomAccessI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-GENOMIC    => $genomicseq,
			   -ESTS       => $ests,
			   -EXONERATE  => $exonerate,
			   -EXONERATE_ARGS => $exargs);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs object
    Returns :   A Bio::EnsEMBL::Pipeline::Runnable::ExonerateESTs object
    Args    :   -genomic:        Bio::PrimarySeqI object (genomic sequence)
                -ests:           Either path to file containing est seqs or reference to array of Bio::Seq
                -exonerate:  path to exonerate executable (optional)
                -exonerate_args: arguments to be passed to exonerate. memory, word size etc
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my( $genomic, $ests, $exonerate, $exargs ) = $self->_rearrange([qw(GENOMIC
								       ESTS
								       EXONERATE
								       EXONERATE_ARGS)],
					      @args);
							   
    $self->throw("No genomic sequence input")           
      unless defined($genomic);
    $self->genomic_sequence($genomic) if defined($genomic);

    $self->throw("No ests specified") 
      unless defined($ests);
    $self->ests($ests) if defined($ests);

    $self->exonerate($exonerate) if defined $exonerate;

    $self->exonerate_args($exargs) if defined($exargs);

    return $self; 
}

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object, or a filename

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
      $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 ests

    Title   :   ests
    Usage   :   $self->ests($ests)
    Function:   Get/set method for ests
    Returns :   
    Args    :   name of a file containing est seq(s), OR reference to an array of Bio::Seq

=cut

sub ests {
    my( $self, $ests ) = @_;   
    if ($ests) { 
      if (ref($ests) eq 'ARRAY') {

	# I'm not at all sure this is right
	my $time = time; chomp($time);
	my $estfile = "/tmp/estfile_.$$.$time.fn";
	$self->estfilename($estfile);

	foreach my $est(@$ests) {
	  $est->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
	}

	$self->{'_ests_sequences'} = $ests;
      }
      else {
	# it's a filename - check the file exists
	$self->throw("[$ests] : file does not exist\n") unless -e $ests;
	$self->estfilename($ests);
	$self->{'_est_sequences'} = $ests;
    }
  }
  
  #NB ref to an array of Bio::Seq
  return $self->{'_est_sequences'};

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


=head2 estfilename

    Title   :   estfilename
    Usage   :   $self->estfilename($filename)
    Function:   Get/set method for estfilename
    Returns :   
    Args    :   

=cut

sub estfilename {
  my ($self, $filename) = @_;
  $self->{'_estfilename'} = $filename if ($filename);
  return $self->{'_estfilename'};
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
  Function: Runs exonerate vs input ests
  Returns : array of Bio::EnsEMBL::FeaturePair
  Args    : 

=cut

sub run_exonerate {
  my ($self) = @_;
  my @res;

  my $estseq  = $self->ests;
  my $estfile = $self->estfilename;

  # do we need to write out the est sequences?
  if(ref($estseq) eq 'ARRAY'){
    eval{
      if (-e $estfile) { $self->throw("alreayd using $estfile\n"); }
      my $estOutput = Bio::SeqIO->new(-file => ">$estfile" , '-format' => 'Fasta')
	or $self->throw("Can't create new Bio::SeqIO from $estfile '$' : $!");
      
      foreach my $eseq(@$estseq) {
	$estOutput->write_seq($eseq);
      }
    };
    
    if($@){
      $self->warn("couldn't run exonerate - problem writing estfile\n");
      return;
    }

  }

  my $exr = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
							     '-exonerate' => $self->exonerate,
							     '-genomic'   => $self->genomic_sequence,
							     '-est'       => $self->estfilename,
							     '-args'      => $self->exonerate_args,
							     '-print'     => 1,
							   );
  
  my $ungapped = 1;
  $exr->run($ungapped);

  my $res = $exr->output; # ref to an array

  # clean up temp files
  if(ref($estseq) eq 'ARRAY'){
    unlink $estfile;
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


