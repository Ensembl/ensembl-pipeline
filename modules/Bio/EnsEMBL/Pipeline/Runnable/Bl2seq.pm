#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Bl2seq

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new(-seq1 => $seq1,
							    -seq2 => $seq2,
							    -program => $program,
							    -alntype => $alntype
							    -min_score => $min_score,
							    -workdir => $workingdir,
							    -results => $results);
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new().

    $seq1 and $seq2 must be Bio::PrimarySeq object.

    $program (optional) must be a sting which locates bl2seq executable.

    $alntype (optional) must be a string which defines the alignment type, blastp, 
              blastn, blastx, tblastx or tblastn (default = blastn).

    $min_score (optional) must be a interger, value for minimum score (bits) for a HSP.

    $workingdir (optional) must be a string which locates the working directory.

    $results (optional) must be a string which locates the results file.

=head1 DESCRIPTION

This runs bl2seq (executable from the ncbi) and provides feature pair output. 
ONLY parsing of 'blastn' output is implemented.

=head2 Methods:

 new,
 seq1,
 seq2,
 program,
 alntype,
 min_score,
 workdir,
 run,
 output

=head1 CONTACT

ensembl-dev@ebi.ec.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::Bl2seq;

use strict;

use Bio::EnsEMBL::Pipeline::Tools::Bl2seq;
use Bio::EnsEMBL::Pipeline::RunnableI;

use vars qw(@ISA);

# Object preamble - inherits from Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

  $self->{'_seq1'} = undef; # location of Bio::Seq object "-i" option 
  $self->{'_seq2'} = undef; # location of Bio::Seq object "-j" option
  $self->{'_program'} = "bl2seq"; # location of bl2seq executable
  $self->{'_alntype'} = "blastn"; # type of alignment "-p" option
  $self->{'_min_score'} = 40; # value for minimum score 
  $self->{'_fplist'} = []; # an array of feature pairs (the output)
  $self->{'_workdir'} = "/tmp"; # location of working directory
  $self->{'_results'} = $self->{'_workdir'}."/results.".$$; # location of result file
  $self->{'_options'} = "";

  my ($seq1, $seq2, $program, $alntype, $options, $min_score, $workdir, $results) =    
    $self->_rearrange([qw(SEQ1 SEQ2 PROGRAM ALNTYPE OPTIONS MIN_SCORE WORKDIR RESULTS)], @args);

  $self->throw("Must pass in both seq1 and seq1 args") if (! defined $seq1 || ! defined $seq2);
  $self->seq1($seq1);
  $self->seq2($seq2);

  $self->program($self->find_executable($program)) if ($program);
  $self->alntype($alntype) if ($alntype);
  $self->options($options) if ($options);
  if ($workdir) {
    $self->workdir($workdir);
    $self->results($self->workdir."/results.".$$);
  }
  $self->results($results) if ($results);
  $self->min_score($min_score) if (defined $min_score);
  
  return $self;
}

=head2 seq1

 Title   : seq1
 Usage   : $obj->seq1($newval)
 Function: 
 Example : 
 Returns : value of seq1
 Args    : newvalue (optional)


=cut

sub seq1 {
  my ($obj,$value) = @_;
  if( defined $value) {
    if (! ref $value || ! $value->isa('Bio::PrimarySeqI')) {
      $obj->throw("$value is not a PrimarySeqI object. Cannot throw");
    }
    my $selfseq = Bio::PrimarySeq->new(-display_id => $value->id,
				       -seq => $value->seq);
    if ($selfseq->length == 0) {
      $obj->throw("attempting to bl2seq seemingly 0 length sequence!");
    }
    $obj->{'_seq1'} = $selfseq;
  }
  return $obj->{'_seq1'};
}

=head2 seq2

 Title   : seq2
 Usage   : $obj->seq2($newval)
 Function: 
 Example : 
 Returns : value of seq2
 Args    : newvalue (optional)


=cut

sub seq2 {
   my ($obj,$value) = @_;
   if( defined $value) {
      if( !ref $value || !$value->isa('Bio::PrimarySeqI') ) {
	  $obj->throw("$value is not a PrimarySeqI object. Cannot throw");
      }
      my $selfseq = Bio::PrimarySeq->new( -display_id => $value->id , -seq => $value->seq);
      if( $selfseq->length == 0 ) {
	  $obj->throw("attempting to bl2seq seemingly 0 length sequence!");
      }
      $obj->{'_seq2'} = $selfseq;

    }
    return $obj->{'_seq2'};

}

=head2 program

    Title   :   program
    Usage   :   $obj->program('/usr/local/ensembl/bin/bl2seq');
    Function:   Get/set method for the location of the bl2seq executable
    Returns :   string
    Args    :   string

=cut

sub program {
  my ($self, $location) = @_;

  if ($location) {
    $self->throw("executable not found at $location: $!\n") unless (-e $location && -x $location);
    $self->{'_program'} = $location ;
  }
  return $self->{'_program'};
}

=head2 alntype

 Title   : alntype
 Usage   : $obj->alntype($newval)
 Function: Get/Set the alntype value defining the alignment type
           (blastp,blastn,blastx,tblastx,tblastn)
 Example : $obj->alntype('blastn')
 Returns : string
 Args    : string (optional)


=cut

sub alntype {
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'_alntype'} = $value;
  }
  return $self->{'_alntype'};
}

=head2 options

 Title   : options
 Usage   : $obj->options($newval)
 Function: Get/Set the options values
           ("-G 1 -E 2")
 Example : $obj->options('-G 1 -E 2')
 Returns : string
 Args    : string (optional)


=cut

sub options {
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'_options'} = $value;
  }
  return $self->{'_options'};
}

=head2 min_score

 Title   : min_score
 Usage   : $obj->min_score($newval)
 Function: Get/Set the min_score value
 Example : 
 Returns : value of min_score used to filter bl2seq results
 Args    : interger (optional)


=cut

sub min_score {
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'_min_score'} = $value;
  }
  return $self->{'_min_score'};
}

=head2 run

 Title   : run
 Usage   : $self->run();
 Function: run bl2seq program (NCBI) on two previously defined sequences
 Example :
 Returns : 1
 Args    : None


=cut

sub run {
  my ($self) = @_;
  
  # dump sequences to work directory
  my $query = $self->workdir."/query.".$$;
  my $sbjct = $self->workdir."/sbjct.".$$;

  $self->filename($query);
  $self->writefile($self->seq1);
  $self->filename($sbjct);
  $self->writefile($self->seq2);

  # run blast and parse results
  $self->run_analysis($query,$sbjct);
  $self->parse_results();   

  # delete sequence files
  unlink($query);
  unlink($sbjct);
  unlink($self->results);

  return 1;
}

sub run_analysis {
  my ($self,$query,$sbjct) = @_;
  print STDERR ("Running bl2seq...\n" . $self->program .
		                     " -i $query" .
		                     " -j $sbjct " .
		                     $self->options .
		                     " -p " . $self->alntype .
		                     " > ". 
		                     $self->results. "\n");

  $self->throw("Failed during bl2seq run, $!\n") unless (system ($self->program .
								 " -i $query" .
								 " -j $sbjct " .
								 $self->options .
								 " -p " . $self->alntype .
								 " > ". 
								 $self->results) == 0);
}

sub parse_results { 
  my ($self) = @_;

  open BL2SEQ, $self->results || 
    $self->throw("Coudn't open file ".$self->results.", $!\n");
  my $filehandle = \*BL2SEQ;

  my $bl2seq_parsing = Bio::EnsEMBL::Pipeline::Tools::Bl2seq->new('-fh' => $filehandle,
					   '-alntype' => $self->alntype,
					   '-min_score' => $self->min_score,
					   '-qname' => $self->seq1->id,
					   '-sname' => $self->seq2->id);

  while (my $DnaDnaAlignFeature = $bl2seq_parsing->nextHSP) {
    $self->_add_fp($DnaDnaAlignFeature);
  }
}

sub _add_fp {
  my ($self,@args) = @_;
  if (@args) {
    push(@{$self->{'_fplist'}},@args);
  } else {
    warn "WARN: Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->_add_fp should have an argument\n";
  }
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   None

=cut

sub output {
  my ($self) = @_;
  return @{$self->{'_fplist'}};
}

=head2 workdir

 Title   : workdir
 Usage   : $obj->workdir($newval)
 Function: 
 Example : 
 Returns : value of workdir
 Args    : newvalue (optional)


=cut

sub workdir{
   my ($self,$value) = @_;
   if( defined $value) {
       $self->{'_workdir'} = $value;
   }
   return $self->{'_workdir'};
}

1;
