# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

 Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp

=head1 SYNOPSIS

 my $seqstream = Bio::SeqIO->new ( -file => $queryfile,
                                   -fmt => 'Fasta',
                                 );
 $seq = $seqstream->next_seq;

 my $signalp = Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp->new ( -QUERY => $seq);
 $signalp->workdir ($workdir);
 $signalp->run;
 my @results = $signalp->output;

=head1 DESCRIPTION

 Signalp takes a Bio::Seq (or Bio::PrimarySeq) object
 and runs signalp on it (detecting signal peptides). 
 The resulting output file is parsed to produce a set of features.

=head1 CONTACT

 Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX
 
 The rest of the documentation details each of the object methods. 
 Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation);




sub multiprotein{
  my ($self) = @_;
  return 0;
}

sub write_protein_file{
  my ($self, $filename, $seq) = @_;

  open (OUT,">".$filename) or throw("couldn't open $filename");
  my $sub_seq = substr ($seq->seq, 0, 50);
	$seq->seq($sub_seq);
  print OUT ">".$seq->display_id."\n".$seq->seq."\n";
  close(OUT) or throw("couldn't close $filename");

  return $filename;
}

=head2 run_program

 Title    : run_program
 Usage    : $self->program
 Function : makes the system call to program
 Example  :
 Returns  : 
 Args     :
 Throws   :

=cut

sub run_analysis {
    my ($self) = @_;
    
    print STDERR "RUNNING: ".$self->program." -t euk ".$self->filename.
      " > ".$self->results."\n";
    
    # run program
    
    $self->throw ("Error running ".$self->program." on ".$self->filename) 
      unless ((system ($self->program." -t euk ".$self->filename. 
                       " > ".$self->results)) == 0); 
  }


=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses program output to give a set of features
 Example  :
 Returns  : 
 Args     : filename (optional, can be filename, filehandle or pipe, not implemented)
 Throws   :

=cut

sub parse_results {
  my ($self) = @_;
  my $filehandle;
  my $resfile = $self->results;
  
  if (-e $resfile) {
    # it's a filename
    if (-z $self->results) {  
	    print STDERR $self->program." didn't find anything\n";
	    return;
    }else {
      open (OUT, "<$resfile") or $self->throw ("Error opening $resfile");
      $filehandle = \*OUT;
    }
  }else {
    # it'a a filehandle
    $filehandle = $resfile;
  }
  
  # parse
  my ($id, $fact1, $fact2, $end);
  while (<$filehandle>) {
    chomp;
    if (/^\>(\S+)/) {
      $id = $1;
    }
    elsif (/max\.\s+Y\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
      $fact1 = $2;
    }
    elsif (/mean\s+S\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
      $fact2 = $2;
      if ($fact1 eq "YES" && $fact2 eq "YES") {
        my $line = <$filehandle>;
        if ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
          $end = $1;
        }
        else {
          $self->throw ("parsing problem in ".$self->program);
        }
        my $fp = $self->create_protein_feature(1, $end, 0, $id, 0, 0,
                                               'Sigp', $self->analysis,
                                               0, 0);
        $self->add_to_output($fp);
	    }
    }
  }
}



1;
