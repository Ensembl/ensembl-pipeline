# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

 Bio::EnsEMBL::Pipeline::Runnable::Protein::Tmhmm

=head1 SYNOPSIS

 my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                   -fmt => 'Fasta',
                                 );
 $seq = $seqstream->next_seq;

 my $tmhmm = Bio::EnsEMBL::Pipeline::Runnable::Protein::Tmhmm->new ( -CLONE => $seq);
 $tmhmm->workdir ($workdir);
 $tmhmm->run;
 my @results = $tmhmm->output;

=head1 DESCRIPTION

 Tmhmm takes a Bio::Seq (or Bio::PrimarySeq) object
 and runs tmhmm on it (detecting transmembrane segments). 
 The resulting output file is parsed to produce a set of features.

=head1 CONTACT

 Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods. 
 Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Tmhmm;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation;


@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation);





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
    # run program
    print STDERR "running ".$self->program." ".$self->filename."\n";
    $self->throw ("Error running ".$self->program." on ".$self->filename) 
        unless ((system ("perl ".$self->program." ".$self->filename.
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
  my $id;
  while (<$filehandle>) {
    #print STDERR;
    chomp;
    next if /^$/;
    if (/^\>(\S+)/) {
      $id = $1;
      print STDERR "have id ".$id."\n";
    }
    elsif (/^%pred/) {
      my ($junk, $values) = split /:/;
      my @tm = split (/,/, $values);
      foreach (@tm) {
        /(\w+)\s+(\d+)\s+(\d+)/;
        my $orien = $1;
        my $start = $2;
        my $end = $3;
        $orien = uc ($orien);
        if ($orien eq "M") {
          my $fp = $self->create_protein_feature($start, $end, 0, $id, 0, 
                                                 0, 'Tmhmm', 
                                                 $self->analysis, 0, 0);
          $self->add_to_output($fp);
        }
	    }
    }
  }
  close $filehandle;   
}

sub multiprotein{
  my ($self) = @_;
  return 1;
}


1;
