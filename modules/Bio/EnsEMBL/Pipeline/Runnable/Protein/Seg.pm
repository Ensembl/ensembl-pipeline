# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg

=head1 SYNOPSIS

  my $seqstream = Bio::SeqIO->new ( -file => $queryfile,
                                    -fmt => 'Fasta',
                                  );
  $seq = $seqstream->next_seq;

  my $seg = Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg->new ( -QUERY => $seq);
  $seg->workdir ($workdir);
  $seg->run;
  my @results = $seg->output;

=head1 DESCRIPTION

  Seg takes a Bio::Seq (or Bio::PrimarySeq) object
  and runs seg on it (detecting low complexity sequences). 
  The resulting output file is parsed to produce a set of features.

=head1 CONTACT
  
  Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

  The rest of the documentation details each of the object methods. 
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;

use vars qw(@ISA);
use strict;
use warnings;

use Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation;

@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation);


sub multiprotein{
  my ($self) = @_;
  return 1;
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

  $self->throw ("Error running ".$self->program." on ".$self->filename) 
    unless ((system ($self->program." ".$self->filename." -l > ".
                     $self->results)) == 0); 
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
	    #print STDERR $self->program." didn't find anything\n";
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
  while (<$filehandle>) {
    chomp;
    next if /^$/;
    if (/^\>/) {
       /^\>(\S+)?\((\d+)\-(\d+)\)\s*complexity=(\S+)/;
       my $tid = $1;
       my $start = $2;
       my $end = $3;
       my $score = $4;
       
       my $fp = $self->create_protein_feature($start, $end, $score, $tid, 
                                              0, 0, 'Seg', 
                                              $self->analysis, 0, 0);
       $self->add_to_output($fp);
     }
  }
  close $filehandle;   
}


=head2 output

 Title    : output
 Usage    : $self->output
 Function : returns an array of feature objects
 Example  :
 Returns  : an array of Bio::EnsEMBL::SeqFeature objects
 Args     :
 Throws   :

=cut

sub output {
    my ($self) = @_;
    my @list = @{$self->{'_flist'}};
    return @{$self->{'_flist'}};
}


=head2 get_low_complexity_length

 Title    : get_low_complexity_length
 Usage    : $len = $self->get_low_complexity_length;
 Function : returns *percentage* low complexity of protein
 Example  :
 Returns  : a percentage_id
 Args     :
 Throws   :
 Notes    : It only makes sense to call this method when the 
    Runnable was created with a single Bio::Seq

=cut



sub get_low_complexity_length {
  my ($self) = @_;

  if ($self->query->length > 0) {    
    my $lc_length = 0;

    foreach my $feat ($self->output) {
      $lc_length += abs($feat->end - $feat->start) + 1;
    }
    
    my $low_complexity = ($lc_length)/($self->query->length);
    
    $low_complexity *= 100;
    
    return $low_complexity;
  }
  else {
    return 0;
  }
}

		
1;
