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
  # run program
  #print STDERR "Running ".$self->program." ".$self->filename." -l > ".
  #  $self->results."\n";
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


=head2 create_feature

 Title    : create_feature
 Usage    : $self->create_feature ($feature)
 Function : creates a Bio::EnsEMBL::SeqFeature object from %feature,
            and pushes it onto @{$self->{'_flist'}}
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub create_feature {
    my ($self, $feat) = @_;

    my $analysis = $self->analysis;


    # create feature object
    my $feat1 = Bio::EnsEMBL::SeqFeature->new ( -seqname     => $feat->{name},
						-start       => $feat->{start},
						-end         => $feat->{end},
						-score       => 0,
						-analysis    => $analysis,
						-percent_id => 0,
						-p_value => 0,
						); 



    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => 0,
					      -end => 0,
					      -analysis => $analysis,
					      -seqname => 'Seg');
    
    
    my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
						-feature2 => $feat2);

    if ($feature) {
	push (@{$self->{'_flist'}}, $feature);
    }
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

sub get_low_complexity_length {
	my ($self) = @_;

	my $lc_length = 0;
  my ($p, $f, $l) = caller;
	foreach my $feat ($self->output) {
		$lc_length += abs($feat->end - $feat->start) + 1;
	}
    
	my $low_complexity = ($lc_length)/($self->query->length);
  #print STDERR "Have lc_length ".$lc_length." and query length ".
  #  $self->query->length." $f:$l\n";
  #print STDERR "Have low complexity ".$low_complexity."\n";
	$low_complexity *= 100;
  #print STDERR "Have low complexity*100 ".$low_complexity."\n";
	return $low_complexity;
}

		
1;
