# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# based on Michelle Clamp's Blast.pm
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => 'Pfam');

  $hmm->workdir ($workdir);
  $hmm->run;
  my @results = $hmm->output;

=head1 DESCRIPTION

  Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs hmmpfam.
  The resulting output file is parsed to produce a set of Bio::EnsEMBL::FeaturePairs.

=head1 CONTACT

   Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;

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
    print STDERR "running ".$self->program." against ".$self->database."\n";

    # some of these options require HMMER 2.2g (August 2001)
    
    my @dbfiles = split(/;/,$self->database);

    print STDERR "FILENAME: ".$self->filename."\n";

    if ($dbfiles[0] =~ /ls/) {
      my $cmd = $self->program .' --acc --cut_ga --cpu 1 '.
        $self->parameters .' '.$dbfiles[0].' '.$self->filename.' > '.
          $self->lsresults;
      print STDERR "$cmd\n";   
      $self->throw ("Error running ".$self->program." on ".
                    $self->filename." against ".$dbfiles[0])unless 
                      ((system ($cmd)) == 0);
    }else {
      die || "ls pfam file has not been provided";
    }

    if ($dbfiles[1] =~ /fs/) { 
       
      my $cmd = $self->program .' --acc --cut_ga --cpu 1 '.
	        $self->parameters.' '.$dbfiles[1]      .' '.$self->filename.
            ' > '.$self->fsresults;
      print STDERR "$cmd\n";   
      $self->throw ("Error running ".$self->program." on ".
                    $self->filename." against ".$dbfiles[1]) unless 
                      ((system ($cmd)) == 0);
    }else {
      die || "fs pfam file has not been provided";
    }
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
  my $fshandle;
  my $id;
  my $resfile = $self->lsresults;
  my $fsfile  = $self->fsresults;
  
		

  if (-e $resfile) {
    # it's a filename
    if ((-z $self->lsresults) && (-z $self->fsresults)) {  
      print STDERR $self->program." didn't find anything\n";
      return;
    }       
    else {
	    open (OUT, "<$resfile") or $self->throw ("Error opening $resfile");
      $filehandle = \*OUT;
      
      
	    open (OUT1, "<$fsfile")  or $self->throw ("Error opening $fsfile");
	    $fshandle = \*OUT1;
    }
  }else {
    # it'a a filehandle
    $filehandle = $resfile;
    $fshandle = $fsfile;
  }
  
  
  
  #First parse what comes from the ls mode matches. 
  #Every match in that case is taken
  while (<$filehandle>) {
    chomp;
    last if /^Alignments of top-scoring domains/;
    next if (/^Model/ || /^\-/ || /^$/);
    if (/^Query sequence:\s+(\S+)/) {
      $id = $1;
    }
    
    if (my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
      
      my $pvalue = sprintf ("%.3e", $evalue);
   
      my $fp = $self->create_protein_feature($start, $end, $score, 
                                             $id, $hstart, $hend, $hid,
                                             $self->analysis, $pvalue,
                                             0);
      $self->add_to_output($fp);
    }
  }
  close FILEHANDLE; 
  
  #Then read all of the fs mode matches. If a match does not overlap with 
  #any ls match thus its taken
  while (<$fshandle>) {
    my ($hid, $start, $end, $hstart, $hend, $score, $evalue);
    my $overlap = undef;
      
    chomp;
      
    last if /^Alignments of top-scoring domains/;
    next if (/^Model/ || /^\-/ || /^$/);
    if (/^Query sequence:\s+(\S+)/) {
      $id = $1;
    }
    if (($hid, $start, $end, $hstart, $hend, $score, $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
      
      
      foreach my $featpair(@{$self->{'_flist'}}) {
        my $lsstart = $featpair->feature1->start;
        my $lsend = $featpair->feature1->end;
        if ((($start >= $lsstart) && ($start <= $lsend)) || (($end >= $lsstart) && ($end <= $lsend))) {
          $overlap = 1;
        }
      }
      
      if (!$overlap) {
        my $pvalue = sprintf ("%.3e", $evalue);
        my $fp = $self->create_protein_feature($start, $end, $score, 
                                               $id, $hstart, $hend, $hid,
                                               $self->analysis, $pvalue,
                                               0);
        $self->add_to_output($fp);
	    }
    }
  }
  close (FS);
}





sub lsresults {
    my ($self) = @_;
    if(!$self->{'_lsresults'}){
      $self->{_lsresults} = $self->results.'.lsout';
    }
    return $self->{_lsresults};
}


sub fsresults {
    my ($self, $results) = @_;
    if(!$self->{'_fsresults'}){
      $self->{_fsresults} = $self->results.'.fsout';
    }
    return $self->{_fsresults};
}

1;
