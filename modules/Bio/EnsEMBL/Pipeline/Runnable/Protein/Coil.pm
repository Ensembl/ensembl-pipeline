# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

 Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil

=head1 SYNOPSIS

 my $seqstream = Bio::SeqIO->new ( -file => $queryfile,
                                   -fmt => 'Fasta',
                                 );
 $seq = $seqstream->next_seq;

 my $ncoils = Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil->new ( -QUERY => $seq);
 $ncoils->workdir ($workdir);
 $ncoils->run;
 my @results = $ncoils->output;

=head1 DESCRIPTION

 Coil takes a Bio::Seq (or Bio::PrimarySeq) object
 and runs ncoils on it (detecting coiled coils). 
 The resulting output file is parsed to produce a set of features.

=head1 CONTACT

 Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods. 
 Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil;

use vars qw(@ISA);
use strict;

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

    my $coilsdir='/usr/local/ensembl/data/coils';
    $ENV{'COILSDIR'}=$coilsdir;
    my $command = $self->program." -f < ".$self->filename." > ".$self->results;
    print STDERR "command ".$command."\n";
    # run program
    print STDERR "running ".$self->program."\n";
    $self->throw ("Error running ".$self->program." on ".$self->filename) 
        unless (system($command) == 0);
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
    }       
    else {
      open (OUT, "<$resfile") or $self->throw ("Error opening $resfile");
      $filehandle = \*OUT;
    }
  }
  else {
    # it'a a filehandle
    $filehandle = $resfile;
  }
  
  # parse
  my %result_hash = _read_fasta ($filehandle);
  close $filehandle;   
  
  foreach my $id (keys %result_hash) {
		
    my $pep = reverse ($result_hash{$id});
    my $count = my $switch = 0;
    my ($start, $end);
    while (my $aa = chop $pep) {
      $count++;
      if (!$switch && $aa eq "x") {
        $start = $count;
        $switch = 1;
      }elsif ($switch && $aa ne "x") {
        $end = $count-1;
        my $fp = $self->create_protein_feature($start, $end, 0, $id, 0, 0,
                                               'ncoils', $self->analysis,
                                               0, 0);
        $self->add_to_output($fp);
        $switch = 0;
	    }
    }
  }
}




#############################
# subroutines
#############################
sub _read_fasta {
    local (*FILE) = @_;
    my ($id , $seq , %name2seq);
    while (<FILE>) {
        chomp;
        if (/^>(\S+)/) {
            my $new_id = $1;
            if ($id) {
                $name2seq{$id} = $seq;
            }
            $id = $new_id ; $seq = "" ;
        } 
        elsif (eof) {
            if ($id) {
                $seq .= $_ ;
                $name2seq{$id} = $seq;
            }
        }
        else {
            $seq .= $_ ;
        }
    }
    return %name2seq;
}


1;
