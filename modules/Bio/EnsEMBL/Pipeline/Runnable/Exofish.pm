# Cared for by Kevin Howe  <klh@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Exofish

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Exofish;


use vars qw(@ISA);
use strict;
# Object preamble

use FileHandle;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Tools::BPlite;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);    
  
  # Now parse the input options and store them in the object
  my( $query, $program, $database, $options, $do_not_project) = 
      $self->_rearrange([qw(QUERY 
                            PROGRAM 
                            DATABASE 
                            OPTIONS
                            DONOTPROJECT)],
                        @args);
  
  if ($query) {
    $self->query($query);
  } else {
    $self->throw("No query sequence input.");
  }
    
  if ($database) {
    $self->database($database);
  } else {
    $self->throw("No database input");
  }

  
  $program = "wutblastx" if not $program;
  $self->program($self->find_executable($program));
  
  my $core_options = "-cpus 1 -compat1.4 -lcfilter -matrix EXOFISH -sort_by_highscore W=5 X=25 T=75 S=89 S2=89";  
  $core_options .= $options;
  $self->options($core_options);

  if (defined($do_not_project)) {
    $self->do_not_project($do_not_project);
  }

  return $self; # success - we hope!
}


###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
  my ($self, $dir) = @_;
  
  $self->workdir('/tmp') unless ($self->workdir($dir));
  $self->checkdir();

  # with -lcfilter, this hardmasks repeat sequence; speeds up search consierably  
  my $seq = $self->query->seq;    
  $seq =~ tr/N/n/;
  $self->query->seq($seq);
  
  # with EXOFISH parameters, tblastx only works if there 
  # is a string of at least 6aa (18nuc) non-Ns somewhere in the sequence 
  # (otherwise minimum score cannot be achieved in any frame)
  
  if ($seq =~ /[ACGT]{18}/) {
    $self->writefile(); 
    $self->run_analysis();
    $self->parse_results;
    $self->deletefiles();
  }    
  
  return 1;
}


=head2 run_analysis

    Title   :   run_analysis
    Usage   :   $obj->run_analysis
    Function:   Runs the blast query
    Returns :   nothing
    Args    :   none

=cut

sub run_analysis {
  my ($self) = @_;
      
  my $database = $self->database;
  my $filename = $self->filename;
  my $command = $self->program . " $database $filename " . $self->options . " 2>&1 > " . $self->results;

  #print STDERR "Running '$command'\n";
  # about to create the output file, so register it for clean-up
  $self->file($self->results);

  open(my $fh, "$command |") || $self->throw("Error opening Blast cmd <$command>.\n");

  # read the STDERR; only report if the message is FATAL  
  while(<$fh>){
    if(/FATAL:(.+)/){
      $self->throw($1);
    }
  }        #print STDERR "-------------- Here ends Blast's STDERR ---------------\n";
  unless(close $fh){
    $self->throw("Error running Blast cmd <$command>. Returned error $? BLAST EXIT: '"
                 . ($? >> 8) 
                 . "', SIGNAL '" 
                 . ($? & 127) 
                 . "', There was " 
                 . ($? & 128 ? 'a' : 'no') 
                 . " core dump");
  }

}


sub parse_results {
  my ($self,$fh) = @_;
  
  my $parser;
  if ($fh) {
    $parser = new Bio::EnsEMBL::Pipeline::Tools::BPlite(-fh => $fh);
  }
  else {
    $parser = new Bio::EnsEMBL::Pipeline::Tools::BPlite(-file => $self->results);
  }

  my @features;

 NAME: 
  while (my $sbjct = $parser->nextSbjct) {
    
  HSP: 
    while( my $hsp = $sbjct->nextHSP) {
      my ($qname) = $hsp->query->seqname =~ /^(\S+)/;
      my ($tname) = $hsp->subject->seqname =~ /^(\S+)/;
      
      my ($qstart, $qend, $qstrand) = ($hsp->query->start, $hsp->query->end, $hsp->query->strand);
      my ($tstart, $tend, $tstrand) = ($hsp->subject->start, $hsp->subject->end, $hsp->subject->strand);
      
      my $pid = $hsp->percent;
      my $score = $hsp->score;
      
      my $aalen = ($qend - $qstart + 1) / 3;
      #print "$qstart $qend $pid $aalen\n";      
      # Exofish filtering strategy. See figure 1 in Roest Crollius et. al. 2000
      
      if ($aalen < 13) {
        next;
      } 
      elsif ($aalen < 15) {
        next if $pid < 95.0;
      } 
      elsif ($aalen < 17) {
        next if $pid < 90.0;
      }
      elsif ($aalen < 19) {
        next if $pid < 85.0;
      }
      elsif ($aalen < 21) {
        next if $pid < 80.0;
      }
      elsif ($aalen < 22) {
        next if $pid < 77.5;
      }
      elsif ($aalen < 23) {
        next if $pid < 75.0;
      }
      elsif ($aalen < 24) {
        next if $pid < 70.0;
      }
      elsif ($aalen < 25) {
        next if $pid < 67.5;
      }
      else {
        # length is 25 amino acids or more. 
        next if $pid < 55.0;
      }
      
      my $fp = Bio::EnsEMBL::FeaturePair->new();
      $fp->seqname($qname);
      $fp->start($qstart);
      $fp->end($qend);
      $fp->strand($qstrand);
      $fp->hseqname($tname);
      $fp->hstart($tstart);
      $fp->hend($tend);
      $fp->hstrand($tstrand);
      $fp->score($score);
      $fp->percent_id($pid);
      
      my $df = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => [$fp]);
      push @features, $df;
    }
  }

  if (not $self->do_not_project) {
    my (@projected_hits);
    
    foreach my $fp (sort {$a->start <=> $b->start} @features) {
      my $reg = {  start    => $fp->start, 
                   end      => $fp->end,
                   score    => [$fp->score],
                   hseqname => [$fp->hseqname],
                   htstart  => [$fp->hstart],
                   hend     => [$fp->hend],

                 };
      
      if (not @projected_hits or $projected_hits[-1]->{end} + 1 < $reg->{start}) {
        push @projected_hits, $reg;
      }
      else {
        # overlap
        push @{$projected_hits[-1]->{hseqname}}, $reg->{hseqname};
        push @{$projected_hits[-1]->{hstart}}, $reg->{hstart};
        push @{$projected_hits[-1]->{hend}}, $reg->{hend};
        push @{$projected_hits[-1]->{score}}, $reg->{score};

        if ($reg->{end} > $projected_hits[-1]->{end}) {
          $projected_hits[-1]->{end} = $reg->{end};
        }	    
      }
    }
    
    @features = ();
    
    foreach my $reg (@projected_hits) {      
      my $fp = Bio::EnsEMBL::FeaturePair->new();
      $fp->seqname($self->query->id);
      $fp->start($reg->{start});
      $fp->end($reg->{end});
      $fp->strand(1);

      my $tot_score = 0.0;
      foreach my $el (@{$reg->{score}}) {
        $tot_score += $el;
      }
      $fp->score($tot_score);

      if (@{$reg->{hseqname}} > 1) {      
        $fp->hseqname(sprintf("Many(%d)", scalar(@{$reg->{hseqname}})));
        $fp->hstart(1);
        $fp->hend($reg->{end} - $reg->{start} + 1);
      }
      else {
        $fp->hseqname($reg->{hseqname}->[0]);
        $fp->hstart($reg->{hstart}->[0]);
        $fp->hend($reg->{hend}->[0]);
      }

      push @features,  Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => [$fp]);
    }    
  }
  
  $self->output(@features);
  
  return $self->output; 
}


################
# get/set
################

sub database {
  my ($self,$database) = @_;

  if ($database) {    
    $self->{'_database'} = $database;
  }

  return $self->{'_database'};
}


sub program {
  my ($self,$program) = @_;

  if ($program) {    
    $self->{'_program'} = $program;
  }

  return $self->{'_program'};
}



sub query {
  my ($self, $seq) = @_;
  if ($seq) {
    unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI"))
    {
      $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
    }
    $self->{'_query'} = $seq ;
    $self->filename($seq->id.".$$.seq");
    $self->results($self->filename.".exofish");
  }
  return $self->{'_query'};
}


sub output {
  my ($self, @arg) = @_;

  if (!defined($self->{'_fplist'})) {
    $self->{'_fplist'} = [];
  }
  
  if (@arg) {
    @{$self->{'_fplist'}} = @arg;
  }

  return @{$self->{'_fplist'}};
}


sub do_not_project {
  my ($self, $do_not_project) = @_;

  if (defined($do_not_project)) {
    $self->{'_do_not_project'} = $do_not_project;
  }
  
  return $self->{'_do_not_project'};
}


1;

