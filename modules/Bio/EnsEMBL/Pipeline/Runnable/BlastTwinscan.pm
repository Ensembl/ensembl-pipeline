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

Bio::EnsEMBL::Pipeline::Runnable::BlastTwinscan

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::BlastTwinscan->new
    ('-genomic'    => $genseq,
     '-features'   => $features,
     '-seqfetcher' => $seqfetcher
    );
    
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::BlastTwinscan;

use vars qw(@ISA);
use strict;

#use Bio::EnsEMBL::Pipeline::Runnable::Twinscan;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::PrimarySeqI;
use Bio::SeqIO;
use Bio::DB::RandomAccessI;
use Bio::EnsEMBL::Pipeline::Tools::BPlite; # NOT Bio::Tools:BPlite - it doesn't support Multi in bioperl 0.7
use Bio::EnsEMBL::Pipeline::Runnable::Twinscan;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_idlist'} = []; #create key to an array of feature pairs
  
  my( $genomic, $ids, $seqfetcher ) = $self->_rearrange([qw(GENOMIC
							    IDS
							    SEQFETCHER)],
							@args);
  
  $self->throw("No genomic sequence input")            unless defined($genomic);
  $self->throw("[$genomic] is not a Bio::PrimarySeqI") unless $genomic->isa("Bio::PrimarySeqI");
  $self->genomic_sequence($genomic) if defined($genomic);
  
  $self->throw("No seqfetcher provided")           
    unless defined($seqfetcher);
  $self->throw("[$seqfetcher] is not a Bio::DB::RandomAccessI") 
    unless $seqfetcher->isa("Bio::DB::RandomAccessI");
  $self->seqfetcher($seqfetcher) if defined($seqfetcher);
  
  if (defined($ids)) {
      $self->ids($ids);
  }
  else {
    $self->throw("Can't run BlastTwinscan without ids!\n");
  }
  
  return $self; # success - we hope!
}

=head2 ids

    Title   :   ids
    Usage   :   $self->ids($ids)
    Function:   Get/set method for ids of sequence to be blasted
    Returns :   array of strings(ids)
    Args    :   array of strings(ids)

=cut

sub ids {
  my ($self, $ids) = @_;
  if (defined($ids)) {
    if (ref($ids) eq "ARRAY") {
      push(@{$self->{'_idlist'}},@$ids);
    } else {
      $self->throw("[$ids] is not an array ref.");
    }
  }

  return @{$self->{'_idlist'}};

}

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::EnsEMBL::Pipeline::SeqFetcher object
    Args    :   Bio::EnsEMBL::Pipeline::SeqFetcher object

=cut

sub seqfetcher {
  my( $self, $value ) = @_;    
  if ($value) {
    #need to check if passed sequence is Bio::DB::RandomAccessI object
    $value->isa("Bio::DB::RandomAccessI") || $self->throw("Input isn't a Bio::DB::RandomAccessI");
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
}

=head2 conservation_sequence

    Title   :   conservation_sequence
    Usage   :   $self->conservation_sequence($conseq)
    Function:   Get/set method for conservation sequence needed by Twinscan
    Returns :   filename
    Args    :   filename

=cut

sub conservation_sequence {
  my( $self, $conseq ) = @_;    
  if ($conseq) {
    my $confile = "/tmp/conseq.$$.dat";
    open (OUT, ">$confile") or $self->throw("cannot open $confile to write conservation sequence:$!\n");
    print OUT "$conseq\n";
    close OUT or $self->throw("cannot close $confile after writing conservation sequence:$!\n");
    $self->{'_conseq'} = $confile;
  }
  return $self->{'_conseq'};
}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs est2genome on each distinct feature id
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;

    my @ids = $self->ids;

    $self->blast_and_conserve(@ids);

    my $conseq = $self->conservation_sequence;
    $self->throw("no conservation seqfile\n") unless -e $conseq;

    my $twin      = new Bio::EnsEMBL::Pipeline::Runnable::Twinscan('-genomic'    => $self->genomic_sequence,
								   '-conseq' => $conseq,
								  );

    $twin->run;
    
    my @f = $twin->output;

    #   foreach my $f (@f) {
    #      print(STDERR "PogAligned output is $f " . $f->seqname . " " . $f->start . "\t" . $f->end . "\t" . $f->score .  "\n");
    #   }
    
    push(@{$self->{'_output'}},@f);

    unlink $conseq;
}

sub blast_and_conserve {
    my ($self,@ids) = @_;

    # fetch sequences & index for blast
    my @seq         = $self->get_Sequences(@ids);
    my @valid_seq   = $self->validate_sequence(@seq);
    my $blastdb     = $self->make_blast_db(@valid_seq);

    # blast genomic vs sequences
    my $blastfile = $self->run_blast($self->genomic_sequence,$blastdb);

    # produce conservation sequence
    my $conseq = $self->make_conservation_seq($blastfile);
    $self->conservation_sequence($conseq);

    unlink $blastfile;
    unlink $blastdb;
    unlink $blastdb.".csq";
    unlink $blastdb.".nhd";
    unlink $blastdb.".ntb";
}

sub run_blast {

    my ($self,$seq,$db) = @_;

    my $blastout = $self->get_tmp_file("/tmp/","blast","swir.msptmp");
    my $seqfile  = $self->get_tmp_file("/tmp/","seq","fa");

    my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
				-file   => ">$seqfile");

    $seqio->write_seq($seq);
    close($seqio->_filehandle);

    my $command  = "wublastn $db $seqfile M=1 N=-1 Q=5 R=1 W=10 X=30 S=30 S2=30 gapS2=30 Z=3000000000 Y=3000000000 B=10000 V=100 -p=1 -warnings > $blastout";

    print (STDERR "Running command $command\n");
    my $status = system($command );

    print("Exit status of blast is $status\n");

    unlink $seqfile;    
    return $blastout;
}

# much of this is adapted directly from Michael & Paul's conseq_new.pl
sub make_conservation_seq {
  my ($self, $blast_file) = @_;

  my $genseq = $self->genomic_sequence->seq;
  my @hsp;
  my @gen_index; # genomic index - contains the conservation symbols
  my %Transition = (
		    A => {'G'=>1, 'R'=>1},
		    C => {'T'=>1, 'Y'=>1},
		    G => {'A'=>1, 'R'=>1},
		    T => {'C'=>1, 'Y'=>1},
		    R => {'A'=>1, 'G'=>1},
		    Y => {'C'=>1, 'T'=>1}
		   );
  
  # note that this order may get changed later
  my %Order = (
	       ':' => 1, # mismatch
	       '|' => 2, # match
	       '.' => 3, # unaligned
	       '/' => 4, # transition
	       '-' => 5, # gap
	      );

  # parse blast results 
  open(BLAST, $blast_file) or $self->throw("could not open $blast_file:$!\n");
  my $blast = new BPlite::Multi(\*BLAST);
  while(my $report = $blast->nextReport) {
    while(my $sbjct = $report->nextSbjct) {
      while(my $hsp = $sbjct->nextHSP) {
	push @hsp, $hsp;
      }
    }
  }
  close BLAST;

  @hsp = sort {$b->score <=> $a->score or $b->percent <=> $a->percent} @hsp;



  foreach my $hsp (@hsp) {
    my ($begin, $end) = ($hsp->qb, $hsp->qe);
    my ($qs, $hs, $as) = ($hsp->qa, $hsp->sa, $hsp->as);
    if ($begin > $end) {
      ($begin, $end) = ($end, $begin);
      $qs = reverse $qs;
      $hs = reverse $hs;
      $as = reverse $as;
    }
    
    my $conseq = "";
    for (my $i=0;$i<length($qs);$i++) {
      my $qb = substr($qs, $i, 1);
      my $hb = substr($hs, $i, 1);
      my $ab = substr($as, $i, 1);
      if    ($qb eq '-') {next}           # skip query gaps, no indexing
      elsif ($hb eq '-') {$conseq .= '-'} # sbjct gap
      elsif ($ab eq ' ') {
	if (exists $Transition{$qb}{$hb}) {$conseq .= '/'}
	else                              {$conseq .= ':'}
      } # mismatch
      else {$conseq .= '|'} # match
    }

    for (my $i=0;$i<length($conseq);$i++) {
      my $index = $i + $begin -1;
      my $prev  = $gen_index[$index];
      my $char = substr($conseq, $i, 1);
      
      # precedence rule - if defined by a better HSP, leave it alone
      if (not defined $prev) {$gen_index[$index] = $char}
    }
  }

  # change undefined values to the unaligned symbol
  for(my $i=0;$i<length($genseq);$i++) {
    $gen_index[$i] = '.' unless defined $gen_index[$i]; # unaligned symbol is .
  }

  # sanity check
  if (length($genseq) != scalar(@gen_index)) {
	$self->throw( "conseq length difference! " . length($genseq) . " != " . scalar(@gen_index) . "\n");
      }

  # output
  my $conseq = join("", @gen_index);
  $conseq =~ tr/\//:/; # change transition to mismatch
  $conseq =~ tr/\-/:/; # change gaps       to mismatch

  my %symbol;
  for(my $i=0;$i<length($conseq);$i++) {
    $symbol{substr($conseq,$i,1)}++;
  }
  
  my $n = 0;

  # convert $conseq from symbols to numbers 
  foreach my $char (sort {$Order{$a} <=> $Order{$b}} keys %symbol) {
    $symbol{$char} = $n;
    $n++;
  }

  my @key   = sort {$Order{$a} <=> $Order{$b}}    keys %symbol;
  my @value = sort {$a <=> $b} values %symbol;
  my $symbols  = join("", @key);
  my $numbers  = join("", @value);
  
  my $code = "\$conseq =~ tr[$symbols][$numbers]";
  eval "$code";
  return $conseq;
}


sub print_FeaturePair {
    my ($self,$pair) = @_;

    print STDERR $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . $pair->score . "\t" .
	$pair->strand . "\t" . $pair->hseqname . "\t" . $pair->hstart . "\t" . $pair->hend . "\t" . $pair->hstrand . "\n";
}

sub make_blast_db {
    my ($self,@seq) = @_;

    my $blastfile = $self->get_tmp_file('/tmp/','blast','fa');
    my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
			       -file   => ">$blastfile");

    print STDERR "Blast db file is $blastfile\n";

    foreach my $seq (@seq) {
	print STDERR "Writing seq " . $seq->id ."\n";
	$seqio->write_seq($seq);
    }

    close($seqio->_filehandle);

    my $status = system("pressdb $blastfile");
    print (STDERR "Status from pressdb $status\n");

    return $blastfile;
}

sub get_tmp_file {
    my ($self,$dir,$stub,$ext) = @_;

    
    if ($dir !~ /\/$/) {
	$dir = $dir . "/";
    }

#    $self->check_disk_space($dir);

    my $num = int(rand(10000));
    my $file = $dir . $stub . "." . $num . "." . $ext;

    while (-e $file) {
	$num = int(rand(10000));
	$file = $stub . "." . $num . "." . $ext;
    }			
    
    return $file;
}
    
sub get_Sequences {
    my ($self,@ids) = @_;

    my @seq;

    foreach my $id (@ids) {
	my $seq = $self->get_Sequence($id);

	if (defined($seq) && $seq->length > 0) {
	    push(@seq,$seq);
	} else {
	    print STDERR "Invalid sequence for $id - skipping\n";
	}
    }

    return @seq;

}

sub validate_sequence {
    my ($self,@seq) = @_;
    my @validated;
    foreach my $seq (@seq)
    {
        print STDERR ("mrna feature $seq is not a Bio::PrimarySeq or Bio::Seq\n") 
                                    unless ($seq->isa("Bio::PrimarySeq") ||
                                            $seq->isa("Bio::Seq"));
        my $sequence = $seq->seq;
        if ($sequence !~ /[^acgtn]/i)
        {
            push (@validated, $seq);
        }
        else 
        {
            $_ = $sequence;
            my $len = length ($_);
            my $invalidCharCount = tr/bB/xX/;

            if ($invalidCharCount / $len > 0.05)
            {
                $self->warn("Ignoring ".$seq->display_id()
                    ." contains more than 5% ($invalidCharCount) "
                    ."odd nucleotide codes ($sequence)\n Type returns "
                    .$seq->moltype().")\n");
            }
            else
            {
                $self->warn ("Cleaned up ".$seq->display_id
                   ." for blast : $invalidCharCount invalid chars \n");
                $seq->seq($_);
                push (@validated, $seq);
            }
        }
    } 
    return @validated;  
}

=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id)
  Function: Fetches sequences with id $id
  Returns : Bio::PrimarySeq
  Args    : none

=cut
    
sub get_Sequence {
    my ($self,$id) = @_;
    my $seqfetcher = $self->seqfetcher;
    my $seq;

    if (!defined($id)) {
      $self->warn("No id input to get_Sequence");
    }  
    
    print(STDERR "Sequence id :  is [$id]\n");

    eval {
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };

    if(!defined($seq) && $@){
      $self->warn("Could not retrieve sequence for [$id]:\n");
    }

    return $seq;
}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
    my ($self) = @_;
    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }
    return @{$self->{'_output'}};
}


sub trim {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{'_trim'} = $arg;
  }
  return $self->{'_trim'};
}

1;


