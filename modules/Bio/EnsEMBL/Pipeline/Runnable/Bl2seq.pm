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
							    -min_eval => $min_eval,
							    -min_score => $min_score,
							    -workdir => $workingdir,
							    -results => $results);
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new().

    $seq1 and $seq2 must be Bio::PrimarySeq object.

    $program (optional) must be a sting which locates bl2seq executable.

    $alntype (optional) must be a string which defines the alignment type, blastp, 
              blastn, blastx, tblastx or tblastn (default = blastn).

    $min_eval (optional) must be a floating, value for minimum Eval for a HSP.

    $min_score (optional) must be a interger, valu for minimum score (bits) for a HSP.

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
 min_eval,
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
use vars qw(@ISA);

# Object preamble - inherits from Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

  $self->{'_seq1'} = undef; # location of Bio::Seq object "-i" option 
  $self->{'_seq2'} = undef; # location of Bio::Seq object "-j" option
  $self->{'_program'} = "bl2seq"; # location of bl2seq executable
  $self->{'_alntype'} = "blastn"; # type of alignment "-p" option
  $self->{'_min_score'} = 40; # value for minimum score 
  $self->{'_min_eval'} = 0.01; # value for minimum E value "-e" option
  $self->{'_fplist'} = []; # an array of feature pairs (the output)
  $self->{'_workdir'} = "/tmp"; # location of temp directory
  $self->{'_results'} = $self->{'_workdir'}."/results.".$$; # location of result file

  my ($seq1, $seq2, $program, $alntype, $min_score, $min_eval, $workdir, $results) = 
    $self->_rearrange([qw(SEQ1 SEQ2 PROGRAM ALNTYPE MIN_SCORE MIN_EVAL WORKDIR RESULTS)], @args);

  $self->throw("Must pass in both seq1 and seq1 args") if (! defined $seq1 || ! defined $seq2);
  $self->seq1($seq1);
  $self->seq2($seq2);

  $self->program($self->find_executable($program)) if ($program);
  $self->alntype($alntype) if ($alntype);
  if ($workdir) {
    $self->workdir($workdir);
    $self->results($self->workdir."/results.".$$);
  }
  $self->results($results) if ($results);
  $self->min_score($min_score) if (defined $min_score);
  $self->min_eval($min_eval) if (defined $min_eval);
  
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

=head2 min_eval

 Title   : min_eval
 Usage   : $obj->min_eval($newval)
 Function: Get/Set the min_eval value
 Example : 
 Returns : value of -e option used by bl2seq
 Args    : floating (optional)


=cut

sub min_eval {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_min_eval'} = $value;
  }
  return $self->{'_min_eval'};
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
  $self->writefile('seq1');
  $self->filename($sbjct);
  $self->writefile('seq2');

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
  my ($g,$W,$G,$E,$X) = qw(T 10 -1 3 10);
  print STDERR ("Running bl2seq\n" . $self->program .
		                     " -i $query" .
		                     " -j $sbjct" .
		                     " -g $g" .
		                     " -W $W" .
		                     " -G $G" .
		                     " -E $E" .
		                     " -X $X" .
		                     " -p " . $self->alntype .
		                     " -e " . $self->min_eval . " > " .
		                     $self->results. "\n");

  $self->throw("Failed during bl2seq run, $!\n") unless (system ($self->program .
								 " -i $query" .
								 " -j $sbjct" .
								 " -g $g" .
								 " -W $W" .
								 " -G $G" .
								 " -E $E" .
								 " -X $X" .
								 " -p " . $self->alntype .
								 " -e " . $self->min_eval . " > " .
								 $self->results) == 0);
}  
  
sub parse_results { 
  my ($self) = @_;
  
  open BL2SEQ, $self->results || 
    $self->throw("Coudn't open file ".$self->results.", $!\n");
  my $filehandle = \*BL2SEQ;

  my $bl2seq_parsing = Bl2seq::Parser->new('-fh' => $filehandle,
					   '-alntype' => 'blastn',
					   '-min_score' => $self->min_score,
					   '-qname' => $self->seq1->id,
					   '-sname' => $self->seq2->id);

  while (my $DnaDnaAlignFeature = $bl2seq_parsing->nextHSP) {
    my @ungapped_features = $DnaDnaAlignFeature->ungapped_features;
    next unless (scalar @ungapped_features);
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

package Bl2seq::Parser;

use strict;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Carp;
use vars qw(@ISA);

# Object preamble - inherits from Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

  $self->{'_fh'} = undef; # filehandle on results file
  $self->{'_alntype'} = undef; # type of blast alignment
  $self->{'_lastline'} = undef; # last line being read in the filehandle
  $self->{'_min_score'} = 0; # last line being read in the filehandle
  $self->{'_qname'} = undef; # filehandle on results file
  $self->{'_sname'} = undef; # filehandle on results file

  my ($fh,$alntype,$min_score,$qname,$sname) = $self->_rearrange([qw(FH ALNTYPE MIN_SCORE QNAME SNAME)], @args);

  $self->throw("Must pass in both fh and alntype args") if (! defined $fh ||
							    ! defined $alntype ||
							    ! defined $qname ||
							    ! defined $sname);
  $self->fh($fh);
  $self->alntype($alntype);
  $self->min_score($min_score) if (defined $min_score);
  $self->qname($qname);
  $self->sname($sname);

  return $self;
}

sub fh {
  my ($self,$value) = @_;
  if(defined $value) {
    if (ref($value) eq "GLOB") {
      $self->{'_fh'} = $value;
    } else {
      $self->throw("Arg for fh method should be a filehandle\n");
    }
  }
  return $self->{'_fh'};
}

sub alntype {
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'_alntype'} = $value;
  }
  return $self->{'_alntype'};
}

sub lastline {
  my ($self,$value) = @_;
  if(defined $value) {
    $self->{'_lastline'} = $value;
  }
  return $self->{'_lastline'};
}

sub min_score {
  my ($self,$value) = @_;
  if(defined $value) {
    $self->{'_min_score'} = $value;
  }
  return $self->{'_min_score'};
}

sub qlength {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_qlength'} = $value;
  }
  return $self->{'_qlength'};

}

sub slength {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_slength'} = $value;
  }
  return $self->{'_slength'};
}

sub qname {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_qname'} = $value;
  }
  return $self->{'_qname'};

}

sub sname {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_sname'} = $value;
  }
  return $self->{'_sname'};
}

sub nextHSP {
  my ($self) = @_;
  my $fh = $self->fh;
  unless (defined $self->lastline) {

    while (defined (my $line = <$fh>)) {
      next if ($line !~ /^\s+\(\d+,*\d*\s+letters\).*$/ &&
	       $line !~ /^\s+Length\s+=\s+\d+.*$/ &&
	       $line =~ /^$/);
      if ($line =~ /^\s+\((\S+)\s+letters\).*$/) {
	my $qlength = $1;
	$qlength =~ s/[^0-9]//g;
	$self->qlength($qlength);
      }
      if ($line =~ /^\s+Length\s+=\s+(\d+).*$/) {
	$self->slength($1);
	$self->lastline($line);
	last;
      }
    }
  }
  my $parse_alntype = "parse_".$self->alntype;
  my $DnaDnaAlignFeature = $self->$parse_alntype($fh,$self->qlength,$self->slength);
  return $DnaDnaAlignFeature
}

sub parse_blastp {
  my ($self) = @_;
  $self->throw("Parsing of 'blastp' bl2seq output not implemented yet\n")
}

sub parse_blastn {
  my ($self,$fh,$qlength,$slength) = @_;
  my ($bits,$score,$expect,$identity,$qhspstart,$qhspend,$qstrand,$shspstart,$shspend,$sstrand);
  my ($qhspseq,$shspseq) = ("","");
  my @blocks = ();
  if (defined $self->lastline) {
    my $line = $self->lastline;
    if ($line =~ /^\s+Score\s+=\s+(\d+\.*\d*)\s+bits\s+\((\d+)\),\s+Expect\s+=\s+(\S+).*$/) {
      ($bits,$score,$expect) = ($1,$2,$3);
    }
  }
  while (defined (my $line = <$fh>)) {
    next if ($line =~ /^$/);
    if ($line =~ /\s+Lambda\s+.*$/) {
      last;
    }
    if ($line =~ /^\s+Score\s+=\s+.*$/ && defined $qhspstart) {
      if ($bits >= $self->min_score) {
	$self->lastline($line);
	last;
      } else {
	$qhspstart = undef;
	$shspstart = undef;
	$qhspseq = "";
	$shspseq = "";
      }
    } 
    if ($line =~ /^\s+Score\s+=\s+(\d+\.*\d*)\s+bits\s+\((\d+)\),\s+Expect\s+=\s+(\S+).*$/) {
      ($bits,$score,$expect) = ($1,$2,$3);
    }
    if ($line =~ /^\s+Identities\s+=\s+\d+\/\d+\s+\((\d+)%\).*$/) {
      $identity = $1;
    }
    if ($line =~ /^\s+Strand\s+=\s+(\S+)\s+\/\s+(\S+)$/) {
      ($qstrand,$sstrand) = ($1,$2);
      if ($qstrand eq "Plus") {
	$qstrand = 1;
      } elsif ($qstrand eq "Minus") {
	$qstrand = -1;
      }
      if ($sstrand eq "Plus") {
	$sstrand = 1;
      } elsif ($sstrand eq "Minus") {
	$sstrand = -1;
      }
    }
    if ($line =~ /^Query:\s+(\d+)\s+(\S+)\s+(\d+)$/) {
      my ($start,$seq,$end) = ($1,$2,$3);
      $qhspstart = $start unless (defined $qhspstart);
      $qhspend = $end;
      $qhspseq .= $seq;
    }	
    if ($line =~ /^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)$/) {
      my ($start,$seq,$end) = ($1,$2,$3);
      $shspstart = $start unless (defined $shspstart);
      $shspend = $end;
      $shspseq .= $seq;
    }
  }

  unless (length($qhspseq) == length($shspseq)) {
    $self->throw("Something goes wrong with the bl2seq parsing\n");
  }
  
  return @blocks unless (defined $qhspstart &&
			 defined $qhspend &&
			 defined $qstrand &&
			 defined $shspstart &&
			 defined $shspend &&
			 defined $sstrand);

  return @blocks if ($bits <= $self->min_score);

  if ($qstrand == -1) {
    my $tmp = $qhspstart;
    $qhspstart = $qhspend;
    $qhspend = $tmp;
  }
  if ($sstrand == -1) {
    my $tmp = $shspstart;
    $shspstart = $shspend;
    $shspend = $tmp;
  }
  my $block = Bl2seq::Block->new();
  $block->bits($bits);
  $block->score($score);
  $block->expect($expect);
  $block->identity($identity);
  $block->qstart($qhspstart);
  $block->qend($qhspend);
  $block->qstrand($qstrand);
  $block->qlength($qlength);
  $block->sstart($shspstart);
  $block->send($shspend);
  $block->sstrand($sstrand);
  $block->slength($slength);
  $block->qseq($qhspseq);
  $block->sseq($shspseq);

  if (abs($block->qend - $block->qstart) == abs($block->send - $block->sstart) &&
      $block->qseq !~ /[acgtnACGTN]+/ &&
      $block->sseq !~ /[acgtnACGTN]+/) {
    push @blocks, $block;
  } else {
    while (my $ungapped_block = $block->nextUngappedBlock) {
      push @blocks, $ungapped_block;
    }
  }
  my @ungapped_features;
  foreach my $ungapped_block (@blocks) {
    my ($qstart,$qend,$qstrand,$sstart,$send,$sstrand,$bits,$perc_id) = ($ungapped_block->qstart,$ungapped_block->qend,$ungapped_block->qstrand,$ungapped_block->sstart,$ungapped_block->send,$ungapped_block->sstrand,$ungapped_block->bits,$ungapped_block->identity);

    my $fp = Bio::EnsEMBL::FeatureFactory->new_feature_pair();

    $fp->start($qstart);
    $fp->end($qend);
    $fp->strand($qstrand);
    $fp->seqname($self->qname);
    $fp->hstart($sstart);
    $fp->hend($send);
    $fp->hstrand($sstrand);
    $fp->hseqname($self->sname);
    $fp->score($bits);
    $fp->percent_id($perc_id);
    
    push @ungapped_features, $fp;
  }
  my $DnaDnaAlignFeature = new Bio::EnsEMBL::DnaDnaAlignFeature('-features' => \@ungapped_features);
  return $DnaDnaAlignFeature;
}

sub parse_blastx {
  my ($self) = @_;
  $self->throw("Parsing of 'blastx' bl2seq output not implemented yet\n")
}

sub parse_tblastx {
  my ($self) = @_;
  $self->throw("Parsing of 'tblastx' bl2seq output not implemented yet\n")
}

sub parse_tblastn {
  my ($self) = @_;
  $self->throw("Parsing of 'tblastn' bl2seq output not implemented yet\n")
}

package Bl2seq::Block;

use strict;
use Carp;
use vars qw(@ISA);

# Object preamble - inherits from Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  
  $self->{'_qlength'} = undef; # length of the query sequence
  $self->{'_slength'} = undef; # length of the subject sequence

  my ($qlength, $slength) = $self->_rearrange([qw(QLENGTH SLENGTH)], @args);

  $self->qlength($qlength) if ($qlength);
  $self->slength($slength) if ($slength);

  return $self;

}

sub bits {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_bits'} = $value;
  }
  return $self->{'_bits'};
}

sub score {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_score'} = $value;
  }
  return $self->{'_score'};
} 

sub expect {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_expect'} = $value;
  }
  return $self->{'_expect'};
}

sub identity {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_identity'} = $value;
  }
  return $self->{'_identity'};
}

sub qstart {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_qstart'} = $value;
  }
  return $self->{'_qstart'};
}

sub qend {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_qend'} = $value;
  }
  return $self->{'_qend'};
}

sub qstrand {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_qstrand'} = $value;
  }
  return $self->{'_qstrand'};
}

sub sstart {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_sstart'} = $value;
  }
  return $self->{'_sstart'};
}

sub send {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_send'} = $value;
  }
  return $self->{'_send'};
}

sub sstrand {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_sstrand'} = $value;
  }
  return $self->{'_sstrand'};
}

sub qlength {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_qlength'} = $value;
  }
  return $self->{'_qlength'};
}

sub slength {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_slength'} = $value;
  }
  return $self->{'_slength'};
}

sub qseq {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_qseq'} = $value;
  }
  return $self->{'_qseq'};
}

sub sseq {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_sseq'} = $value;
  }
  return $self->{'_sseq'};
}

sub nextUngappedBlock {
  my ($self) = @_;
  my $block_initialized = 0;
  my ($newqseq,$newsseq) = ("","");
  my $block = Bl2seq::Block->new();

# Assigning here to each ungapped block the bits, score, expect and identity
# of the original gapped block ($self). Not very conventional...

  $block->bits($self->bits);
  $block->score($self->score);
  $block->expect($self->expect);
  $block->identity($self->identity);

  for (my $i = 0; $i <= length($self->qseq); $i++) {

    if (substr($self->qseq,$i,1) eq "-" ||
	substr($self->sseq,$i,1) eq "-" ||
	$i == length($self->qseq)) {

      if (substr($self->sseq,$i,1) eq "-" && ! $block_initialized) {
	if ($self->qstrand == 1) {
	  $self->qstart($self->qstart + 1);
	} elsif ($self->qstrand == -1) {
	  $self->qend($self->qend - 1);
	}
	$self->qseq(substr($self->qseq,$i + 1));
	$self->sseq(substr($self->sseq,$i + 1));
	$i = -1;
	next;
      }
      if (substr($self->qseq,$i,1) eq "-" && ! $block_initialized) {
	if ($self->sstrand == 1) {
	  $self->sstart($self->sstart + 1)
	} elsif ($self->sstrand == -1) {
	  $self->send($self->send - 1);
	}
	$self->qseq(substr($self->qseq,$i + 1));
	$self->sseq(substr($self->sseq,$i + 1));
	$i = -1;
	next;
      }

      next unless ($block_initialized);
      
      $self->qseq(substr($self->qseq,$i));
      $self->sseq(substr($self->sseq,$i));
      

      $newqseq .= substr($self->qseq,$i,1) if ($i == length($self->qseq));
      $newsseq .= substr($self->sseq,$i,1) if ($i == length($self->qseq));

      if ($self->qstrand == 1) {
	$block->qend($self->qstart + $i - 1);
      } elsif ($self->qstrand == -1) {
	$block->qstart($self->qend - $i + 1);
      }
      if ($self->sstrand == 1) {
	$block->send($self->sstart + $i - 1);
      } elsif ($self->sstrand == -1) {
	$block->sstart($self->send - $i + 1);
      }
      
      unless ($i == length($self->qseq)) {
	if ($self->qstrand == 1) {
	  $self->qstart($block->qend + 1);
	} elsif ($self->qstrand == -1) {
	  $self->qend($block->qstart - 1);
	} 
	if ($self->sstrand == 1) {
	  $self->sstart($block->send + 1)
	} elsif ($self->sstrand == -1) {
	  $self->send($block->sstart - 1)
	} 
      }

      $block->qstrand($self->qstrand);
      $block->qlength($self->qlength);


      $block->sstrand($self->sstrand);
      $block->slength($self->slength);

      $block->qseq($newqseq);
      $block->sseq($newsseq);
      
      return $block;

    } elsif (! $block_initialized) {

      $newqseq .= substr($self->qseq,$i,1);
      $newsseq .= substr($self->sseq,$i,1);

      if ($self->qstrand == 1) {
	$block->qstart($self->qstart + $i);
      } elsif ($self->qstrand == -1) {
	$block->qend($self->qend - $i);
      }
      if ($self->sstrand == 1) {
	$block->sstart($self->sstart + $i)
      } elsif ($self->sstrand == -1) {
	$block->send($self->send - $i);
      }
      $block_initialized = 1;

    } else {

      $newqseq .= substr($self->qseq,$i,1);
      $newsseq .= substr($self->sseq,$i,1);

    }
  }
}

1;
