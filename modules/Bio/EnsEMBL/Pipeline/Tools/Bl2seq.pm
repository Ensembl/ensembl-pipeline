
package Bio::EnsEMBL::Pipeline::Tools::Bl2seq;

use strict;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Pipeline::Tools::Block;
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
#  print STDERR "parse_alntype: $parse_alntype\n";
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
  my $block = Bio::EnsEMBL::Pipeline::Tools::Block->new();
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
    while (my $ungapped_block = $block->nextUngappedBlock("blastn")) {
      push @blocks, $ungapped_block;
    }
  }
  my @ungapped_features;
  foreach my $ungapped_block (@blocks) {
    my ($qstart,$qend,$qstrand,$sstart,$send,$sstrand,$bits,$perc_id) = ($ungapped_block->qstart,$ungapped_block->qend,$ungapped_block->qstrand,$ungapped_block->sstart,$ungapped_block->send,$ungapped_block->sstrand,$ungapped_block->bits,$ungapped_block->identity);

    my $fp = new Bio::EnsEMBL::FeaturePair;
#    my $seqone = new Bio::EnsEMBL::SeqFeature;
#    my $seqtwo = new Bio::EnsEMBL::SeqFeature;
     
#    $fp->feature1($seqone);
#    $fp->feature2($seqtwo);

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
  # set here only the score, bits, percent_id,expect
  return $DnaDnaAlignFeature;
}

sub parse_blastx {
  my ($self) = @_;
  $self->throw("Parsing of 'blastx' bl2seq output not implemented yet\n")
}

sub parse_tblastx {
#  my ($self) = @_;
#  $self->throw("Parsing of 'tblastx' bl2seq output not implemented yet\n")
  my ($self,$fh,$qlength,$slength) = @_;
  my ($bits,$score,$expect,$identity,$positivity,$qhspstart,$qhspend,$qstrand,$shspstart,$shspend,$sstrand);
  my ($qhspseq,$shspseq) = ("","");
  my @blocks = ();
  if (defined $self->lastline) {
    my $line = $self->lastline;
    if ($line =~ /^\s+Score\s+=\s+(\d+\.*\d*)\s+bits\s+\((\d+)\),\s+Expect\(*\S*\)*\s+=\s+(\S+).*$/) {
      ($bits,$score,$expect) = ($1,$2,$3);
    }
  }
  while (defined (my $line = <$fh>)) {
    next if ($line =~ /^$/);
    if ($line =~ /\s+Lambda\s+.*$/) {
      last;
    }
    if ($line =~ /^\s+Score\s+=\s+.*$/ && defined $qhspstart) {
#      print STDERR "bits: $bits\n";
#      print STDERR "min_score: ",$self->min_score,"\n";
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
    if ($line =~ /^\s+Score\s+=\s+(\d+\.*\d*)\s+bits\s+\((\d+)\),\s+Expect\(*\S*\)*\s+=\s+(\S+).*$/) {
      ($bits,$score,$expect) = ($1,$2,$3);
    }
    if ($line =~ /^\s+Identities\s+=\s+\d+\/\d+\s+\((\d+)%\),\s+Positives\s+=\s+\d+\/\d+\s+\((\d+)%\).*$/) {
      ($identity,$positivity) = ($1,$2);
    }
    if ($line =~ /^\s+Frame\s+=\s+([+-][123])\s+\/\s+([+-][123])$/) {
      ($qstrand,$sstrand) = ($1,$2);
      if ($qstrand > 0) {
	$qstrand = 1;
      } elsif ($qstrand < 0) {
	$qstrand = -1;
      }
      if ($sstrand > 0) {
	$sstrand = 1;
      } elsif ($sstrand < 0) {
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
  my $block = Bio::EnsEMBL::Pipeline::Tools::Block->new();
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
      $block->qseq !~ /[a-ik-np-tv-zA-IK-NP-TV-Z*]+/ &&
      $block->sseq !~ /[a-ik-np-tv-zA-IK-NP-TV-Z*]+/) {
    push @blocks, $block;
  } else {
    while (my $ungapped_block = $block->nextUngappedBlock("tblastx")) {
      push @blocks, $ungapped_block;
    }
  }
  my @ungapped_features;
  foreach my $ungapped_block (@blocks) {
    my ($qstart,$qend,$qstrand,$sstart,$send,$sstrand,$bits,$perc_id) = ($ungapped_block->qstart,$ungapped_block->qend,$ungapped_block->qstrand,$ungapped_block->sstart,$ungapped_block->send,$ungapped_block->sstrand,$ungapped_block->bits,$ungapped_block->identity);

    my $fp = new Bio::EnsEMBL::FeaturePair;
    my $seqone = new Bio::EnsEMBL::SeqFeature;
    my $seqtwo = new Bio::EnsEMBL::SeqFeature;
     
    $fp->feature1($seqone);
    $fp->feature2($seqtwo);


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

sub parse_tblastn {
  my ($self) = @_;
  $self->throw("Parsing of 'tblastn' bl2seq output not implemented yet\n")
}

1;
