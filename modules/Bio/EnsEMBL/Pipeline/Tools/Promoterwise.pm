package Bio::EnsEMBL::Pipeline::Tools::Promoterwise;

use vars qw(@ISA);
use Bio::EnsEMBL::DnaDnaAlignFeature;
use strict;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;
  
  $self->{'_fh'} = undef; # filehandle on results file
  $self->{'_file'} = undef; # path for a results file
  $self->{'_query_id'} = undef; # query id 
  $self->{'_target_id'} = undef; # target id
  $self->{'_eof'} = 0; # indicate if end of file and fh closed
  $self->{'_parsing_initialized'} = 0;
  $self->{'_score_in_mind'} = undef;
  $self->{'_command_line'} = "";
  $self->{'_matrix'} = "";
  $self->{'_options'} = "";
  
  
  my ($fh, $file, $query_id, $target_id) = $self->_rearrange([qw(FH FILE QUERY_ID TARGET_ID)], @args);

  $self->query_id($query_id) if (defined $query_id);
  $self->target_id($target_id) if (defined $target_id);
  
  if ((defined $fh && defined $file) ||
      !(defined $fh || defined $file)){ 
    $self->throw("Must pass in either fh or file argument");
  }
  if (defined $fh) {
    $self->{'_fh'} = $fh;
  } else {
    $self->file($file);
    open F, $self->file ||
      die "Can not open $file\n";
    $self->{'_fh'} = \*F;
  }
  $self->_initialize;
  return $self;
}

sub query_id {
  my $self = shift;
  $self->{'_query_id'} = shift if(@_);
  return $self->{'_query_id'};
}

sub target_id {
  my $self = shift;
  $self->{'_target_id'} = shift if(@_);
  return $self->{'_target_id'};
}

sub fh {
  my ($self,$value) = @_;
  if(defined $value) {
    if (ref($value) eq "GLOB") {
      $self->{'_fh'} = $value;
    } else {
      $self->throw("value for fh method should be a filehandle\n");
    }
  }
  return $self->{'_fh'};
}

sub file {
  my ($self,$value) = @_;
  if(defined $value) {
    if (-e $value) {
      $self->{'_file'} = $value;
    } else {
      $self->throw("file $value not found\n");
    }
  }
  return $self->{'_file'};
}

sub _initialize {
  my ($self) = @_;

  return undef if ($self->eof);

  my $fh = $self->fh;

  while (my $line = <$fh>) {
    last if ($line =~ /^>\S+$/);
  }
  return $self->_parsing_initialized(1);
}

sub eof {
  my ($self,$value) = @_;
  if(defined $value) {
    $self->{'_eof'} = $value;
  }
  return $self->{'_eof'};
}

sub _parsing_initialized {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'_parsing_initialized'} = $value;
  }
  return $self->{'_parsing_initialized'};
  
}

sub nextAlignment {
  my ($self) = @_;

  return undef if ($self->eof);
  
  my $fh = $self->fh;
  
  my $look_for_model_state = 0;
  my ($score, $start, $end, $hstart, $hend, $hstrand, $model_state);
  my $query_seq = "";
  my $target_seq = "";

  while (my $line = <$fh>) {
    next if ($line =~ /^$/);
    if ($line =~ /^\s+Score\s+=\s+(\d+\.\d+)\s+bits\s+\(\d+\),.*$/) {
      $score = $1;
      next;
    }
    if ($line =~ /^Query:\s+(\d+)\s+(\S+)\s+(\d+)\s+$/) {
      my ($s,$a,$e) = ($1,$2,$3);
      $start = $s unless (defined $start);
      $end = $e;
      $query_seq .= $a;
      $look_for_model_state = 1;
    }

    if ($line =~ /^\s+([A-D])\s+.*$/ && $look_for_model_state) {
      my $model_state = $1;
      $look_for_model_state = 0;
    }

    if ($line =~ /^Sbjct:\s+(-?\d+)\s+(\S+)\s+(-?\d+)\s+$/) {
      my ($s,$a,$e) = ($1,$2,$3);
      if ($s > 0) {
        $hstart = $s unless (defined $hstart);
        $hend = $e;
        $hstrand = 1 unless (defined $hstrand);
      } elsif ($s < 0) {
        $hend = abs($s) unless (defined $hend);
        $hstart = abs($e);
        $hstrand = -1 unless (defined $hstrand);
      }
      $target_seq .= $2;
    }
    
    if ($line =~ /^>\S+$/) {
      my $cigar_string = cigar_gen($query_seq, $target_seq);
      my $alignment = new Bio::EnsEMBL::DnaDnaAlignFeature(-cigar_string => $cigar_string);
      $alignment->seqname($self->query_id);
      $alignment->start($start);
      $alignment->end($end);
      $alignment->strand(1);
      $alignment->hseqname($self->target_id);
      $alignment->hstart($hstart);
      $alignment->hend($hend);
      $alignment->hstrand($hstrand);
      if (defined $self->score_in_mind) {
        $alignment->score($self->score_in_mind);
        $self->score_in_mind(undef);
      } else {
        $alignment->score($score);
      }
      return $alignment
    }

    if ($line =~ /^>--.*$/) {
      $self->score_in_mind($score);
      my $cigar_string = cigar_gen($query_seq, $target_seq);
      my $alignment = new Bio::EnsEMBL::DnaDnaAlignFeature(-cigar_string => $cigar_string);
      $alignment->seqname($self->query_id);
      $alignment->start($start);
      $alignment->end($end);
      $alignment->strand(1);
      $alignment->hseqname($self->target_id);
      $alignment->hstart($hstart);
      $alignment->hend($hend);
      $alignment->hstrand($hstrand);
      $alignment->score($score);
      return $alignment
    }
  }

  my $cigar_string = cigar_gen($query_seq, $target_seq);
  my $alignment = new Bio::EnsEMBL::DnaDnaAlignFeature(-cigar_string => $cigar_string);
  $alignment->seqname($self->query_id);
  $alignment->start($start);
  $alignment->end($end);
  $alignment->strand(1);
  $alignment->hseqname($self->target_id);
  $alignment->hstart($hstart);
  $alignment->hend($hend);
  $alignment->hstrand($hstrand);
  if (defined $self->score_in_mind) {
    $alignment->score($self->score_in_mind);
    $self->score_in_mind(undef);
  } else {
    $alignment->score($score);
  }
  
  $self->eof(1);

  return $alignment

}

sub score_in_mind {
  my $self = shift;
  $self->{'_score_in_mind'} = shift if(@_);
  return $self->{'_score_in_mind'};
}

sub cigar_gen {
  my ($q,$s) = @_;
  my @q = split //,$q;
  my @s = split //,$s;
  my $i = 0;
  my @ret = ();
  for (; $i <= $#q; $i++) {
    my $q = $q[$i];
    my $s = $s[$i];
    if($q eq "\-") {
      push @ret,"D";
      next;
    }
    if($s eq "\-") {
      push @ret,"I";
      next;
    }
    push @ret,"M";
  }
  my $c = 0;
  my $ret = "";
  for ($i=1; $i <= $#ret; $i++) {
    if ($ret[$i] eq $ret[$i-1]) {
      $c++;
      next;
    }
    if($c == 0) {
      $ret .= $ret[$i-1];
      next;
    }
    $ret .= sprintf "%d$ret[$i-1]",++$c;
    $c = 0;
  }
  if($c == 0) {
    $ret .= $ret[$i-1];
  } else {
    $ret .= sprintf "%d$ret[$i-1]",++$c;
  }
  return $ret;
}

__END__;

#
# score information :
#this subrouting repeat the previous score for entry that are <------------> (see promoterwise output).
#

sub get_all_scores {
    my $self = shift;
     if ($self->{RESULT_HEADERS})
     {
	 my @score = @{$self->{RESULT_HEADERS}};#get all lines query seq.
	     
	     my $array_ref = get_value (\@score, 3,3);
	 my @result;
	 my $tmp;
	 foreach my $seq (@$array_ref)
	 {
	     if (defined $$seq[0]){ $tmp = $$seq[0];}
	     push @result, $tmp;
	 }
	 if (@result){return \@result};
     }
}

#ithis subrouting do not return the score for entries like <--------->
#which is part of the previous alignment. (DBA output)

sub get_scores {
    my $self = shift;
     if ($self->{RESULT_HEADERS})
     {
	 my @score = @{$self->{RESULT_HEADERS}};#get all lines query seq.
	     
	     my $array_ref = get_value (\@score, 3,3);
	 my @result;
	 my $tmp;
	 foreach my $seq (@$array_ref)
	 {
	     
	     push @result, $$seq[0]; #will put seq[0] even is it's not defined.
	 }
	 if (@result){return \@result};
     }
}






















































