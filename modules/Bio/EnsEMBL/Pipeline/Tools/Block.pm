package Bio::EnsEMBL::Pipeline::Tools::Block ;

use strict;
use Bio::EnsEMBL::Pipeline::Tools::Block;
use Carp;
use vars qw(@ISA);

# Object preamble - inherits from Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

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
  my ($self,$alnprog) = @_;
  my $block_initialized = 0;
  my ($newqseq,$newsseq) = ("","");
  my $block = Bio::EnsEMBL::Pipeline::Tools::Block->new();

  my $offset;
  if ($alnprog eq "blastn") {
    $offset = 1;
#    $hoffset = 1;
  } elsif ($alnprog eq "tblastx") {
    $offset = 3;
#    $hoffset = 3;
  }# elsif ($alnprog eq "tblastn") {
#    $offset = 1;
#    $hoffset = 3;
#  } elsif ($alnprog eq "blastx") {
#    $offset = 3;
#    $hoffset = 1;
#  } 

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
	  $self->qstart($self->qstart + 1 * $offset);
	} elsif ($self->qstrand == -1) {
	  $self->qend($self->qend - 1 * $offset);
	}
	$self->qseq(substr($self->qseq,$i + 1));
	$self->sseq(substr($self->sseq,$i + 1));
	$i = -1;
	next;
      }
      if (substr($self->qseq,$i,1) eq "-" && ! $block_initialized) {
	if ($self->sstrand == 1) {
	  $self->sstart($self->sstart + 1 * $offset)
	} elsif ($self->sstrand == -1) {
	  $self->send($self->send - 1 * $offset);
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
	$block->qend($self->qstart + $i  * $offset - 1);
      } elsif ($self->qstrand == -1) {
	$block->qstart($self->qend - $i  * $offset + 1);
      }
      if ($self->sstrand == 1) {
	$block->send($self->sstart + $i * $offset - 1);
      } elsif ($self->sstrand == -1) {
	$block->sstart($self->send - $i * $offset + 1);
      }
      
      if ($self->qstrand == 1) {
	$self->qstart($block->qend + 1 * $offset);
      } elsif ($self->qstrand == -1) {
	$self->qend($block->qstart - 1 * $offset);
      } 
      if ($self->sstrand == 1) {
	$self->sstart($block->send + 1 * $offset)
      } elsif ($self->sstrand == -1) {
	$self->send($block->sstart - 1 * $offset)
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
	$block->qstart($self->qstart + $i * $offset);
      } elsif ($self->qstrand == -1) {
	$block->qend($self->qend - $i * $offset);
      }
      if ($self->sstrand == 1) {
	$block->sstart($self->sstart + $i * $offset)
      } elsif ($self->sstrand == -1) {
	$block->send($self->send - $i * $offset);
      }
      $block_initialized = 1;

    } else {

      $newqseq .= substr($self->qseq,$i,1);
      $newsseq .= substr($self->sseq,$i,1);

    }
  }
}

1;
