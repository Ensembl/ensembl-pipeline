package Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::Seq;

@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;
  
  my ($db) = $self->_rearrange([qw(DB)],@args);

  $self->db($db) if $db;

  return $self
}


sub _fetch {
  my ($self, $id) = @_;

  my $command;

  if (ref $id eq 'ARRAY') {
    $command = $self->db->seqfetch_command . " @$id";
  } else {
    $command = $self->db->seqfetch_command . " " . $id;
  }

  open(CMD, "$command |") or die "Can't execute fetch command";

  my %seqs;
  my $id_line;

  while (<CMD>){

    if (/^>/){
      $id_line = $_;
      $id_line =~ s/^>//;

      if ($self->db->index_type =~ /wu/){
	$id_line =~ s/([\w\_\.]+).*\n/$1/;
      }
      if ($self->db->index_type eq 'ncbi'){
	$id_line =~ s/[^\|]+\|([\w\_\.]+).*\n/$1/;
      }

      next
    }

    $seqs{$id_line} .= $_;
  }

  close CMD;

  my @bioseqs;

  foreach my $seq_id (keys %seqs) {

    $seqs{$seq_id} =~ s/\n//g;
    
    my $bioseq = Bio::Seq->new(-display_id => $seq_id,
			       -seq        => $seqs{$seq_id}
			      );
    
    push (@bioseqs, $bioseq)
  }

  return \@bioseqs;
}

### Some aliae

sub fetch {
  my ($self, $id) = @_;

  my $seqs = $self->_fetch($id);

  return shift @$seqs
}

sub get_Seq_by_acc {
  my ($self, $id) = @_;

  return $self->fetch($id)
}

sub batch_fetch {
  my ($self, $ids) = @_;

  return $self->_fetch($ids);
}


sub db {
  my $self = shift;

  if (@_) {
    
    my $value = shift;
    
    unless ($value->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB")){
      $self->throw("Blast database object is not a "
		   ."Bio::EnsEMBL::Pipeline::Runnable::BlastDB.\n")
    }
    
    $self->{_db} = $value;
  }
  
  return $self->{_db}
}

