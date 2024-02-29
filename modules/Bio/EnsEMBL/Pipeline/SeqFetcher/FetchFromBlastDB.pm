=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB -

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB;

use warnings ;
use strict;
use Bio::Seq;

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
    for (my $i = 0; $i < scalar @$id; $i++){
      $id->[$i] =~ s/(\S+).*/$1/;
    }

    $command = $self->db->seqfetch_command . " @$id";
  } else {
    $id =~ s/(\S+).*/$1/;
    $command = $self->db->seqfetch_command . " " . $id;
  }

  open(CMD, "$command |") or die "Can't execute fetch command";

  my %seqs;
  my %descs;
  my $id_line;
  my $desc;

  while (<CMD>){

    if (/^>/){
      $id_line = $_;
      $id_line =~ s/^>//;

      if ($self->db->index_type =~ /wu/){
	$id_line =~ s/([\w\_\.]+)\s*(.*)\n/$1/;
	$desc = $2;
      }
      if ($self->db->index_type eq 'ncbi'){
	$id_line =~ s/[^\|]+\|([\w\_\.]+)\s*(.*)\n/$1/;
	$desc = $2;
      }

      $desc =~ s/\t/ /g;

      $descs{$id_line} = $desc;

      next
    }

    $seqs{$id_line} .= $_;
  }

  close CMD;

  my @bioseqs;

  foreach my $seq_id (keys %seqs) {

    $seqs{$seq_id} =~ s/\n//g;

    my $bioseq = Bio::Seq->new(-display_id => $seq_id,
			       -seq        => $seqs{$seq_id},
			       -desc       => $descs{$seq_id}
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
