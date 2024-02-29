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

Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs -

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs->new(
							      -executable => $exe,
							      -db    => $db
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using getseqs by James Cuff. Sequences are fetched from a
database previously formatted with makeindex

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;

use warnings ;
use strict;
use Bio::DB::RandomAccessI;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $db) = $self->_rearrange(['EXECUTABLE', 'DB'], @args);

  if (!defined $exe) {
    $exe = 'getseqs';
  }
  $self->executable($exe);

  # expect an array of dbs
  $self->throw("Expected a reference to an array of db\n") unless ref($db) eq 'ARRAY';
  if (defined $db) {
    $self->db($db);
  }

  return $self; # success - we hope!
}

=head2 executable

  Title   : executable
  Usage   : $self->executable('/path/to/executable');
  Function: Get/set for the path to the executable being used by the module. If not set, the executable is looked for in $PATH.
  Returns : string
  Args    : string

=cut

sub executable {
  my ($self, $exe) = @_;
  if ($exe)
    {
      $self->{'_exe'} = $exe;
    }
  return $self->{'_exe'};
}

=head2 db

  Title   : db
  Usage   : $self->db('/data/blastdb/dbname');
  Function: Get/set for dbs to be searched. Checks that the database
            appropriate files are present, but nothing else.
  Returns : string
  Args    : string

=cut

sub db {

  my ($self, $dbs) = @_;
  if (!defined($self->{'_db'})) {
    $self->{'_db'} = [];
  }
  if (defined $dbs){
    if (ref($dbs) eq 'ARRAY') {
      foreach my $db(@$dbs){
	$self->throw("are you sure that $db has been formatted with makeindex?\n")
	  unless ( -e "$db.jidx");
	push (@{$self->{'_db'}},$db);
      }
    }
  }
  return (@{$self->{'_db'}});

}

=head2 get_Seq_byacc

  Title   : get_Seq_by_acc
  Usage   : $self->get_Seq_by_acc($accession);
  Function: Does the sequence retrieval via getseqs
  Returns : Bio::Seq
  Args    :

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;

  if (!defined($acc)) {
    $self->throw("No accession input");
  }

  my $seqstr;
  my $seq;
  my $getseqs = $self->executable;
  my @seqdb  = $self->db;

 SEQDB:
  while(scalar(@seqdb) && !(defined $seq)){
    my $database = pop(@seqdb);

    last SEQDB unless defined $database;
    my $cmd = "$getseqs '$acc' $database";

    open(IN,"$cmd 2>/dev/null |") or $self->throw("Error forking getseqs for accession [$acc]: getseqs");
    my $seqstr;
    while(<IN>){
      chomp;
      $seqstr .= $_;
    }

    close IN or $self->throw("Error running getseqs for $acc: $!\n");

    if(defined $seqstr && !($seqstr=~/Sorry/)) {
    chomp($seqstr);
    $seq = new Bio::Seq('-seq'               => $seqstr,
			'-accession_number'  => $acc,
			'-desc'              => "",
			'-display_id'        => $acc);

    }
  }
  $self->throw("Could not getseqs sequence for [$acc]\n") unless defined $seq;
  return $seq;
}

sub  get_Seq_by_id {
  my $self = @_;
  return undef;
}

1;
