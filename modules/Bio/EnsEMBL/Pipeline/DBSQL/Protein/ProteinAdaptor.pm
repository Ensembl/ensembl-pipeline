# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

  Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinAdaptor - Object representing an adaptor
                                                           for proteins to an EnsEMBL DB

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::Protein::DBAdaptor;
  use Bio::EnsEMBL::DBSQL::Protein::Protein_Adaptor;

  my $db = Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor->new ( -host   => 'ics1e',
                                                                    -dbname => 'test',
                                                                    -user   => 'marc',
                                                                  );               

  my $proteinAdaptor = $db->get_ProteinAdaptor;
  my $protein = $proteinAdaptor->fetch_Protein_by_dbid ($id);

=head1 DESCRIPTION

  This Object inherits from BaseAdaptor, following the new adaptor rules.
  The main method is fetch_Protein_by_dbid, which returns a complete protein object.


=head1 CONTACT

  Marc Sohrmann (ms2@sanger.ac.uk)

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::PrimarySeq;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_all_Peptide

 Title    : fetch_all_Peptide
 Usage    : $self->fetch_all_Peptide
 Function : get all the Peptides from the peptide table, building PrimarySeq objects
 Example  : 
 Returns  : array of Bio::PrimarySeq objects, or empty array if no Peptides are found
 Args     : 
 Throws   :

=cut

sub fetch_all_Peptide {
    my ($self) = @_;
    my $sth = $self->prepare ( q{ SELECT proteinId, peptide
                                    FROM peptide
                                } );
    $sth->execute;
    my $rows = $sth->fetchall_arrayref;   # returns undef if there is no data
                                          # in the statement handle
    my @array = ();
    foreach my $entry (@$rows) {
        unless (defined $entry->[0] && defined $entry->[1]) {
            $self->throw ("problems getting the protein id and seq from the DB");
	}
        my $pep = Bio::PrimarySeq->new ( -id      => $entry->[0],
                                         -seq     => $entry->[1],
                                         -moltype => 'protein',
             		               );
        push (@array, $pep);
    }
    $sth->finish;
    return @array;
}

=head2 fetch_Peptide_by_dbid

 Title    : fetch_Peptide_by_dbid
 Usage    : $self->fetch_Peptide_by_dbid
 Function : get a Peptide from the peptide table, building a PrimarySeq object
 Example  : 
 Returns  : a Bio::PrimarySeq object
 Args     : id
 Throws   :

=cut

sub fetch_Peptide_by_dbid {
    my ($self, $id) = @_;
    my $sth = $self->prepare ( q{ SELECT proteinId,  peptide
                                    FROM peptide
                                   WHERE proteinId = ?
                                } );
    $sth->execute ($id);
    my @row = $sth->fetchrow_array;   # returns undef if there is no data
                                      # in the statement handle
    unless (defined $row[0] && defined $row[1]) {
        $self->throw ("problems getting the protein id and seq from the DB");
    }
    my $pep = Bio::PrimarySeq->new ( -id      => $row[0],
                                     -seq     => $row[1],
                                     -moltype => 'protein',
         		           );
    $sth->finish;
    return $pep;
}


=head2 fetch_Protein_by_dbid

 Title    : fetch_Protein_by_dbid
 Usage    : $self->fetch_Protein_by_dbid
 Function : get a Protein from the peptide table, building a PrimarySeq object
 Example  : 
 Returns  : a Bio::PrimarySeq object
 Args     : id
 Throws   :

=cut

sub fetch_Protein_by_dbid {
   my ($self, $id) = @_;
   $self->fetch_Peptide_by_dbid ($id);
} 

=head2 fetch_Protein_ligth

 Title    : fetch_Protein_ligth
 Usage    : $self->fetch_Protein_ligth
 Function : an alias for a function called by Emmanual in Hmmpfam.pm,
            get a Protein from the peptide table, building a PrimarySeq object
 Example  : 
 Returns  : a Bio::PrimarySeq object
 Args     : id
 Throws   :

=cut

sub fetch_Protein_ligth {
   my ($self, $id) = @_;
   $self->fetch_Peptide_by_dbid ($id);
} 


=head2 fetch_PepSeq_by_dbid

 Title    : fetch_PepSeq_by_dbid
 Usage    : $self->fetch_PepSeq_by_dbid
 Function : get a Peptide sequence from the peptide table
 Example  : 
 Returns  : a sequence string
 Args     : id
 Throws   :

=cut

sub fetch_PepSeq_by_dbid {
    my ($self, $id) = @_;
    my $sth = $self->prepare ( q{ SELECT peptide
                                    FROM peptide
                                   WHERE proteinId = ?
                                } );
    $sth->execute ($id);
    my @row = $sth->fetchrow_array;   # returns undef if there is no data
                                      # in the statement handle
    unless (defined $row[0]) {
        $self->throw ("problems getting the protein seq from the DB");
    }
    return $row[0];
}


=head2 update_swallId

 Title    : update_swallId
 Usage    : $self->update_swallId
 Function : updates the swallId in the protein table
 Example  : 
 Returns  : nothing
 Args     : id, swallId
 Throws   :

=cut

sub update_swallId {
    my ($self, $proteinId, $swallId) = @_;
    my $sth = $self->prepare ( q{ UPDATE protein
                                     SET swallId = ?
                                   WHERE proteinId = ?
                                } );
    $sth->execute ($swallId, $proteinId);
    $sth->finish;
    return;
}


=head2 submit

 Title    : submit
 Usage    : $self->submit
 Function : populate the protein, peptide and InputIdAnalysis tables
 Example  : 
 Returns  : nothing
 Args     : Bio::Seq or Bio::PrimarySeq
 Throws   :

=cut

sub submit {
    my ($self, $prot, $org) = @_;
    unless ($prot->isa ('Bio::Seq') || $prot->isa ('Bio::PrimarySeq')) {
        $self->throw ("only Bio::Seq and Bio::PrimarySeq objects can be submitted");
    }
    my $sth1 = $self->prepare ( q{ INSERT INTO protein
                                               (proteinId, created, swallId, length, organismId)
                                        VALUES (?, now(), ?, ?, ?)
                                 } );

    my $sth2 = $self->prepare ( q{ INSERT INTO peptide
                                               (proteinId, peptide)
                                        VALUES (?, ?)
                                 } );

    my $sth3 = $self-> prepare ( q{ INSERT INTO InputIdAnalysis
                                                (inputID, class, analysisId, created)
	                                 VALUES (?, ?, ?, now())
                                  } );

    $sth1->execute ($prot->id, (($prot->accession_number) && ($prot->accession_number ne "unknown")) ? $prot->accession_number : undef, $prot->length, $org);
    $sth2->execute ($prot->id, $prot->seq);
    $sth3->execute ($prot->id, (defined $prot->moltype) ? $prot->moltype : 'protein', 1); 

    $sth1->finish;
    $sth2->finish;
    $sth3->finish;
    return;
}



