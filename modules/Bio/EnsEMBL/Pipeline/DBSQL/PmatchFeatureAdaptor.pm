=head1 NAME

=head1 SYNOPSIS

  
=head1 DESCRIPTION


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::PmatchFeatureAdaptor;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::DBSQL::BaseAdaptor
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


use Bio::EnsEMBL::Pipeline::PmatchFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


sub write_PmatchFeatures {
  my ($self,@features) = @_;


  my %protein;

  foreach my $f (@features) {
    if (!defined($protein{$f->protein_id})) {
      $protein{$f->protein_id} =  $self->write_protein($f->protein_id,$f->cdna_id);
    }
    
    my $protein_internal_id = $protein{$f->protein_id};

    if (!defined($protein_internal_id)) {
      $self->throw("No internal id found for " . $f->protein_id . "\n");
    }

    my $query = "insert into pmatch_feature values(null," . 
        $protein_internal_id. ",'" . 
	$f->chr_name . "'," .  
	$f->start   . ","  . 
	$f->end     . ","  .
	$f->coverage . ")";

    my $sth = $self->prepare($query);

    my $res = $sth->execute;

  }
}

sub write_protein {
  my ($self,$protein_id,$cdna_id) = @_;

  my $tmpcdna = $self->get_cdna_id($protein_id);

  if (defined ($tmpcdna)) {
    if ($tmpcdna ne $cdna_id) {
      $self->throw("ERROR: Protein $protein_id already exists with different cdna $tmpcdna : $cdna_id\n");
    }

    return $self->get_protein_internal_id($protein_id);

  } else {

    my $query = "insert into protein values(null,'$protein_id','$cdna_id')";
    my $sth = $self->prepare($query);
    my $res = $sth->execute;

    $sth = $self->prepare("select LAST_INSERT_ID()");
    $res = $sth->execute;

    my ($id) = $sth->fetchrow  or $self->throw("Failed to get last insert id");
    
    return $id;
  }
}

sub get_protein_internal_id {
  my ($self,$protein_id) = @_;

  my $query = "select protein_internal_id from protein where protein_id = '$protein_id'";
  my $sth = $self->prepare($query);
  my $res = $sth->execute;

  if ($sth->rows > 0) {
    my $row = $sth->fetchrow_hashref;
    my $id = $row->{'protein_internal_id'};
    return $id;
  }
}

sub delete_protein {
  my ($self,$protein_id) = @_;

  my $internal_id = $self->get_protein_internal_id($protein_id);

  if ($internal_id eq "") {
    $self->throw("Protein $protein_id does not exist in the database");
  }

  my $query = "delete from protein where protein_id = '$protein_id'";

  my $sth = $self->prepare($query);
  my $res = $sth->execute;

  $query  = "delete from pmatch_feature where protein_internal_id = $internal_id";

  $sth    = $self->prepare($query);
  $res    = $sth->execute;

}
    

sub get_PmatchFeatures_by_protein_id {
  my ($self,$prot_id) = @_;

  my $query = "select * from pmatch_feature,protein where protein.protein_internal_id = pmatch_feature.protein_internal_id and protein.protein_id = '$prot_id'";

  my $sth = $self->prepare($query);

  my $res  = $sth->execute;

  my @pmatch;
  my %cdnas;

  while (my $row = $sth->fetchrow_hashref) {
    my $prot_internal_id = $row->{protein_internal_id};
    my $chr_name         = $row->{chr_name};
    my $start            = $row->{start};
    my $end              = $row->{end};
    my $coverage         = $row->{coverage};

    if (!defined($cdnas{$prot_internal_id})) {

      $cdnas{$prot_internal_id}    = $self->get_cdna_id($prot_internal_id);

    }

    my $feature = new Bio::EnsEMBL::Pipeline::PmatchFeature(-protein_id  => $prot_id,
							    -start    => $start,
							    -end      => $end,
							    -chr_name => $chr_name,
							    -cdna_id  => $cdnas{$prot_internal_id},
							    -coverage => $coverage,
							    );

    push(@pmatch,$feature);

  }
  return @pmatch;
}

sub get_cdna_id {
  my ($self,$protein_id) = @_;
  
  my $query = "select cdna_id from protein where protein_id = '$protein_id'";

  my $sth = $self->prepare($query);

  my $res  = $sth->execute;
  my $cdna_id;

  while (my $row = $sth->fetchrow_hashref) {
    $cdna_id = $row->{cdna_id};
  }
  return $cdna_id;
}


sub get_protein_id {
  my ($self,$prot_internal_id) = @_;

  my $query = "select protein_id from protein where protein_internal_id = $prot_internal_id";

  my $sth = $self->prepare($query);

  my $res  = $sth->execute;
  my $prot_id;

  while (my $row = $sth->fetchrow_hashref) {
     $prot_id = $row->{protein_id};
  }
  return $prot_id;
}


sub get_cdna_id_by_internal_protein_id {
  my ($self,$prot_internal_id) = @_;

  my $query = "select cdna_id from protein where protein_internal_id = $prot_internal_id";

  my $sth = $self->prepare($query);

  my $res  = $sth->execute;
  my $cdna_id;

  while (my $row = $sth->fetchrow_hashref) {
    $cdna_id = $row->{cdna_id};
  }
  return $cdna_id;
}


sub exists_protein {
  my ($self,$protein_id) = @_;

  my $query = "select * from protein where protein_id = '$protein_id'";
  
  my $sth = $self->prepare($query);
  my $res = $sth->execute;

  if ($sth->rows > 0) {
    return 1;
  } else {
    return 0;
  }
}

sub create_sql {
  my($self,$dbname) = @_;

  my $create = "create database $dbname";

  my $protein = "create table protein ( protein_internal_id int(10) unsigned NOT NULL auto_increment,
					  protein_id          varchar(40) NOT NULL,
					  cdna_id             varchar(40) NOT NULL,
					  
					  PRIMARY KEY (protein_internal_id),
					  UNIQUE  protein_id(protein_id)
					);";
  
  
  my $feature = "CREATE TABLE pmatch_feature (feature_internal_id int(10) unsigned NOT NULL auto_increment,
						protein_internal_id int(10) unsigned NOT NULL,
						chr_name            varchar(40) NOT NULL, 
						start               int(10) unsigned NOT NULL,
						end                 int(10) unsigned NOT NULL,
				 		coverage            double(16,4) NOT NULL,
						 
						PRIMARY KEY(feature_internal_id),
                                                UNIQUE(protein_internal_id,chr_name,start,end),
						KEY(protein_internal_id));   );";

  return $create . "\n" . $protein . "\n" . $feature;
}

1;      
    

