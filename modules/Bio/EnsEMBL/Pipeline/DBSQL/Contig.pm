#
# BioPerl module for Contig
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL:Contig - Handle onto a database stored contig

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::DBSQL::Contig;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

BEGIN {
    print STDERR "Warning - file still using DBSQL::Contig. Should use RawContig!\n";
};

use Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::RawContig;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::RawContigI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my ($id,$dbobj) = $self->_rearrange([qw(ID
					  DBOBJ
					  )],@args);

  $id    || $self->throw("Cannot make contig db object without an id");
  $dbobj || $self->throw("Cannot make contig db object without db object");

  $dbobj->isa('Bio::EnsEMBL::Pipeline::DBSQL::Obj') || $self->throw("Cannot make contig db object with a [$dbobj] object");

  $self->id($id);
  $self->_dbobj($dbobj);

  $self->{_exon_pairs} = [];
  $self->{_genes}      = [];
  return $make; # success - we hope!
}


sub id {
    my ($self,$value) = @_;
    
    if( defined $value) {
      $self->{'_id'} = $value;
    }
    return $self->{'_id'};

}

sub _dbobj {
    my ($self,$value) = @_;
    
    if ( defined $value && $value->isa("Bio::EnsEMBL::Pipeline::DB::ObjI")) {
	$self->{'_dbobj'} = $value;
    } else {
	$self->throw("[$value] is not a Bio::EnsEMBL::Pipeline::DB::ObjI");
    }

    return $self->{'_dbobj'};

}

sub add_SimilarityFeature {
    my ($self,$arg) = @_;

    if (defined($arg) && $arg->isa("Bio::EnsEMBL::SeqFeatureI")) {
	if (!(defined($self->{_similarity_features}))) {
	    $self->{_similarity_features} = [];
	}
	push(@{$self->{_similarity_features}},$arg);
    }
}
sub get_all_SimilarityFeatures {
    my ($self) = @_;


    if (defined($self->{_similarity_features})) {
	return @{$self->{_similarity_features}};
    } 
}
    
sub make_ExonPairs {
    my ($self,@exons) = @_;

    my $gap = 5;

    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
    } else {
	@exons = sort {$b->start <=> $a->start} @exons;
    }

    my %pairhash;

    for (my $i = 0; $i < scalar(@exons)-1; $i++) {

	my %idhash;
	my $exon1 = $exons[$i];

	    for (my $j = $i+1 ; $j < scalar(@exons); $j++) {

		my $exon2 = $exons[$j];

		my %doneidhash;

		F1: foreach my $f1 ($exon1->each_Supporting_Feature) {
		    
		    F2: foreach my $f2 ($exon2->each_Supporting_Feature) {

			next F1 if (!($f1->isa("Bio::EnsEMBL::FeaturePair")));
			next F2 if (!($f2->isa("Bio::EnsEMBL::FeaturePair")));

			if ($f1->hseqname eq $f2->hseqname &&
			    $f1->strand   == $f2->strand   &&
			    !(defined($idhash{$f1->hseqname})) &&
			    !(defined($pairhash{$exon1}{$exon2}))) {

			    my $ispair = 0;

			    if ($f1->strand == 1) {
				if (abs($f2->hstart - $f1->hend) < $gap) {
				    if (!(defined($doneidhash{$f1->hseqname}))) {
					$ispair = 1;
				    }
				}
			    } elsif ($f1->strand == -1) {
				if (abs($f1->hend - $f2->hstart) < $gap) {
				    if (!(defined($doneidhash{$f1->hseqname}))) {
					$ispair = 1;
				    }
				}
			    }

			    if ($ispair == 1) {
				eval {
				    print(STDERR "Making new pair " . $exon1->id . "\t" .  $exon2->id . "\n");
				    
				    my $pair = $self->makePair($exon1,$exon2);
				    
				    $idhash    {$f1->hseqname} = 1;
				    $doneidhash{$f1->hseqname} = 1;

				    if ($pair->coverage > 1) {
					$pairhash{$exon1}{$exon2}  = 1;
				    }

				};
				if ($@) {
				    warn("Error making ExonPair from [" . $exon1->id . "][" .$exon2->id ."] $@");
				}
			    }
			    
			}
		    }
		}
	    }
    }
}


sub makePair {
    my ($self,$exon1,$exon2) = @_;

    my $tmppair = new Bio::EnsEMBL::Pipeline::ExonPair(-exon1 => $exon1,
						       -exon2 => $exon2,
						       -type  => "ABUTTING",
						       );
    my $found = 0;

    foreach my $p ($self->get_all_ExonPairs) {
	if ($p->compare($tmppair) == 1) {
	    $p->add_coverage;
	    $tmppair = $p;
	    $found = 1;
	}
    }

    if ($found == 0) {
	$self->add_ExonPair($tmppair);
    }

    return $tmppair;
}

sub get_all_ExonPairs {
    my ($self) = @_;

    return @{$self->{_exon_pairs}};
}

sub add_ExonPair {
    my ($self,$arg) = @_;

    if (defined($arg) && $arg->isa("Bio::EnsEMBL::Pipeline::ExonPair")) {
	push(@{$self->{_exon_pairs}},$arg);
    } else {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::ExonPair");
    }
}

sub get_Heads {
    my ($self) = @_;
    
    my $sth = $self->_dbobj->prepare("select e1.exon1_id " .
				     "from exon_pair as e1 " .
				     "left join exon_pair as e2 on e1.exon1_id = e2.exon2_id " .
				     "where e2.exon1_id is null");
    
    my $res = $sth->execute;
    
    my $rowhash = $sth->fetchrow_hashref;
    
}
	
sub get_Tails {
    my ($self) = @_;
    
    my $sth = $self->_dbobj->prepare("select e1.exon2_id " .
				     "from exon_pair as e1 " .
				     "left join exon_pair as e2 on e1.exon2_id = e2.exon1_id " .
				     "where e2.exon1_id is null");
    
    my $res     = $sth->execute;
    my $rowhash = $sth->fetchrow_hashref;
    
}
		


sub add_Gene {
    my ($self,$gene) = @_;

    if (!(defined($self->{_genes}))) {
	$self->{_genes} = [];
    }

    push(@{$self->{_genes}},$gene);

}

sub get_all_Genes {
    my ($self) = @_;

    return @{$self->{_genes}};

}
