## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..2\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}


use Bio::EnsEMBL::MappedExon;
use Bio::EnsEMBL::Pipeline::GeneComp;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::DB::ObjI;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

$exon1 = Bio::EnsEMBL::MappedExon->new();
$exon1->start(10);
$exon1->end(20);
$exon1->strand(1);
$exon1->contig_id('AC00013.1');
$exon1->id('id-1');
$exon1->has_changed_version(1);

$exon2 = Bio::EnsEMBL::MappedExon->new();
$exon2->start(30);
$exon2->end(40);
$exon2->strand(1);
$exon2->contig_id('AC00013.1');
$exon2->id('id-2');
$exon2->has_changed_version(1);


$exon3 = Bio::EnsEMBL::MappedExon->new();
$exon3->start(50);
$exon3->end(60);
$exon3->strand(1);
$exon3->contig_id('AC00013.1');
$exon3->id('id-3');
$exon3->has_changed_version(1);


$trans = Bio::EnsEMBL::Transcript->new();
$trans->add_Exon($exon1);
$trans->add_Exon($exon2);
$trans->add_Exon($exon3);
$trans->id('old-trans-id-1');

$gene = Bio::EnsEMBL::Gene->new();
$gene->add_Transcript($trans);
$gene->id('old-gene-id-1');



$trans2 = Bio::EnsEMBL::Transcript->new();
$trans2->add_Exon($exon1);
$trans2->add_Exon($exon2);
$trans2->add_Exon($exon3);
$trans2->id('new-trans-id-1');

$gene2 = Bio::EnsEMBL::Gene->new();
$gene2->add_Transcript($trans);
$gene2->id('new-gene-id-1');

@tempgenes = ($gene2);
@oldgenes  = ($gene);

$dummydb = DummyDB->new();


&Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Genes_to_real_Genes($dummydb,\@tempgenes,\@oldgenes);
@arr = $gene2->each_Transcript();

my $ttrans = shift(@arr);

if( $gene2->id ne 'old-gene-id-1' || $ttrans->id ne 'old-trans-id-1' ) {
  print "not ok 2\n";
} else {
  print "ok 2\n";
}


$gene2->_dump(\*STDERR);



package DummyDB;

@DummyDB::ISA = ('Bio::EnsEMBL::Pipeline::DB::ObjI');

$id = 1;

sub new {
    my $class = shift;
    my $self = {};

    bless $self,$class;
    return $self;
}

sub get_new_GeneID {
    return "generated-gene-id" . $id++ ;
}

sub get_new_ExonID {
    return "generated-exon-id" . $id++ ;
}

sub get_new_TranscriptID {
    return "generated-transcript-id" . $id++ ;
}

