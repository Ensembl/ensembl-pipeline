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
BEGIN { $| = 1; print "1..3\n"; 
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

$exon4 = Bio::EnsEMBL::MappedExon->new();
$exon4->start(80);
$exon4->end(90);
$exon4->strand(1);
$exon4->contig_id('AC00013.1');
$exon4->id('id-4');
$exon4->has_changed_version(1);

$exon5 = Bio::EnsEMBL::MappedExon->new();
$exon5->start(100);
$exon5->end(110);
$exon5->strand(1);
$exon5->contig_id('AC00013.1');
$exon5->id('id-5');
$exon5->has_changed_version(1);

$trans = Bio::EnsEMBL::Transcript->new();
$trans->add_Exon($exon1);
$trans->add_Exon($exon2);
$trans->add_Exon($exon3);
$trans->id('old-trans-id-1');

$gene = Bio::EnsEMBL::Gene->new();
$gene->add_Transcript($trans);
$gene->id('old-gene-id-1');

# build a merge target

$exon6 = Bio::EnsEMBL::MappedExon->new();
$exon6->start(30);
$exon6->end(40);
$exon6->strand(1);
$exon6->contig_id('AC00013.1');
$exon6->id('id-6');
$exon6->has_changed_version(1);

$exon7 = Bio::EnsEMBL::MappedExon->new();
$exon7->start(50);
$exon7->end(60);
$exon7->strand(1);
$exon7->contig_id('AC00013.1');
$exon7->id('id-7');
$exon7->has_changed_version(1);

$exon8 = Bio::EnsEMBL::MappedExon->new();
$exon8->start(80);
$exon8->end(90);
$exon8->strand(1);
$exon8->contig_id('AC00013.1');
$exon8->id('id-8');
$exon8->has_changed_version(1);

$exon9 = Bio::EnsEMBL::MappedExon->new();
$exon9->start(100);
$exon9->end(110);
$exon9->strand(1);
$exon9->contig_id('AC00013.1');
$exon9->id('id-9');
$exon9->has_changed_version(1);

$exon10 = Bio::EnsEMBL::MappedExon->new();
$exon10->start(140);
$exon10->end(190);
$exon10->strand(1);
$exon10->contig_id('AC00013.2');
$exon10->id('id-10');
$exon10->has_changed_version(1);


$transa = Bio::EnsEMBL::Transcript->new();
$transa->add_Exon($exon10);
$transa->add_Exon($exon6);
$transa->add_Exon($exon7);
$transa->id('old-trans-id-merge-1');

$genea = Bio::EnsEMBL::Gene->new();
$genea->add_Transcript($transa);
$genea->id('old-gene-id-merge-1');

$transb = Bio::EnsEMBL::Transcript->new();
$transb->add_Exon($exon8);
$transb->add_Exon($exon9);
$transb->id('old-trans-id-merge-2');

$geneb = Bio::EnsEMBL::Gene->new();
$geneb->add_Transcript($transb);
$geneb->id('old-gene-id-merge-2');

# new genes

$trans2 = Bio::EnsEMBL::Transcript->new();
$trans2->add_Exon($exon1);
$trans2->add_Exon($exon2);
$trans2->add_Exon($exon3);
$trans2->id('new-trans-id-1');

$gene2 = Bio::EnsEMBL::Gene->new();
$gene2->add_Transcript($trans);
$gene2->id('new-gene-id-1');

$trans3 = Bio::EnsEMBL::Transcript->new();
$trans3->add_Exon($exon4);
$trans3->add_Exon($exon5);
$trans3->id('new-trans-id-2');

$gene3 = Bio::EnsEMBL::Gene->new();
$gene3->add_Transcript($trans3);
$gene3->id('new-gene-id-2');

$trans4 = Bio::EnsEMBL::Transcript->new();
$trans4->add_Exon($exon10);
$trans4->add_Exon($exon6);
$trans4->add_Exon($exon7);
$trans4->add_Exon($exon8);
$trans4->id('new-trans-id-4');

$gene4 = Bio::EnsEMBL::Gene->new();
$gene4->add_Transcript($trans4);
$gene4->id('new-gene-id-4');


@tempgenes = ($gene2,$gene3,$gene4);
@oldgenes  = ($gene,$genea,$geneb);

$dummydb = DummyDB->new();


&Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Genes_to_real_Genes($dummydb,\@tempgenes,\@oldgenes);

@arr = $gene2->each_Transcript();
my $ttrans = shift(@arr);

if( $gene2->id ne 'old-gene-id-1' || $ttrans->id ne 'old-trans-id-1' ) {
  print "not ok 2\n";
} else {
  print "ok 2\n";
}


@arr = $gene3->each_Transcript();
$ttrans = shift(@arr);

if( $gene3->id !~ /generated/ || $ttrans->id !~ /generated/ ) {
  print "not ok 3\n";
} else {
  print "ok 3\n";
}


$gene2->_dump(\*STDERR);
print STDERR "\n\n";
$gene3->_dump(\*STDERR);
print STDERR "\n\n";
$gene4->_dump(\*STDERR);

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

