use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 12;}

use Bio::EnsEMBL::Pipeline::ExonPair;
use Bio::EnsEMBL::Exon;

ok(1);

ok(my $exon1 = Bio::EnsEMBL::Exon->new);

ok($exon1->start(10));

ok($exon1->end(20));

ok($exon1->strand(1));

ok($exon1->contig_id('AC00013.1'));


ok(my $exon2 = Bio::EnsEMBL::Exon->new);

ok($exon2->start(30));

ok($exon2->end(40));

ok($exon2->strand(1));

ok($exon2->contig_id('AC00013.1'));

ok(my $ep = Bio::EnsEMBL::Pipeline::ExonPair->new( -exon1 => $exon1,
						-exon2 => $exon2,
						-type  => 'silly'));


