
use strict;
use Bio::EnsEMBL::Pipeline::Tools::Bl2seq;
use Bio::EnsEMBL::Pipeline::Runnable::Bl2seq;
use Bio::SeqIO;
use lib 't';
use Test;

BEGIN { $| = 1; plan test => 5;}

ok(1);

ok(open(RESULT, '<t/data/bl2seq.out'));

my $tools_bl2seq 
   = Bio::EnsEMBL::Pipeline::Tools::Bl2seq->new(
                  '-fh'        => \*RESULT,
		  '-alntype'   => 'blastn',
		  '-min_score' => 0,
		  '-qname'     => 'Hsa_RLN1',
		  '-sname'     => 'Hsa_RLN2');

ok($tools_bl2seq->isa("Bio::EnsEMBL::Pipeline::Tools::Bl2seq"));

ok(my $dna_align_feature = $tools_bl2seq->nextHSP);

ok($dna_align_feature->isa("Bio::EnsEMBL::DnaDnaAlignFeature"));

close RESULT;
