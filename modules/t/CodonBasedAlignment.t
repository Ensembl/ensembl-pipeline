
use strict;
use lib 't';
use Test;
use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;

BEGIN { $| = 1; plan test => 5;}

ok(1);

my $cba = 
  Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment->new(
         -genetic_code => 1);

ok($cba->isa("Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment"));

my $aligned_seqs;


{
  $cba->filename('t/data/relaxins.fa');
  $cba->sequences_from_file;
  ok($cba->sequences);

  ok($aligned_seqs = $cba->run_alignment);
}

ok(@$aligned_seqs);
