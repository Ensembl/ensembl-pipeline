use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 14;};

use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

my $pwd = `pwd`; chomp($pwd);

ok(my $dnaseqio = new Bio::SeqIO('-file' => "$pwd/t/data/AP000074.fa", 
				 '-format' => 'fasta'));

ok(my $dnaseq   = $dnaseqio->next_seq->seq);

ok(my $clone    =  Bio::PrimarySeq->new(  '-seq'         => $dnaseq,
					  '-id'          => 'AP000074',
					  '-accession'   => 'AP000074',
					  '-moltype'     => 'protein'));

ok(my $pepseqio = new Bio::SeqIO('-file' => "$pwd/t/data/AP000074.pep",  
				 '-format' => 'fasta'));

ok(my $pepseq   = $pepseqio->next_seq->seq);

ok(my $peptide  =  Bio::PrimarySeq->new(  '-seq'         => $pepseq,
					  '-id'          => 'AP000074',
					  '-accession'   => 'AP000074',
					  '-moltype'     => 'protein'));

# First lets test the dna-dna blast
ok(my $blastn = Bio::EnsEMBL::Pipeline::Runnable::Blast->new
    ('-query'    => $clone,
     '-program'  => 'wublastn',
     '-database' => "$pwd/t/data/AI053588.fa",
     '-threshold' => 1e-6,
     ));

if(!defined $blastn->get_regex("$pwd/t/data/AI053588.fa")){
  print STDERR "No regex defined for $pwd/t/data/AI053588.fa - adding one\n";
  #>AI053588 AI053588 qi68h04.x1 NCI_CGAP_Ov26 ...
  $blastn->add_regex("$pwd/t/data/AI053588.fa", '^(\w+)\s+');
}

ok($blastn->run);

foreach my $out ($blastn->output) {
    print $out->gffstring . "\n";
}

# Now the dna-pep blast
ok(my $blastx = Bio::EnsEMBL::Pipeline::Runnable::Blast->new 
    ('-query'    => $clone,
     '-program'  => 'wublastx',
     '-database' => "$pwd/t/data/AP000074.pep",
     '-threshold' => 1e-6,
     ));

if(!defined $blastn->get_regex("$pwd/t/data/AP000074.pep")){
  print STDERR "No regex defined for $pwd/t/data/AP000074.pep - adding one\n";
  # >AP000074|GENSCAN_predicted_peptide_1|956_aa
  $blastn->add_regex("$pwd/t/data/AP000074.pep", '^(\w+)\s*');
}

ok($blastx->run);

foreach my $out ($blastx->output) {
    print $out->gffstring . "\n";
}

# Now the pep-dna blast

ok(my $tblastn = Bio::EnsEMBL::Pipeline::Runnable::Blast->new 
    ('-query'    => $peptide,
     '-program'  => 'wutblastn',
     '-database' => "$pwd/t/data/AI053588.fa",
     '-threshold' => 1e-6,
     ));

ok($tblastn->run);

foreach my $out ($tblastn->output) {
    print $out->gffstring . "\n";
}

ok(1);

