#! /usr/local/bin/perl
#
#all things about bsub jobs
#


my $start = 250;
my $length = 250;
my $end = 56250;
my $count = 161;
my $queue = "ecslong";
for ($start=250; $start<=$end; $start+=$length) {
  print "job_num is ", ++$count, "start is $start\n";
  if ($count>161 and $count<=225) {
    $queue = "ecslong";
  }
  elsif ($count>225) {
    exit;
    $queue = "acarilong";
  }
  #`./snp_ssaha_match.pl $start`;
  #`cat /nfs/acari/yuan/ensembl/test/SNP_SSAHA_FROM_NO_CLONE_seq_$start | /usr/local/badger/bin/ssahaClient tcs1a 50000 20 0 0 0 DNA 10000 3 size>/nfs/acari/yuan/ensembl/test/NO_CLONE_seq_$start`;
  `bsub -q $queue -o /work4/yuan/OUT1/ssaha_out_$start  -e /work4/yuan/ERROR1/ssaha_error_$start /nfs/acari/yuan/ensembl/test/snp_ssaha_match.pl $start`;
}

#for (my $start=944286; $start <= 944402; $start++) {
# `bkill $start`;
#}























































































