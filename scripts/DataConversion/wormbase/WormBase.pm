package WormBase;
require Exporter;


our @ISA = qw(Exporter);
our @EXPORT = qw(get_seq_ids get_sequences agp_parse);


sub get_seq_ids{
  my ($fh) = @_;

  my @seq_ids;
  
  while(<$fh>){
   chomp;
   #I	47490	107680	3	F	AC024796.1	1	60191	+
   #print;
   #print "\n";
   my ($status, $contig) =
    (split)[4, 5];
   if($status eq 'N'){
     next;
   }
   if($contig eq '.'){
     next;
   }
   push(@seq_ids, $contig)
  }
  return \@seq_ids;
}


sub get_sequences{
  my ($seq_ids, $seqfetcher) = @_;
  my %seqs;
  foreach my $id(@$seq_ids){
    my $seq = $seqfetcher->get_Seq_by_acc($id);
    $seqs{$id} = $seq;
  }
  return(\%seqs);
}


sub agp_parse{
  my ($fh, $chr_id, $agp_type) = @_;
  my $chr_hash = {};
  while(<$fh>){
    chomp;
    #I	47490	107680	3	F	AC024796.1	1	60191	+
    #print;
    #print "\n";
    my ($chr, $chr_start, $chr_end, $gap,  $contig, $raw_start, $raw_end, $raw_ori) =
      (split)[0, 1, 2, 4, 5, 6, 7, 8];
    if($gap eq 'N'){
      next;
    }
    if($contig eq '.'){
      next;
    }
    if ($raw_ori eq '+') {
      $raw_ori = 1;
    }
    elsif ($raw_ori eq '-') {
      $raw_ori = -1;
    }
    else {
      $raw_ori = 1;
      #print "$chr Contig  $contig  $chr_start \n";
      #print "Warning assumed default orientation for $contig\n";
    }
    
    if($chr_hash->{$contig}){
      die "contig ".$contig." has been found twice something odd is going on\n";
    }else{
      $chr_hash->{$contig} = { 'chromosome_id' => $chr_id,
			       'chr_start' => $chr_start,
			       'chr_end' => $chr_end,
			       'superctg_name' => $chr,
			       'superctg_start' => $chr_start,
			       'superctg_end' => $chr_end,
			       'superctg_ori' => 1,
			       'contig_start' => $raw_start,
			       'contig_end' => $raw_end,
			       'contig_ori' => $raw_ori,
			       'type' => $agp_type,
			     }

    }
  }
  return $chr_hash;
}
