#!/usr/local/ensembl/bin/perl -w


use strict;
use Bio::EnsEMBL::Pipeline::Config::Protein_Annotation::General qw(
								   PA_PEPTIDE_FILE
								   PA_CHUNKS_DIR
								   PA_CHUNK_SIZE
								  );



&chunk_pepfile($PA_PEPTIDE_FILE, $PA_CHUNKS_DIR, $PA_CHUNK_SIZE);



sub chunk_pepfile {
  my ($pep_file, $scratchdir, $size) = @_;
  
  #Chunk the peptide file
  open (PEPFILE, "$pep_file") or die "couldn't open $pep_file $!";
  my $count = 0;
  my $chunk = 1;
  #print STDERR "chunking peptide file\n";
  
  
  $/ = "\>";
  #print "have opened ".$pep_file."\n";
  while(<PEPFILE>){
    #print $_."\n";
    if ($_ ne "\>") {
      if ($count == 0) {
	open (CHUNK,">".$scratchdir."/chunk.$chunk") or die "couldn't opne ".$scratchdir."/chunk.$chunk";
	#print "have opened ".$scratchdir."/chunks/chunk.$chunk\n";
      }
      
      $_ =~ s/\>$//;  
      
      print CHUNK ">$_";
      $count++;
      if ($count == $size) {
	$count = 0;
	$chunk++;
      }
    }
  }
  $/ = "\n";
}
