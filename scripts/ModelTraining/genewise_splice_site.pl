#!/usr/local/ensembl/bin/perl -w
#This script will produce 2 files (5' and 3' splice site regions).
#These 2 files can then be used to create a genewise gene.stat model
#Contact ensembl-dev@ebi.ac.uk

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dbpass    = undef;
my $dnadbhost;
my $dnadbname;
my $up;
my $down;

my $genetype; # default genetype


GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$dbhost,
	   'dnadbhost:s' => \$dnadbhost,
	   'dnadbname:s' => \$dnadbname,
	   'genetype:s'  => \$genetype,
	   'up:s'=> \$up,
	   'down:s' =>\$down, 
	  );

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost ){
  print STDERR "Usage: $0 -dbname -dbhost -dnadbhost -dnadbname -up -down\n";
  exit(0);
}



my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dbuser,
					       '-dbname' => $dnadbname,
					       );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-dnadb'  => $dnadb,
					   );

print STDERR "connected to $dbname : $dbhost\n";

print STDERR "checking genes of type $genetype\n";

open (SPLICEUP,">$up") || die;
open (SPLICEDOWN,">$down") || die;

my @gene_ids;

if ($genetype) {
    my $query = "select gene_id from gene where type = '$genetype'";
    my $sth = $db->prepare($query);
    $sth->execute();
    while (my $id = $sth->fetchrow) {
	push (@gene_ids,$id);
    }
}

else {
    @gene_ids      = @{$db->get_GeneAdaptor->list_geneIds};
}

my $slice_adaptor = $db->get_SliceAdaptor;

 GENE: foreach my $gene_id ( @gene_ids){
     
     my $gene;
    
     
     $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id,1);
     
     if ($genetype){
	 next GENE unless ( $gene->type eq $genetype );
     }
     
	
    my @transcripts = @{$gene->get_all_Transcripts};
     
     foreach my $transcript(@transcripts) { 
	 
	 $transcript->sort;
	 
	 my $strand = $transcript->start_Exon->strand;
	 my @exons  = @{$transcript->get_all_Exons};
	 
	 my $introns  = scalar(@exons) - 1 ; 
	 
	 my $correct        = 0;
	 
	 # all exons in the transcripts are in the same seqname coordinate system:
	 my $slice = $transcript->start_Exon->contig;
	 
	 if ($strand == 1 ){
	     
	   INTRON:
	     for (my $i=0; $i<$#exons; $i++ ){
		 
		 my $upstream_exon   = $exons[$i];
		 
		 my $downstream_exon = $exons[$i+1];
		 
		 my $upstream_site;
		 my $downstream_site;
		 
		 my $upstream_site1;
		 my $downstream_site1;
		 
		
		 eval{
		     $upstream_site1 = 
			 $slice->subseq( ($upstream_exon->end     - 3), ($upstream_exon->end     + 10 ) );
		     
		     $downstream_site1 = 
			 $slice->subseq( ($downstream_exon->start - 10), ($downstream_exon->start + 3 ) );
		     
		 };
		 unless ( $upstream_site1 && $downstream_site1 ){
		     print STDERR "problems retrieving sequence for splice sites\n$@";
		     next INTRON;
		 }
		 
		 eval{
		     $upstream_site = 
			 $slice->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
		     $downstream_site = 
			 $slice->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );
		 };
		 unless ( $upstream_site && $downstream_site ){
		     print STDERR "problems retrieving sequence for splice sites\n$@";
		     next INTRON;
		 }
		 
		 #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
		 ## good pairs of upstream-downstream intron sites:
		 ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
	    
		 ## bad  pairs of upstream-downstream intron sites (they imply wrong strand)
		 ##...###CT...AC###...   ...###GT...AT###...   ...###CT...GC###...
      
		 if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		       ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		       ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
		     print SPLICEUP "dummy".$upstream_exon->dbID. "\t".$upstream_site1."\n";
		     print SPLICEDOWN "dummy".$upstream_exon->dbID."\t".$downstream_site1."\n";
		     $correct++;
		 }
	     } # end of INTRON
	 } # end of strand
	 elsif ( $strand == -1 ){
	     
	     #  example:
	     #                                  ------CT...AC---... 
	     #  transcript in reverse strand -> ######GA...TG###... 
	     # we calculate AC in the slice and the revcomp to get GT == good site
	     
	   INTRON:
	     for (my $i=0; $i<$#exons; $i++ ){
		 my $upstream_exon   = $exons[$i];
		 my $downstream_exon = $exons[$i+1];
		 my $upstream_site;
		 my $downstream_site;
		 my $up_site;
		 my $down_site;
	    
		 my $upstream_site1;
		 my $downstream_site1;
		 my $up_site1;
		 my $down_site1;
		 
		 eval{
		     $up_site1 = 
			 $slice->subseq( ($upstream_exon->start - 10), ($upstream_exon->start + 3) );
		     $down_site1 = 
			 $slice->subseq( ($downstream_exon->end - 3), ($downstream_exon->end + 10 ) );
		 };
		 unless ( $up_site1 && $down_site1 ){
		     print STDERR "problems retrieving sequence for splice sites\n$@";
		     next INTRON;
		 }
		 ( $upstream_site1   = reverse(  $up_site1  ) ) =~ tr/ACGTacgt/TGCAtgca/;
		 ( $downstream_site1 = reverse( $down_site1 ) ) =~ tr/ACGTacgt/TGCAtgca/;
		 

		 eval{
		     $up_site = 
			 $slice->subseq( ($upstream_exon->start - 2), ($upstream_exon->start - 1) );
		     $down_site = 
			 $slice->subseq( ($downstream_exon->end + 1), ($downstream_exon->end + 2 ) );
	    };
		 unless ( $up_site && $down_site ){
		     print STDERR "problems retrieving sequence for splice sites\n$@";
		     next INTRON;
		 }
		 ( $upstream_site   = reverse(  $up_site  ) ) =~ tr/ACGTacgt/TGCAtgca/;
		 ( $downstream_site = reverse( $down_site ) ) =~ tr/ACGTacgt/TGCAtgca/;
		 
		 #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
		 if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		       ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		       ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
		     $correct++;
		     
		     print SPLICEUP "dummy".$upstream_exon->dbID. "\t".$upstream_site1."\n";
		     print SPLICEDOWN "dummy".$upstream_exon->dbID."\t".$downstream_site1."\n";
		     
		 }
		 
	     } # end of INTRON
	 } # end of strand
     
     }# end of transcript

}# end of GENE








