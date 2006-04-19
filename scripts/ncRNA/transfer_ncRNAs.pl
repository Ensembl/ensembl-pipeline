#! /usr/local/bin/perl

use strict;
use ncRNA::ncRNA_update_config;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Getopt::Long;

# some stuff to think about....
# want to delete ncRNAs from databases where they already exist
# want to merge genes with human (chicken?)
# dont really want to do this automatically incase it goes mad esp since you are deleting ncRNA genes...
# best way to delete the ncRNA genes is to load them and then remove using the api

# so lets use the existing config info to get the species so there is no room for error there
# but go one at a time and enter the db to write / delete from / to by hand

my $pass;
my $verbose;
my $write;
my $delete;
my $merge;
my $final_dbname;
my $final_port;
my $final_host;
my @species_list;
my $logic_name;
my $whitelist;
my %whitehash;
my $xrefs;
my $skipchecks;
my $increment;

$| = 1;

&GetOptions(
	    'pass=s'     => \$pass,
	    'write!'     => \$write,
	    'verbose!'   => \$verbose,
	    'delete!'    => \$delete,
	    'merge!'     => \$merge,
	    'dbname=s'   => \$final_dbname,
	    'dbhost=s'   => \$final_host,
	    'dbport=s'   => \$final_port,
	    'species=s'  => \@species_list,
	    'analysis=s' => \$logic_name,
	    'whitelist=s'=> \$whitelist,
	    'xrefs=s'    => \$xrefs,
	    'skipchecks!'=> \$skipchecks,
	    'increment!' => \$increment,
	   );

print "transfer_ncRNAs\n-pass *\n-write \n-verbose \n-delete \n-merge \n-dbname *(final db) \n-dbhost * \n-dbport * \n-species *(1 at at time)
-analysis *(analysis of genes to delete ie: ncRNA) \n-xrefs (file to dunp xref data in) \n-skipchecks (dont do checks) \n-increment (just add ncRNAs where they are new)
-whitelist list of ids to keep
\n* = essential\n" 
  unless ($pass && $final_port && $final_host && $final_dbname);

if ($whitelist) {
  print "Loading Whitelist of genes to keep\n" if $verbose;
  open (LIST,$whitelist) or die "Cannot open whitelist file $whitelist\n";
  while (<LIST>) {
    chomp;
    $whitehash{$_} = 1;
  }
}
print "Looking at the config file\n" if $verbose;;
if ($xrefs) {
  open (XREFS,">$xrefs") or die "Cannot open xref file $xrefs\n";
}


my @speciess;
if (scalar(@species_list)) {
  @species_list = split(/,/,join(',',@species_list));
  foreach my $species (@species_list) {
    if ($CONFIG->{$species}) {
      push @speciess, $species;
    } else {
      print "Skipping species $species\n";
    }
  }
} else {
  @speciess = (keys %$CONFIG);
}

die "Only 1 species at a time please\n" unless (scalar(@speciess) == 1) ;


foreach my $species (@speciess) {

  my $host    = $CONFIG->{$species}->{"WRITEHOST"};
  my $user    = 'ensro';
  my $dbname  = $CONFIG->{$species}->{"WRITENAME"};
  my $port    = $CONFIG->{$species}->{"WRITEPORT"};
    
  # This needs to be the final core db you want to write the genes to
  # could be the genebuild database or an existing database on 3364....hmmm...
    
  my @overlapping_genes;
  my @genes_to_copy;
  my @duplicted;
  my $total;
  my %blacklist;
    
  my $sdb = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     -host   => $host,
     -user   => $user,
     -port   => $port,
     -dbname => $dbname,
    );
    
    
  my $final_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     -host   => $final_host,
     -user   => 'ensadmin',
     -port   => $final_port,
     -dbname => $final_dbname,
     -pass   => $pass,
    );
  die ("Cannot find databases ") unless $final_db  && $sdb;
  my $final_ga = $final_db->get_GeneAdaptor;
  my $sga = $sdb->get_GeneAdaptor;
  my $ssa = $sdb->get_SliceAdaptor;
  my $final_sa = $final_db->get_SliceAdaptor;
  my $total_ncRNAs;
  my @ncRNAs;
  my %failed;
    
  print "Fetching Predicted Genes\n" if $verbose;
  my @predicted_genes = @{$sga->fetch_all};
    
  my ($desc,$str,$mat)=0;

  unless ($skipchecks){
    print "Tests\n" if $verbose;;
    foreach my $gene (@predicted_genes) {
      # some genes lie over the end of the seq region, get rid of them
      if ($gene->seq_region_start <= 0 or $gene->seq_region_end > $gene->seq_region_length){
	print "Dropping ".$gene->dbID." as it falls of the edge of the slice\n";
	$blacklist{$gene->dbID} = 1;
      }
      # gene must have description
      unless ($gene->description =~ /\S+\s+\[Source.+\]/){
	$failed{$gene->dbID."\tdescription"} = $gene;
	$desc++;
      }
      # gene must have full length structure line
      # miRNAs must have a mature coord attributes
      foreach my $trans (@{$gene->get_all_Transcripts}) {
	my $attributes = $trans->get_all_Attributes;
	my $structure = &decode_str($attributes);
	unless ($trans->length == length($structure)){
	  $failed{$gene->dbID."\tstructure"} = $gene;
	  $blacklist{$gene->dbID} =1;
	  $str++;
	}
	if ($gene->analysis->logic_name eq 'miRNA') {
	  my $mature = 0;
	  foreach my $attrib (@{$attributes}) {
	    $mature = 1 if $attrib->name eq 'Micro RNA';
	  }
	  unless ($mature){
	    $failed{$gene->dbID."\tmature attribute"} = $gene;
	    $mat++;
	  }
	}
      }
      # duplicate genes
      my @duplications = @{$sga->fetch_all_by_Slice($gene->feature_Slice)};
      if (scalar(@duplications > 1)) {
	print "Genes ";
	@duplications = sort {$a->dbID <=> $b->dbID} @duplications;
	print $duplications[0]->dbID." and ";
	for (my $i =1 ; $i < scalar(@duplications) ; $i++) {
	  my $dbid = $duplications[$i]->dbID;
	  print $dbid;
	  $blacklist{$dbid} = 1;
	  print " are duplicted keeping ".$duplications[0]->dbID." and dumping the other one \n" if $verbose;;
	}
      }
    }
    foreach my $gene (keys %failed) {
      print "$gene\t".$failed{$gene}->description."\n"  if ($gene =~/description/ &&   $verbose);
      print "$gene\tstructure is wrong \n" if ($gene =~/structure/ &&   $verbose);
      print "$gene\tno mature coords \n" if ($gene =~/mature/   &&   $verbose);
    }
 

    print "All descriptions pass tests\n" unless($desc) ;
    print "All structures pass tests\n" unless($str);
    print "All mature miRNA attributes pass tests\n" unless($mat);
  }
  my $analysis = $final_db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
  # if the analysis exists see if you can get hold of any ncRNA genes....
  if ($analysis && $delete) {
    my @old_ncRNAs = @{$final_ga->generic_fetch(" analysis_id = ".$analysis->dbID)};
    if (scalar(@old_ncRNAs) > 0) {
      print "Found ".scalar(@old_ncRNAs)." genes of type ".$analysis->logic_name." from $final_dbname\nshall I delete them? (Y/N) ";
      my $reply = <>;
      chomp $reply;
      if ($reply eq "Y" or $reply eq "y") {
	# do the delete
	foreach my $old_ncRNA (@old_ncRNAs) {
	  eval {
	    $final_ga->remove($old_ncRNA);
	  };
	  if ($@) {
	    print "Error removing gene ".$old_ncRNA->dbID."\n$@\n" if $verbose;;
	  }
	}
      }
    }
  }

  

  unless ($skipchecks){
    print "looking for overlaps with coding exons and ncRNAs\n" if $verbose;;
    # check for overlaps
    NCRNA : foreach my $pg (@predicted_genes){
      my $slice = $final_sa->fetch_by_region
	(
	 'toplevel',
	 $pg->seq_region_name,
	 $pg->start,
	 $pg->end,
	 $pg->strand,
	);
      unless ($slice){
	print "NOs slice found for\n".
	  $pg->seq_region_name,"\n",
	    $pg->start,"\n",
	      $pg->end."\n" if $verbose;
	die "HELP\n";
      }
      my @real_genes = @{$final_ga->fetch_all_by_Slice_constraint($slice)};
      
    GENE:  foreach my $real_gene (@real_genes) {
	# gene is a non coding gene already annotated - delete the old one
	if ($merge or $increment) {
	  if ($real_gene->analysis->logic_name eq "ncRNA" && $real_gene->strand == 1) {
	    # delete the gene
	    if ($merge){
	      print "Replacing ".$real_gene->dbID."\t" if $verbose;
	      $final_ga->remove($real_gene) if $merge;
	    }
	    if ($increment){
	      print "Keeping overlapping gene ".$real_gene->dbID."\n" if $verbose;
	      push @overlapping_genes,$pg;
	      next NCRNA;
	    }
	    if ($real_gene->stable_id) {
	      print "gene ".$pg->dbID."\t".$pg->description." overlaps real_gene ".$real_gene->stable_id.
		" in gene ".$real_gene->stable_id."\n" if $verbose;
	    } else {
	      print "gene ".$pg->dbID."\t".$pg->description." overlaps real_gene ".$real_gene->dbID.
		" in gene ".$real_gene->dbID."\n" if $verbose;
	    }
	  }
	}
	# exon is coding
	foreach my $trans (@{$real_gene->get_all_Transcripts}) {
	  foreach my $exon (@{$trans->get_all_translateable_Exons}) {
	    my $overlaps = ($slice->start+$exon->end() >= $pg->start() &&
			    $slice->start+$exon->start() <= $pg->end());
	    if ($overlaps) {
	      if ($exon->stable_id) {
		print "gene ".$pg->dbID."\t".$pg->description." overlaps exon ".$exon->stable_id.
		  " in gene ".$real_gene->stable_id."\n" if $verbose;;
	      } else {
		print "gene ".$pg->dbID."\t".$pg->description." overlaps exon ".$exon->dbID.
		  " in gene ".$real_gene->dbID."\n" if $verbose;
	      }
	      unless ($whitehash{$pg->dbID}) {
		push @overlapping_genes,$pg;
		next NCRNA;
	      } else {
		print $pg->dbID." is protected by the whitelist - keeping it\n" if $verbose;;
	      }		
	    }
	  }
	}	
      }
      if ($blacklist{$pg->dbID}) {
	push @duplicted,$pg;
      } else {
	push @genes_to_copy,$pg;
      }
    }
  } else {
    @genes_to_copy = @predicted_genes;
  }

  print scalar(@overlapping_genes) if $verbose;;
  print " overlap " if $verbose;
  print scalar(@duplicted) if $verbose;;
  print " are duplicated " if $verbose;;
  print scalar(@genes_to_copy) if $verbose;;
  print " are therefore availible for copying\n" if $verbose;;

# foreach my $g (@genes_to_copy){
#    print $g->description."\n";
#  }

  if ($write) {
    print "Copying genes to ".$final_dbname."\n" if $verbose;;
    my $ncRNA_analysis = $final_db->get_AnalysisAdaptor->fetch_by_logic_name("ncRNA");
    unless ($ncRNA_analysis){
      $ncRNA_analysis = Bio::EnsEMBL::Analysis->new
	(
	 -logic_name      => 'ncRNA',
	 -gff_source      => 'ensembl',
	 -gff_feature     => 'gene',
	);
      $final_db->get_AnalysisAdaptor->store($ncRNA_analysis);
    }
    foreach my $gene (@genes_to_copy) {
      # lazyload
      foreach my $trans (@{$gene->get_all_Transcripts}) {
	foreach my $exon (@{$trans->get_all_Exons}) {
	}
	foreach my $xref (@{$trans->get_all_DBEntries}) {
	  #	  print "$gene\t$trans\t$xref\n";
	  #	  print  $xref->description."\n";
	}
      }
      # adjust  analysis
      $gene->analysis($ncRNA_analysis);
      $gene->dbID(undef);
      $gene->adaptor(undef);
      #Biotypes
      print "Storing gene ".$gene->dbID."\t" if $verbose;;
      # store gene
      eval {
	my $new_id = $final_ga->store($gene);
	if ($xrefs) {
	  # dump out all the xref info
	  my $new_gene = $final_ga->fetch_by_dbID($new_id);
	  unless ($new_gene){
	    print "Warning newly written gene not found $new_id\n" if $verbose;;
	  }
	  next unless($gene->analysis->logic_name eq "ncRNA");
	  print XREFS "$new_id\t";
	  foreach my $trans (@{$new_gene->get_all_Transcripts}) {
	    print XREFS $trans->dbID."\t";
	    foreach my $xref (@{$trans->get_all_DBEntries}) {
	      print XREFS $xref->dbname."\t"; 
	      print XREFS $xref->primary_id."\t";
	      print XREFS $xref->display_id."\t";
	      print XREFS $new_gene->description."\n";
	    }
	  }
	}
	print   "wrote gene " . $gene->dbID . "\n" if $verbose;;
      };
      if ( $@ ) {
	die("UNABLE TO WRITE GENE:\n$@");
      }
    }
  } else {
    print "Not writing - write protect on!\n" if $verbose;;
  }
}


exit;

# sub for decoding gene structure

sub decode_str{
  my ($attributes)= @_;
  my $code;
  my $offset = 0;
  my $chr;
  my $string;
  my $start;
  my $end;
  my $num = 0;
  my $str = "";
  my $complete_str;
  my %str_hash;
  foreach my $attribute (@$attributes){
    next unless $attribute->name eq "Structure";
    my $value = $attribute->value;
    # remove string header;
    if ($value =~ /(\d+):\d+\t(.+)/){
      $str_hash{$1} = $2;
    } else {
      throw("Cannot parse encoded string $attribute\n");
    }
  }
  my @sorted_attributes = sort { $str_hash{$a} cmp $str_hash{$b} } keys %str_hash;
  foreach my $order (@sorted_attributes){
    $string = $str_hash{$order};
    $num =0;
    $str = "";
    for (my $pos =0; $pos <= length($string) ; $pos++){
      $chr =  substr($string,$pos,1);
      if ($chr =~ /\D/){
	$complete_str .= $str unless($num);
	if ($num){
	  for (my$i =1 ; $i <= $num ; $i++){
	    $complete_str .= $str;
	  }
	}
	$str = $chr;
	$num = "";
	next;
      }
      if ($chr =~ /\d/){
	$num .= $chr;
      }
    }
    # empty array 
    if ($num){
    for (my$i =1 ; $i <= $num ; $i++){
      $complete_str .= $str;
    }
  }
    else{
      $complete_str .= $str;
    }
  }
return $complete_str;

}

__END__
