#! /usr/local/bin/perl

use strict;
use ncRNA_update_config;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Getopt::Long;

# This code is intended to transfer the ncRNA genes into the staging server db
# it deletes old ncRNAs
# copies new ones with stable ids and xrefs
# it transfers stable ids and keeps track of the stable_id mappping data

my $pass;
my $write;
my $delete;
my $final_dbname;
my $final_port = 3306;
my $final_host;
my $species;
my $list;
my $xrefs;
my $skipchecks;
my $biotype;
my @new_ncRNAs;
my @old_ncRNAs;
my $coding_overlaps;
my $noncoding_overlaps;
my $noncoding_nonoverlapping;
my $sids;
my $whitelist;
my $count;
my $fbs;
my $dump;
my $no_stable_ids;

$| = 1;

&GetOptions(
	    'pass=s'     => \$pass,
	    'write!'     => \$write,
	    'delete!'    => \$delete,
	    'dbname=s'   => \$final_dbname,
	    'dbhost=s'   => \$final_host,
	    'dbport=s'   => \$final_port,
	    'species=s'  => \$species,
	    'whitelist=s'=> \$list,
	    'xrefs=s'    => \$xrefs,
	    'biotype=s'  => \$biotype,
	    'stable:s'   => \$sids,
	    'slice_fetch!'=> \$fbs,
	    'dump!'      => \$dump,
	    'no_ids!'    => \$no_stable_ids,
	   );

die("transfer_ncRNAs\n-pass *\n-write \n-delete \n-dbname *(final db) \n-dbhost * \n-dbport * \n-species *(1 at at time)
-xrefs (file to dump xref data in) 
-whitelist list of ids to keep
-biotype biotype of genes to transfer
-stable* (file to put stable id mapping data in )
-slice_fetch (fetch the genes slice at a time (quicker in some cases)
-dump (skip all the rest and just dump the xrefs)
-no_ids (do the load without any stable ids)
* = essential\n")
  unless ($pass && $final_port && $final_host && $final_dbname );

die ("transfer_ncRNAs need a file to put stable ids in \n")  unless ($sids or $no_stable_ids);
# get whitelist
if ($list){
  open (LIST,$list) or die "Cannot open whitelist file $list\n";
  while (<LIST>) {
    chomp;
    $whitelist->{$_} = 1;
  }
}

# open xref file
if ($xrefs) {
  open (XREFS,">$xrefs") or die "Cannot open xref file $xrefs\n";
}
# stable_id mapping file
if ($sids) {
  open (SIDS,">$sids") or die "Cannot open stable_id file $sids\n";
}
# open the db connections

my $host    = $CONFIG->{$species}->{"WRITEHOST"};
my $user    = 'ensro';
my $dbname  = $CONFIG->{$species}->{"WRITENAME"};
my $port    = $CONFIG->{$species}->{"WRITEPORT"};

print "$species: Using data in $dbname\@$host:$port\n" ; 

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

print "$species: Using data in $final_dbname\@$final_host:$final_port\n" ; 
die ("Cannot find databases ") unless $final_db  && $sdb;
my $final_ga = $final_db->get_GeneAdaptor;
my $sga = $sdb->get_GeneAdaptor;
my $ssa = $sdb->get_SliceAdaptor;
my $saa = $sdb->get_AnalysisAdaptor;
my $final_sa = $final_db->get_SliceAdaptor;
my $final_aa = $final_db->get_AnalysisAdaptor;
my $analysis;
# just dump the xrefs
if ($dump){
  print "Skipping loading and dumping xrefs\n";
  dump_xrefs($final_ga) if $xrefs;
  exit;
}
# check that the final db has the external db table loaded
check_exdb($final_db);
check_meta($sdb,$final_db);

print "fetching and lazy-loading new predictions...\n";
my $new_hash = fetch_genes($sga,$biotype,$fbs);
print "\nfetching and lazy-loading old predictions...\n";
my $old_hash = fetch_genes($final_ga,$biotype,$fbs);

print "\nFound ".scalar(keys %$new_hash)." new predictions\n";
print "Found ".scalar(keys %$old_hash)." old predictions\n";

$analysis = $final_aa->fetch_by_logic_name("ncRNA");
unless ($analysis){
  print "$final_dbname needs a ncRNA analysis object, loading one\n";
  $analysis = $saa->fetch_by_logic_name("ncRNA");
  die ("ncRNA analysis not found\n") unless $analysis;
  $final_aa->store($analysis);
}

print "Checks\n";
# checks
# blacklist genes over the end of a seq region...
print "Falls off slice\n";
my $overhangs = drop_overhangs( $new_hash);

# duplicte genes
print "Duplicated\n";
my $duplications = duplicates( $new_hash, $sga );

print "Overlaps\n";
# non-coding overlaps
($noncoding_overlaps,$coding_overlaps) = overlaps($new_hash, $final_sa , $final_ga);

# make blacklist of genes to drop
my $blacklist = blacklist($overhangs,$duplications,$coding_overlaps);

print "Transferring stable ids\n" unless ($no_stable_ids);
# transfer stable_ids
my $mapping_session = stable_id_mapping($noncoding_overlaps,$old_hash,$new_hash,$blacklist) unless ($no_stable_ids);

print "Generating stable ids for new predictions\n" unless ($no_stable_ids);
generate_new_ids($new_hash,$final_ga,$blacklist,$mapping_session) unless ($no_stable_ids);

# delete
delete_genes($old_hash,$final_ga ) if $delete;

# write
write_genes($new_hash,$blacklist,$final_ga) if $write;

# dump xrefs
print "Dumping xrefs...\n" if $xrefs;
dump_xrefs($final_ga) if $xrefs;

exit;



sub check_exdb {
  my($db) = @_;
  # test final db for external db table 
  my $query = "SELECT external_db_id FROM external_db WHERE db_name = 'RFAM'";
  die("Cannot find RFAM in external db table\n") unless sql($query,$db)->[0] == 4200;
  $query = "SELECT external_db_id FROM external_db WHERE db_name = 'miRBase'";
  die("Cannot find miRBase in external db table\n") unless sql($query,$db)->[0]  == 3300;
  return;
}

sub drop_overhangs {
  my ($genes) = @_;
  my %blacklist;
  foreach my $key (keys %$genes) {
    my $gene = $genes->{$key};
    # some genes lie over the end of the seq region, get rid of them
    if ($gene->seq_region_start <= 0 or $gene->seq_region_end > $gene->seq_region_length) {
      print "Dropping ".$gene->dbID." as it falls of the edge of the slice\n";
      $blacklist{$gene->dbID} = 1;
    }
  }
  return \%blacklist;
}

sub duplicates {
  my ($genes,$ga) = @_;
  my %blacklist;
  foreach my $key (keys %$genes) {
    my $gene = $genes->{$key};
    # duplicate genes
    my @duplications = @{$ga->fetch_all_by_Slice($gene->feature_Slice)};
    if (scalar(@duplications > 1)) {
      print "Genes ";
      @duplications = sort {$a->dbID <=> $b->dbID} @duplications;
      print $duplications[0]->dbID." and ";
      for (my $i =1 ; $i < scalar(@duplications) ; $i++) {
	my $dbid = $duplications[$i]->dbID;
	print $dbid;
	$blacklist{$dbid} = 1;
	print " are duplicted keeping ".$duplications[0]->dbID." and dumping the other one \n";
      }
    }
  }
  return \%blacklist;
}

sub overlaps {
  my ($genes , $sa, $ga)  = @_;
  # check for overlaps
  my %coding;
  my %noncoding;
 NCRNA: foreach my $key (keys %$genes) {
    my $gene = $genes->{$key};
    my $slice = $sa->fetch_by_region
      (
       'toplevel',
       $gene->seq_region_name,
       $gene->start,
       $gene->end,
       $gene->strand,
      );
    unless ($slice){
      die("NOs slice found for\n".
	  $gene->seq_region_name,"\n",
	  $gene->start,"\n",
	  $gene->end."\n"
	 );
    }
    my @overlaps = @{$ga->fetch_all_by_Slice_constraint($slice,"seq_region_strand = ".$gene->strand)};
  GENE:  foreach my $overlap (@overlaps) {
      # store the overlapping gene in a hash keyed on the predicted genes dbID
      if ($overlap->analysis->logic_name eq "ncRNA") {
	# just check its one of our non coding genes
	next unless scalar(@{$overlap->get_all_Exons}) == 1;
	next if $overlap->biotype =~ /^Mt_/;
	if($biotype && $overlap->biotype ne $biotype){
	  warn("The non coding gene overlapping yours is of a different type, so it won't get deleted: $biotype vs ".$overlap->biotype." not transferring this stable id \n");
	  next GENE;
	}
	# catch problem where you have multiple overlapping ncRNAs, maybe you can manually delete one
	# so that it transfers the correct stable_id
	if ($noncoding{$gene->dbID}){
	  warn("Something fishy is going on here I have more than one non coding overlap for this gene".
	       $gene->dbID." ".$gene->seq_region_name." ".$gene->start.":".$gene->end.":".$gene->strand.
	       " overlaps\nA: ".$noncoding{$gene->dbID}->dbID." ".$noncoding{$gene->dbID}->biotype."\nB: ".
	       $overlap->dbID." ".$overlap->biotype."\n");
	  print "Do you want to transfer the stable id of A or B?";
	  my $reply = <>;
	  chomp $reply;
	  next GENE if $reply eq "A" or $reply eq "a";
	}
	$noncoding{$gene->dbID} = $overlap;
      } else {
	# overlapping gene is coding
	# want to know if it actually overlaps a coding exon	
	# exon is coding
	foreach my $trans (@{$overlap->get_all_Transcripts}) {
	  foreach my $exon (@{$trans->get_all_translateable_Exons}) {
	    my $codingexon = ($slice->start+$exon->end() >= $gene->start() &&
			    $slice->start+$exon->start() <= $gene->end());
	    if ($codingexon) {
	      print "gene ".$gene->dbID."\t".$gene->description." overlaps coding exon "
		.$exon->stable_id." in  real_gene ".$overlap->stable_id."\t";
	      if ($whitelist->{$gene->dbID}){
		print "Gene protected by whitelist\n";
	      } else {
		$coding{$gene->dbID} = $overlap;
		print "\n";
	      }
	      next GENE;
	    }
	  }
	}
      }
    }
  }
  return (\%noncoding, \%coding);
}

sub blacklist {
  my ($overhangs,$duplications,$coding_overlaps) = @_;
  my %list;
  foreach my $key (keys %$overhangs){
    $list{$key} = 1 unless $whitelist->{$key};
  }
  foreach my $key (keys %$duplications){
    $list{$key} = 1 unless $whitelist->{$key};
  }
  foreach my $key (keys %$coding_overlaps){
    $list{$key} = 1 unless $whitelist->{$key};
  }
  return \%list;
}

sub stable_id_mapping {
  my ($non_coding_overlaps,$old_hash,$new_hash,$blacklist) = @_;
  # get the assembly information
  my $last_session = sql('SELECT max(mapping_session_id) from mapping_session',$final_db)->[0];
  $last_session++;
  my $old_release = sql('SELECT meta_value from meta where meta_key = "schema_version"',$sdb)->[0];
  my $new_release = sql('SELECT meta_value from meta where meta_key = "schema_version"',$final_db)->[0];
  my $old_assembly = sql('SELECT meta_value from meta where meta_key = "assembly.default"',$sdb)->[0];
  my $new_assembly = sql('SELECT meta_value from meta where meta_key = "assembly.default"',$final_db)->[0];
  print SIDS "INSERT INTO  mapping_session(mapping_session_id,old_db_name,new_db_name,old_release,new_release,old_assembly,new_assembly,created) ".
    " VALUES($last_session,\'$dbname\',\'$final_dbname\',\'$old_release\',\'$new_release\',\'$old_assembly\',\'$new_assembly\',now());\n";
  # trasfer them where you have overlaps with non coding genes
  # put the appropriate entries in the stable id mapping table
  my %done;
  foreach my $new_ncRNA_id (keys %$non_coding_overlaps){
    next if $blacklist->{$new_ncRNA_id};
    # fetch gene, exon and transcript objects ( 1 of each)
    my $new_gene =  $new_hash->{$new_ncRNA_id};
    my $new_trans = $new_gene->get_all_Transcripts->[0];
    my $new_exon = $new_trans-> get_all_Exons->[0];
    my $old_gene = $non_coding_overlaps->{$new_ncRNA_id};
    my $old_trans = $old_gene->get_all_Transcripts->[0];
    my $old_exon = $old_trans-> get_all_Exons->[0];
    $new_gene->stable_id($old_gene->stable_id);
    $new_gene->version($old_gene->version);
    $new_trans->stable_id($old_trans->stable_id);
    $new_trans->version($old_trans->version);
    $new_exon->stable_id($old_exon->stable_id);
    $new_exon->version($old_exon->version);
    my $version = $old_trans->version;
    $version++;
    # transferring
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES(\'".
      $new_gene->stable_id."\',".$old_gene->version.",\'".$new_gene->stable_id."\',$version,$last_session,\'gene\',0);\n";
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES(\'".
      $new_trans->stable_id."\',".$old_trans->version.",\'".$new_trans->stable_id."\',$version,$last_session,\'transcript\',0);\n";
    $done{$old_gene->dbID} = 1;
  }
  # deleting
  foreach my $old_gene (keys %$old_hash){
    next if $done{$old_gene};
    my $gene =  $old_hash->{$old_gene};
    my $trans = $gene->get_all_Transcripts->[0];
    # old gene does not need stable id transferring stable id is dead 
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES('".
      $gene->stable_id."',".$gene->version.",null,0,$last_session,'gene',0);\n";
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES('".
      $trans->stable_id."',".$trans->version.",null,0,$last_session,'transcript',0);\n";
    # need gene archive entries also...
    print SIDS "INSERT INTO gene_archive(gene_stable_id,gene_version,transcript_stable_id,transcript_version,translation_stable_id,translation_version,peptide_archive_id,mapping_session_id) ";
    print SIDS "VALUES('". $gene->stable_id."',".$gene->version.",'". $trans->stable_id."',".$trans->version.",'',0,0,$last_session);\n";
  }
  # new stable ids need to be added with the appropirate code...
  return $last_session;
}

sub generate_new_ids{
  my ($ncRNAs,$ga,$list,$last_session) = @_;
  my ($gsp,$gsi,$tsp,$tsi,$esp,$esi);
  # ensembl_ids all have 11 numbers at the end and an indeterminate number of letters at the start
  
  if (sql("SELECT max(stable_id) from gene_stable_id ;",$ga)->[0] =~ /^(\D+)(\d+)$/){
    $gsp = $1;
    $gsi = $2;
  }
  if ( sql("SELECT max(stable_id) from transcript_stable_id ;",$ga)->[0] =~ /^(\D+)(\d+)$/){
    $tsp = $1;
    $tsi = $2;
  }
  if ( sql("SELECT max(stable_id) from exon_stable_id ;",$ga)->[0] =~ /^(\D+)(\d+)$/){
    $esp = $1;
    $esi = $2;
  }
  # check stable ids of dead genes are not higher than the maximum in the gene stable id table
  if (sql("SELECT max(gene_stable_id) from gene_archive ;",$ga)->[0] =~ /^(\D+)(\d+)$/){
    if ($2 > $gsi){
      print "dead gene with higher id $gsp$gsi\n";
      $gsp = $1;
      $gsi = $2;
      print " becomes $gsp$gsi\n";
    }
  }
  if ( sql("SELECT max(transcript_stable_id) from gene_archive ;",$ga)->[0] =~ /^(\D+)(\d+)$/){
    if ($2 > $tsi){
      print "dead trans with higher id $tsp$tsi";
      $tsp = $1;
      $tsi = $2;
      print " becomes $tsp$tsi\n";
    }
  }

    unless ($gsp && $gsi && $tsp && $tsi && $esp && $esi){
    print "Cannot figure out how to make new stable ids\nHave got :";
    print sql("SELECT max(stable_id) from gene_stable_id ;",$ga)->[0];
    print " ".sql("SELECT max(stable_id) from transcript_stable_id ;",$ga)->[0];
    print " ".sql("SELECT max(stable_id) from exon_stable_id ;",$ga)->[0]."\n";;
    die();
  }
  foreach my $ncRNA_id (keys %$ncRNAs){
    next if $list->{$ncRNA_id};
    my $gene = $ncRNAs->{$ncRNA_id};
    next if $gene->stable_id;
    my $trans = $gene->get_all_Transcripts->[0];
    my $exon = $trans->get_all_Exons->[0];
    $gsi++;
    $tsi++;
    $esi++;
    $gene->stable_id($gsp.$gsi);
    $trans->stable_id($tsp.$tsi);
    $exon->stable_id($esp.$esi);
    $gene->version(1);
    $trans->version(1);
    $exon->version(1);
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES(";
    print SIDS "null,0,'".$gene->stable_id."',1,$last_session,'gene',0);\n";
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES(";
    print SIDS "null,0,'".$trans->stable_id."',1,$last_session,'transcript',0);\n";
  }
  return;
}


sub delete_genes {
  my ($old_hash,$ga) = @_;
  return if scalar(keys %$old_hash) == 0;
  print STDERR "Warning you have not specified a biotype *ALL* ncRNAs will be deleted\n" unless ($biotype);
  print STDERR "Found ".scalar(keys %$old_hash)." genes  $final_dbname\nshall I delete them? (Y/N) ";
  my $reply = <>;
  chomp $reply;
  if ($reply eq "Y" or $reply eq "y") {
    foreach my $key (keys %$old_hash){
      my $gene = lazy_load($old_hash->{$key});
      unless ($gene->biotype =~ /RNA/ && $gene->analysis->logic_name eq 'ncRNA'){
	throw("Gene to be deleted is not a non coding gene ".$gene->dbID."\t".$gene->stable_id."\t".$gene->biotype."\n");
      }
      print "Deleting gene ".$gene->dbID."\t".$gene->stable_id."\t".$gene->biotype."\n";
      $ga->remove($gene);
    }
  }
  return;
}

sub write_genes {
  my  ($new_hash,$blacklist,$ga) = @_;
  foreach my $key (keys %$new_hash){
    if ($blacklist->{$key}){
      print "Skipping $key - blacklisted\n";
      next;
    }
    my $gene = lazy_load($new_hash->{$key});
    unless ($gene->biotype =~ /RNA/ && ( $gene->analysis->logic_name eq 'ncRNA' or $gene->analysis->logic_name eq 'miRNA') ){
      throw("Gene to be written is not a non coding gene ".$gene->dbID."\t".$gene->stable_id."\t".$gene->biotype."\n");
    }
    # copy the analysis
    $gene->analysis($analysis);
    my $trans = $gene->get_all_Transcripts->[0];
    $trans->analysis($analysis);
    $trans->biotype($gene->biotype);
    $trans->status($gene->status);
    print "Storing gene ".$gene->dbID."\t".$gene->stable_id."\t".$gene->biotype."\n" ;
    $ga->store($gene);
  }
  return;
}

sub lazy_load {
  my ($gene) = @_;
  $gene->stable_id;
 $gene->get_all_Transcripts->[0]->stable_id;
 $gene->get_all_Exons->[0]->stable_id;
 $gene->get_all_DBEntries;
 $gene->get_all_Transcripts->[0]->get_all_supporting_features;
 return $gene;
}

sub dump_xrefs {
  my ($ga)= @_;
  # get all non coding genes to dump xrefs of
  my $gene_hash = fetch_genes($ga,undef,$fbs);
  foreach my $key (keys %$gene_hash){
    my $gene = $gene_hash->{$key};
    next if $gene->biotype =~ /Mt_/;
    next unless $gene->analysis->logic_name eq 'ncRNA';
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      my @xrefs = @{$trans->get_all_DBEntries};
      if (@xrefs){
	foreach my $xref (@xrefs) {
	  next unless ($xref->dbname eq 'miRNA_Registry' or $xref->dbname eq 'RFAM');
	  print XREFS $gene->dbID."\t";
	  print XREFS $trans->dbID."\t";
	  print XREFS $xref->dbname."\t"; 
	  print XREFS $xref->primary_id."\t";
	  print XREFS $xref->display_id."\t";
	  print XREFS $gene->description."\t";
	  print XREFS $gene->status."\n";
	}
      }
    }
  }
  return;
}

sub sql {
  my ($query,$db) = @_;
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  my @array = $sth->fetchrow_array;
  return \@array;
}

sub fetch_genes {
  my ($ga,$biotype,$slice) = @_;
  my %ncRNA_hash;
  my @ncRNAs;
  throw("Cannot fetch genes without gene adaptor $ga") unless $ga;
  if ($slice){
    my @slices = @{$final_sa->fetch_all('toplevel')};
    my $inc = scalar(@slices) / 20;
    print STDERR "|------------------|\r|";
    foreach my $slice (@slices){
      $count++;
      if ($count >= $inc){
	$count = 0;
	print STDERR "=";
      }
      if ($biotype) {
	@ncRNAs =  @{$ga->fetch_all_by_Slice_constraint($slice,"biotype = '".$biotype."'")};
      } else {
	@ncRNAs =  @{$ga->fetch_all_by_Slice_constraint($slice,'biotype like "%RNA"')};
      }
      foreach my $ncRNA (@ncRNAs) {
	next unless ($ncRNA->analysis->logic_name eq 'ncRNA' or $ncRNA->analysis->logic_name eq 'miRNA') ;
	next if  $ncRNA->biotype =~ /Mt_/;
	next if  $ncRNA->description =~ /RNAI/;	
	$ncRNA_hash{$ncRNA->dbID} = lazy_load($ncRNA);
      }
    }
  } else {
      if ($biotype) {
	@ncRNAs =  @{$ga->generic_fetch("biotype = '".$biotype."'")};
      } else {
	@ncRNAs =  @{$ga->generic_fetch('biotype like "%RNA"')};
      }      
      foreach my $ncRNA (@ncRNAs) {
	next unless ($ncRNA->analysis->logic_name eq 'ncRNA' or $ncRNA->analysis->logic_name eq 'miRNA') ;
	next if  $ncRNA->biotype =~ /Mt_/;
	next if  $ncRNA->description =~ /RNAI/;
	$ncRNA_hash{$ncRNA->dbID} = lazy_load($ncRNA);
      }
    }
  print STDERR  "\n";
  return \%ncRNA_hash;
}

sub check_meta {
  my ($db1,$db2) = @_;
  my $m1 = $db1->get_MetaContainer;
  my $m2 = $db2->get_MetaContainer;
  my $c1 = sql("SELECT meta_value from meta where meta_key = 'assembly.default'",$db1)->[0];
  my $c2 = sql("SELECT meta_value from meta where meta_key = 'assembly.default'",$db2)->[0];

  unless ($m1->get_Species->common_name eq $m2->get_Species->common_name){
    throw("Species do not agree ".$m1->get_Species->common_name." != ". $m2->get_Species->common_name."\n");
  }
  unless ($m1->get_taxonomy_id eq $m2->get_taxonomy_id ){
    throw("Tax ids do not agree ".$m1->get_taxonomy_id." != ". $m2->get_taxonomy_id."\n");
  }
  unless ($c1 eq $c2){
    throw("Coord systems do not agree $c1 != $c2 \n");
  }
  return;
}


__END__
