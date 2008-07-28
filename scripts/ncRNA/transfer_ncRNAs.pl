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
# adding a merge option for chicken and human

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
my $biotype_to_skip;
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
my $merge;
my $merge_set;
#list of gene-descriptions to ignore
my @genestoignore = ();
my $makewhitelist = 0;
my $new_release;
my $use_old_ncRNAs;

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
	    'biotype_to_skip=s'  => \$biotype_to_skip,
	    'stable:s'   => \$sids,
	    'slice_fetch!' => \$fbs,
	    'dump!'      => \$dump,
	    'no_ids!'    => \$no_stable_ids,
	    'makewhitelist!' => \$makewhitelist,
	    'release:s'  => \$new_release,
	    'merge!'     => \$merge,
	    'use_old_ncRNAs!'  => \$use_old_ncRNAs
	   );

die("transfer_ncRNAs\n-pass $pass  * 
-write $write
-delete $delete
-dbname  $final_dbname *(final db)
-dbhost $final_host*
-dbport $final_port * 
-species $species *(1 at at time)
-xrefs $xrefs(file to dump xref data in) 
-whitelist $list list of ids to keep
-biotype $biotype biotype of genes to transfer
-biotype_to_skip  $biotype_to_skip biotype of genes to not transfer
-stable $sids (file to put stable id mapping data in )
-slice_fetch $fbs (fetch the genes slice at a time (quicker in some cases)
-dump $dump (skip all the rest and just dump the xrefs)
-no_ids $no_stable_ids (do the load without any stable ids)
-release $new_release ( the number of the release ie 44 )*
-merge ( special case for human where there are Sean Eddys genes that we want to keep except for where we have a better prediction )
-use_old_ncRNAs  Fetch the old ncRNAs from the dna db on livemirror ( usually ) useful if the core db on staging has no ncRNAs in it
* = essential\n")
  unless ($pass && $final_port && $final_host && $final_dbname &&  $new_release);

die ("transfer_ncRNAs need a file to put stable ids in \n")  unless ($sids or $no_stable_ids);
# get whitelist
if ($list){
  open (LIST,$list) or die "Cannot open whitelist file $list\n";
  while (<LIST>) {
    chomp;
    $whitelist->{$_} = 1;
  }
}
if($makewhitelist){
  my $whitelist_file = $species.".whitelist";
  open(WHITE, ">$whitelist_file") or die "cant create whitelist file $whitelist_file.\n";
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
my $dnahost    = $CONFIG->{$species}->{"DBHOST"};
my $dnadbname  = $CONFIG->{$species}->{"DBNAME"};
my $dnaport    = $CONFIG->{$species}->{"DBPORT"};


print "$species: Using data in $dbname\@$host:$port\n" ; 

my $sdb = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $host,
   -user   => $user,
   -port   => $port,
   -dbname => $dbname,
  );

# using the staging server db as the dna db because it has the same assembly right?
# and the server is quieter generally

# old data base on livemirror - useful to have sometimes
my $olddb = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $dnahost,
   -user   => 'ensro',
   -port   => 3306,
   -dbname => $dnadbname,
  );


my $final_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $final_host,
   -user   => 'ensadmin',
   -port   => $final_port,
   -dbname => $final_dbname,
   -pass   => $pass,
  );
throw("No dna database found\n") unless $final_db;

$sdb->dnadb($final_db);

print "$species: Using data in $final_dbname\@$final_host:$final_port\n" ; 
die ("Cannot find databases ") unless $final_db  && $sdb;
my $final_ga = $final_db->get_GeneAdaptor;
my $old_ga =   $olddb->get_GeneAdaptor;
my $sga = $sdb->get_GeneAdaptor;
my $ssa = $sdb->get_SliceAdaptor;
my $saa = $sdb->get_AnalysisAdaptor;
my $final_sa = $final_db->get_SliceAdaptor;
my $old_sa = $olddb->get_SliceAdaptor;
my $old_ga = $olddb->get_GeneAdaptor;
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

print "fetching and lazy-loading new predictions from $dbname @ $host....\n";
my $new_hash = fetch_genes($sga,$biotype,$biotype_to_skip,$fbs,\@genestoignore,$ssa);
my $old_hash;
if ($use_old_ncRNAs){
  print "Fetching old ncRNAs from previous release $dnadbname @ $dnahost....\n";
  $old_hash = fetch_genes($old_ga,$biotype_to_skip,$biotype,$fbs,\@genestoignore,$old_sa);
} else {
  print "Fetching old ncRNAs from previous release $final_dbname @ $final_host....\n";
  $old_hash = fetch_genes($final_ga,$biotype,$biotype_to_skip,$fbs,\@genestoignore,$final_sa);
}

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

# genes with high AT content are blacklisted because they are just wrong...
my $repeats = at_content( $new_hash );

print "Overlaps\n";
# non-coding overlaps
if ($use_old_ncRNAs){
($noncoding_overlaps,$coding_overlaps,$merge_set) = overlaps($new_hash, $old_sa , $old_ga);
} else {
($noncoding_overlaps,$coding_overlaps,$merge_set) = overlaps($new_hash, $final_sa , $final_ga);
}

# make blacklist of genes to drop
my $blacklist = blacklist($overhangs,$duplications,$coding_overlaps,$repeats);

print "Transferring stable ids\n" unless ($no_stable_ids);
# transfer stable_ids
my $mapping_session = stable_id_mapping($noncoding_overlaps,$old_hash,$new_hash,$blacklist,$merge_set) unless ($no_stable_ids);

print "Generating stable ids for new predictions\n" unless ($no_stable_ids);
generate_new_ids($new_hash,$final_ga,$blacklist,$mapping_session) unless ($no_stable_ids);

# delete
delete_genes($old_hash,$final_ga,$merge_set ) if $delete;
 
# need to clear out the xrefs as miRBase complicted things by changing the xref names without 
# changing the accessions which causes the API to go a little nuts
print "Deleting old unused ncRNA xrefs\n" if $write;
delete_sql("DELETE xref FROM xref LEFT JOIN object_xref ON xref.xref_id = object_xref.xref_id
WHERE external_db_id in(3300, 4200) and object_xref.xref_id is null",$final_db) if $write;

# write
write_genes($new_hash,$blacklist,$final_ga) if $write;

# dump xrefs
print "Dumping xrefs...\n" if $xrefs;
dump_xrefs($final_ga) if $xrefs;

close(WHITE) if($makewhitelist);

unless($no_stable_ids or !$sids){
  print STDERR "\n# mysql -uensadmin -p".$pass." -h".$final_host." -P".$final_port." -D".$final_dbname." < ".$sids."\n\n";
}

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
    foreach my $exon ( @{$gene->get_all_Exons} ) {
      # some genes lie over the end of the seq region, get rid of them
      if ($exon->seq_region_start <= 0 or $exon->seq_region_end > $gene->seq_region_length) {
	print "Dropping ".$gene->dbID." as it falls of the edge of the slice\n";
	$blacklist{$gene->dbID} = 1;
      }
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
    my @duplications = @{$ga->fetch_all_by_Slice_constraint($gene->feature_Slice,'biotype like "%RNA%" ')};
    if (scalar(@duplications > 1)) {
      print "Genes ";
      @duplications = sort {$a->dbID <=> $b->dbID} @duplications;
      print $duplications[0]->dbID." and ";
      for (my $i =1 ; $i < scalar(@duplications) ; $i++) {
	my $dbid = $duplications[$i]->dbID;
	print $dbid;
	$blacklist{$dbid} = 1;
	print " are duplicated keeping ".$duplications[0]->dbID." and dumping the other one \n";
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
  my %merge_set;
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
      warn("NO slice found for\n".
	  $gene->seq_region_name,"\n",
	  $gene->start,"\n",
	  $gene->end."\n"
	 );
      $coding{$gene->dbID} = 1;
      next NCRNA;
    }
    my @overlaps = @{$ga->fetch_all_by_Slice_constraint($slice,"seq_region_strand = ".$gene->strand)};
  GENE:  foreach my $overlap (@overlaps) {
      # store the overlapping gene in a hash keyed on the predicted genes dbID
      if ($overlap->analysis->logic_name eq "ncRNA") {
	# just check its one of our non coding genes
	next unless scalar(@{$overlap->get_all_Exons}) == 1;
	next if $overlap->biotype =~ /^Mt_/;
	# used in a gene merge to remove overlapping ncRNAs from Sean Eddys set
	$merge_set{$overlap->dbID} = 1 unless $overlap->description =~ /RFAM/;
	if( $overlap->biotype ne $gene->biotype  && !$biotype){
	  print "The non coding gene overlap is of a different type, ".  $gene->biotype . " vs ".$overlap->biotype." not transferring this stable id \n";
	  next GENE;
	}
	if($biotype && $overlap->biotype ne $biotype ){
	  print "The non coding gene overlapping yours is of a different type, so it won't get deleted: $biotype vs ".$overlap->biotype." not transferring this stable id \n";
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
	# ignore overlaps if they are miRNAs we always end up letting these through anyway!
	next GENE if $gene->biotype eq 'miRNA';
	foreach my $trans (@{$overlap->get_all_Transcripts}) {
	  foreach my $exon (@{$trans->get_all_translateable_Exons}) {
	    my $codingexon = ($slice->start+$exon->end() >= $gene->start() &&
			    $slice->start+$exon->start() <= $gene->end());
	    if ($codingexon) {
	      print "gene ".$gene->dbID."\t".$gene->description." overlaps coding exon "
		.$exon->stable_id." in  real_gene ".$overlap->stable_id."\t";
	      if($makewhitelist){
		print WHITE $gene->dbID."\n";
	      }
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
  return (\%noncoding, \%coding, \%merge_set);
}

sub at_content {
  my ($genes) = @_;
    my %blacklist;
  foreach my $key (keys %$genes) {
    my $perc_at;
    my $longest_at;
    my $count;
    my $gene = $genes->{$key};
    # we want to eliminate genes with high 'AT' content, so we have 23 tests based on analysis of data over all genomes..
    # not miRNAs though, just infernal predictions...
    next if $gene->biotype eq 'miRNA';
    my $seq = $gene->get_all_Transcripts->[0]->seq->seq;
    while ($seq =~ /(AT)/g) { $count++ }
    $perc_at = int($count / length($seq) * 200);
    if ($seq =~ /((AT)+)/){
      $longest_at = length($1);
    }
    # tests 
    if ( $perc_at > 40 or $longest_at > 10 ) {
      $blacklist{$gene->dbID} = 1 ;
    }
  }
  return \%blacklist;
}

sub blacklist {
  my ($overhangs,$duplications,$coding_overlaps,$repeats) = @_;
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
  foreach my $key (keys %$repeats){
    $list{$key} = 1 unless $whitelist->{$key};
  }
  return \%list;
}

sub stable_id_mapping {
  my ($non_coding_overlaps,$old_hash,$new_hash,$blacklist,$merge_set) = @_;
  # get the assembly information
  my $last_session = sql('SELECT max(mapping_session_id) from mapping_session',$final_db)->[0];
  my $new_assembly = sql('SELECT meta_value from meta where meta_key = "assembly.default"',$final_db)->[0];
  my ($last_db, $old_release, $old_assembly);
  if($last_session){
    $last_db = sql("SELECT new_db_name from mapping_session where mapping_session_id = $last_session",$final_db)->[0];
    $old_release = sql("SELECT new_release from mapping_session where mapping_session_id = $last_session",$final_db)->[0];
    $old_assembly = sql("SELECT new_assembly from mapping_session where mapping_session_id = $last_session",$final_db)->[0];
  }
  else{
    $last_session = 0;
    $last_db      = $final_dbname;
    $old_release  = $new_release;
    $old_assembly = $new_assembly;
  }
  my $new_session =  $last_session + 1;
  print SIDS "INSERT INTO  mapping_session(mapping_session_id,old_db_name,new_db_name,old_release,new_release,old_assembly,new_assembly,created) ".
    " VALUES($new_session,\'$last_db\',\'$final_dbname\',\'$old_release\',\'$new_release\',\'$old_assembly\',\'$new_assembly\',now());\n";
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
    $done{$old_gene->dbID} = 1;
  }
  # deleting
  foreach my $old_gene (keys %$old_hash){
    next if $done{$old_gene};
    if ( $merge ){
      # we might not delete all the old genes so we dont want stable id mappings for them
      next unless ( $old_hash->{$old_gene}->description =~ /RFAM/ or $old_hash->{$old_gene}->biotype eq 'miRNA' or $merge_set->{$old_gene} )
    }
    my $gene =  $old_hash->{$old_gene};
    my $trans = $gene->get_all_Transcripts->[0];
    # old gene does not need stable id transferring stable id is dead 
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES('".
      $gene->stable_id."',".$gene->version.",null,0,$new_session,'gene',0);\n";
    print SIDS "INSERT INTO stable_id_event(old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES('".
      $trans->stable_id."',".$trans->version.",null,0,$new_session,'transcript',0);\n";
    # need gene archive entries also...
    print SIDS "INSERT INTO gene_archive(gene_stable_id,gene_version,transcript_stable_id,transcript_version,translation_stable_id,translation_version,peptide_archive_id,mapping_session_id) ";
    print SIDS "VALUES('". $gene->stable_id."',".$gene->version.",'". $trans->stable_id."',".$trans->version.",'',0,0,$new_session);\n";
  }
  # new stable ids need to be added with the appropirate code...
  return $new_session;
}

sub generate_new_ids{
  my ($ncRNAs,$ga,$list,$last_session) = @_;
  my ($gsp,$gsi,$tsp,$tsi,$esp,$esi);
  # ensembl_ids all have 11 numbers at the end and an indeterminate number of letters at the start
  # The stable id not like SIN% is to fix a problem with the old stable ids that were retired in fugu
  if (sql("SELECT max(stable_id) from gene_stable_id WHERE left(stable_id,3) not in ('SIN','NEW') ;",$ga)->[0] =~ /^(\D+0+)(\d+)$/){
    $gsp = $1;
    $gsi = $2;
  }
  if ( sql("SELECT max(stable_id) from transcript_stable_id WHERE left(stable_id,3) not in ('SIN','NEW') ;",$ga)->[0] =~ /^(\D+0+)(\d+)$/){
    $tsp = $1;
    $tsi = $2;
  }
  if ( sql("SELECT max(stable_id) from exon_stable_id WHERE left(stable_id,3) not in ('SIN','NEW');",$ga)->[0] =~ /^(\D+0+)(\d+)$/){
    $esp = $1;
    $esi = $2;
  }
  # check stable ids of dead genes are not higher than the maximum in the gene stable id table
  if (sql("SELECT max(gene_stable_id) from gene_archive WHERE left(gene_stable_id,3) not in ('SIN','NEW') ;",$ga)->[0] =~ /^(\D+0+)(\d+)$/){
    if ($2 > $gsi){
      print "dead gene with higher id $gsp$gsi\n";
      $gsp = $1;
      $gsi = $2;
      print " becomes $gsp$gsi\n";
    }
  }
  if ( sql("SELECT max(transcript_stable_id) from gene_archive  WHERE  left(transcript_stable_id,3) not in ('SIN','NEW') ;",$ga)->[0] =~ /^(\D+0+)(\d+)$/){
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
  my ($old_hash,$ga,$merge_set) = @_;
  return if scalar(keys %$old_hash) == 0;
  print STDERR "Warning you have not specified a biotype *ALL* ncRNAs will be deleted\n" unless ($biotype);
  print STDERR "Found ".scalar(keys %$old_hash)." genes  $final_dbname\nshall I delete them? (Y/N) ";
  my $reply = <>;
  chomp $reply;
  if ($reply eq "Y" or $reply eq "y") {
    foreach my $key (keys %$old_hash){
      my $gene;
      if ( $merge ){
	# delete ALL the RFAM  / miRBase genes and any overlapping Sean Eddy genes
	if ( $old_hash->{$key}->description =~ /RFAM/ or $old_hash->{$key}->biotype eq 'miRNA' or $merge_set->{$key} ){
	  print "MERGE: " . $old_hash->{$key}->description . " ";

	  print $old_hash->{$key}->seq_region_name . " " . 
	    $old_hash->{$key}->biotype . " " .
	      $old_hash->{$key}->seq_region_start . " " .
		$old_hash->{$key}->seq_region_end . " " .
		  $old_hash->{$key}->seq_region_strand . "\n";
	  $gene = lazy_load($old_hash->{$key});
	}
      } else {
	$gene = lazy_load($old_hash->{$key});
      }
      next unless (defined($gene));
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

    #set dates
    $gene->created_date(time());
    $gene->modified_date(time());
    $trans->created_date(time());
    $trans->modified_date(time());
    foreach my $exon (@{$trans->get_all_Exons}){
      $exon->created_date(time());
      $exon->modified_date(time());
    }

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
	  next unless ($xref->dbname eq 'miRBase' or $xref->dbname eq 'RFAM' or $xref->dbname eq 'HGNC');
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

sub delete_sql {
  my ($query,$db) = @_;
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  return;
}

sub fetch_genes {
  my ($ga,$biotype,$biotype_to_skip,$slice,$genestoignore,$sa) = @_;
  my %ncRNA_hash;
  my @ncRNAs;
  my %ignored_ncRNAs;
  throw("Cannot fetch genes without gene adaptor $ga") unless $ga;
  if ($slice){
    my @slices = @{$sa->fetch_all('toplevel')};
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
      } elsif ($biotype_to_skip) {
	@ncRNAs =  @{$ga->fetch_all_by_Slice_constraint($slice,"biotype != '".$biotype_to_skip."'")};
      } else {
	@ncRNAs =  @{$ga->fetch_all_by_Slice_constraint($slice,'biotype like "%RNA%" ')};
      }
      foreach my $ncRNA (@ncRNAs) {
	next unless ($ncRNA->analysis->logic_name eq 'ncRNA' or $ncRNA->analysis->logic_name eq 'miRNA') ;
	unless ($ncRNA->biotype eq 'miRNA' or
		$ncRNA->biotype eq 'misc_RNA' or
		$ncRNA->biotype eq 'snRNA' or
		$ncRNA->biotype eq 'snoRNA' or
		$ncRNA->biotype eq 'rRNA') {
	  $ignored_ncRNAs{$ncRNA->biotype} ++;
	  next;
	}
	next if  $ncRNA->description =~ /RNAI/;
	#skip selected genes
	foreach my $genestoignore (@$genestoignore){
	  next if  $ncRNA->description =~ /$genestoignore/;
	}
	$ncRNA_hash{$ncRNA->dbID} = lazy_load($ncRNA);
      }
    }
  } else {
      if ($biotype) {
	@ncRNAs =  @{$ga->generic_fetch("biotype = '".$biotype."'")};
      } elsif ($biotype_to_skip) {
	@ncRNAs =  @{$ga->generic_fetch("biotype != '".$biotype_to_skip."'")};
      } else {
	@ncRNAs =  @{$ga->generic_fetch('biotype like "%RNA%"')};
      }      
      foreach my $ncRNA (@ncRNAs) {
	next unless ($ncRNA->analysis->logic_name eq 'ncRNA' or $ncRNA->analysis->logic_name eq 'miRNA') ;
	unless ($ncRNA->biotype eq 'miRNA' or
		$ncRNA->biotype eq 'misc_RNA' or
		$ncRNA->biotype eq 'snRNA' or
		$ncRNA->biotype eq 'snoRNA' or
		$ncRNA->biotype eq 'rRNA') {
	  $ignored_ncRNAs{$ncRNA->biotype} ++;
	  next;
	}
	next if  $ncRNA->description =~ /RNAI/;
	foreach my $genestoignore (@$genestoignore){
	  next if  $ncRNA->description =~ /$genestoignore/;
	}
	$ncRNA_hash{$ncRNA->dbID} = lazy_load($ncRNA);
      }
    }
  print STDERR  "\n";
  print "Ignoring the following ncRNAs:\n";
  foreach my $key ( keys %ignored_ncRNAs ) {
    print "\t$key\t" . $ignored_ncRNAs{$key} . "\n";
  }
  return \%ncRNA_hash;
}

sub check_meta {
  my ($db1,$db2) = @_;
  my $m1 = $db1->get_MetaContainer;
  my $m2 = $db2->get_MetaContainer;
  my $c1 = sql("SELECT meta_value from meta where meta_key = 'assembly.default'",$db1)->[0];
  my $c2 = sql("SELECT meta_value from meta where meta_key = 'assembly.default'",$db2)->[0];

  unless ($m1->get_taxonomy_id eq $m2->get_taxonomy_id ){
    throw("Tax ids do not agree ".$m1->get_taxonomy_id." != ". $m2->get_taxonomy_id."\n");
  }
  unless ($c1 eq $c2){
    throw("Coord systems do not agree $c1 != $c2 \n");
  }
  return;
}


__END__
