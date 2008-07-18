#!/usr/local/ensembl/bin/perl -w 

#script to gather evidence from the various analyses I have run to determine why cDNAs have no alignments in the db

#perl identify_hit_failings.pl 

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::KillList::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases;


use Getopt::Long;

my $record_separator = $/;

my ($gss, $seq_file, $user, $host);
my ($port, $dbname, $species, $vertrna, $refseq);
my ($vertrna_update, $infile, $outfile, $findN_prog);
my ($reasons_file);

&GetOptions(
            'gss=s'      	=> \$gss,
            'seq_file=s'      	=> \$seq_file,
            'user=s'            => \$user,
            'host=s'    	=> \$host,
            'port=s'    	=> \$port,
            'dbname=s'    	=> \$dbname,
            'species=s'    	=> \$species,
            'vertrna=s'    	=> \$vertrna,
            'refseq=s'    	=> \$refseq,
            'vertrna_update=s' 	=> \$vertrna_update,
            'infile=s' 		=> \$infile,
            'outfile=s' 	=> \$outfile,
            'findN_prog=s' 	=> \$findN_prog,
       	    'reasons_file=s'    => \$reasons_file
	   );

my $db1 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host    => $host,
						     -port    => $port,
						     -user    => $user,
						     -dbname  => $dbname
						    );


my (%cdnas);
my $cov_threshold = 90;
my $pid_threshold = 97;

open (OUT, ">", "$outfile") or die("can t open $outfile $!\n");
open(IN, "<", "$infile") or die("can t read $infile $!\n");

while(<IN>){ 
	chomp;
	my @tmp = split /\s+/, $_;
	if ($tmp[0]=~/rpp/){
		#rejected because thought to be a processed pseudogene
		$cdnas{$tmp[1]}{"Processed pseudogene"} = 1;

	}elsif($tmp[0]=~/only/){
		#rejected because only contains introns > 250000bp
		$cdnas{$tmp[3]}{"All long introns"} = 1;

	}elsif($tmp[0]=~/reject/){
		#rejected because contains an intron 250000-400000bp long
		#but coverage and pid are not both > 98%
		#row looks like: reject: intron 350849 coverage 99.09 %id 94.80 AY182732.1
		if ($tmp[4] < 98){
			$cdnas{$tmp[7]}{"Low coverage with long intron"} = $tmp[4];
		}elsif ($tmp[6] < 98){
			#store the best match - ie highest pid - sorted in ExonerateTranscriptFilter so just want first hit
			if (!exists $cdnas{$tmp[7]}{"Low percent_id with long intron"}){
				$cdnas{$tmp[7]}{"Low percent_id with long intron"} = $tmp[6];
			}
		}

	}elsif($tmp[0]=~/max_coverage/ && scalar @tmp == 7){
		#rejected because coverage or %id did not reach the thresholds, 
		#current threshold = 90% coverage and 97% id
		#or there was a match with better coverage

		#row looks like: max_coverage 66.18 coverage 66.18 %id 98.61 AK134284.1

		#store the max coverage
		my $max_cov = $tmp[1];

		if ($tmp[3] < $tmp[1]){
			$cdnas{$tmp[6]}{"other_hit_with_higher_cov"} = $tmp[3];	
		}
		elsif ($tmp[3] < 90){ #so best coverage too low:
			$cdnas{$tmp[6]}{"Low coverage"} = $tmp[3];	
		}
		elsif ($tmp[5] < 97){ #do we only want to record highest pid hit?
			if ((!defined $cdnas{$tmp[6]}{"low percent_id"}) || ($tmp[5] > $cdnas{$tmp[6]}{"Low percent_id"})){
				$cdnas{$tmp[6]}{"Low percent_id"} = $tmp[5];
			}	
		}

	}else{
		print "strange line $_\n";
	}
}
close IN;


#get entries for species of interest, combine base file & update file
#read update file

my @files = ($vertrna_update, $vertrna);
my (%EMBL_ids);


foreach my $file(@files){
	open(IN, "<", $file) or die("can t read $file\n");
	while (my $entry = <IN>){
		if($entry =~ m/$species/){
			#extract & save id
			$entry =~ /^>[\w\d]+\s([\w\.\d]+)\s.+/;
			if(!$1){ die "\n$file: unmatched id pattern:\n$entry\n"; }
			$EMBL_ids{$1} = 1;
		}
	}
	close IN;
}


#read RefSeq file - has diff format so have to do separately
open(IN, "<", $refseq) or die("can t read $refseq.\n");

while (my $entry = <IN>){
	my $id = "";
	if($entry =~ m/^>gi.+ref\|(NM_.+)\| $species.*/){
		$id = $1;
	}
	elsif($entry =~ m/^>gi.+ref\|(NR_.+)\| $species.*/){
		$id = $1;
	}
	else{
		next;
	}
	if($id){
	  #reduce header to accession number
	  $EMBL_ids{$id} = 1;
	}
}
close IN;   

#list those which did not give an output from ExonerateTranscriptFilter:
#this may be because they are in the database...!
foreach my $k (keys %EMBL_ids){
	if (!exists $cdnas{$k}){	
		$cdnas{$k}{"No output from Exonerate"} = 1;
	}
}

#now run the findN.pl to see if have large strings of Ns:
`perl $findN_prog $seq_file > many_n.out`;

open IN, "many_n.out";
while (<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	$cdnas{$tmp[0]}{">10% N-strings"} = 1;
}
close IN;

#check the database to see if any of thse cdnas has a hit in the sd3_cDNA_update_tmp db
#might have if aligned after polyA processing or final run with rpp rule + maxintron + no softmask...

#want to find the list of hit_names with a match in the db
my $sql = ("select distinct hit_name from dna_align_feature"); 

my $q1 = $db1->dbc->prepare($sql) or die "sql error";
$q1->execute();

my $n = keys %cdnas;
#print "had $n cdnas...";
while (my $hit = $q1->fetchrow_array){
	#if exists in cdnas - remove because it has an alignment in the db
	if ($cdnas{$hit}){
		#print "$hit\n";
		delete $cdnas{$hit}; 
	}
}	

#now have done other analyses need to check whether the reason we have no output is because
#cdnas are on the kill list

#re-create the kill_list hash, from the file dumped by new_cDNA_upadte.pl 
my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => 'cDNA_update');
my %kill_list = %{$kill_list_object->get_kill_list()};

foreach my $k (keys %kill_list){
  foreach my $cdna_acc (keys %cdnas) {
    my $truncated_acc;
    if ($cdna_acc =~ /^(\w+)\./) {
      $truncated_acc = $1;
    }
    if ($k eq $truncated_acc) {
      delete $cdnas{$cdna_acc};
      $cdnas{$cdna_acc}{"See kill list"} = 1;
    }
  }
}


open(LIST, "<", $gss) or die("can't open gss list $gss");
my %gss_acc_list;
while (<LIST>){
	my @tmp = split/\s+/, $_;
	$gss_acc_list{$tmp[1]} = 1;
}
close LIST;

my %gss_hits;
foreach my $k (keys %cdnas){ #because looping over cdnas - shouldn't delete as you go along - have to do after
	if ($k=~/(\w+)\./){
		my $acc = $1;
		if (exists $gss_acc_list{$acc}){
			$gss_hits{$k} = 1;
		}
	}
}		
#now loop through ones on gss_list:
foreach my $sv (keys %gss_hits){
	if (exists $cdnas{$sv}){
		delete $cdnas{$sv}; #gss taking precedence over kill list	
		$cdnas{$sv}{"GSS sequence"} = 1;
	}
}	 



#print out details about these missing sequences - possible reasons they are missing:

#sort the hit_names to make the output file easier to browse:
my @sorted = sort keys %cdnas;

my $type = 'cDNA';
my $target_score = 'NULL';
my $ensembl_id = 'NULL';
my $ensembl_object_type = 'NULL';
my $analysis_id;
my $external_db_id;
my $id;
my $reason_id;
my $query_score;


#set the analysis_id:
$sql = ("select analysis_id from analysis where logic_name = 'cDNA_update'"); 
$q1 = $db1->dbc->prepare($sql) or die "sql error";
$q1->execute();
$analysis_id = $q1->fetchrow_array;
print "unmapped_analysis_id = $analysis_id\n";

#set the refseq and embl_ids
$sql = ("select external_db_id from external_db where db_name = 'RefSeq_dna'"); 
$q1 = $db1->dbc->prepare($sql) or die "sql error";
$q1->execute();
my $refseq_id = $q1->fetchrow_array;

$sql = ("select external_db_id from external_db where db_name = 'EMBL'"); 
$q1 = $db1->dbc->prepare($sql) or die "sql error";
$q1->execute();
my $embl_id = $q1->fetchrow_array;
print "RefSeq = $refseq_id and EMBL = $embl_id\n";

#read reasons file into database
my %reason_list;
open (IN, "$reasons_file") or die("Can't open $reasons_file $!\n");
while(<IN>){
	my @tmp = split/\t/, $_;
	$reason_list{$tmp[1]} = $tmp[0];
}
close IN;


for my $id (@sorted){
	#only want to print out ones with maximum coverage because using 'best in genome'
	
	my %reasons = (); 
	
	#decide whether cDNA is from is embl or refseq
	if ($id=~/^N[MR]/){
		$external_db_id = $refseq_id;
	}else{
		$external_db_id = $embl_id;
	}
	
	if (exists $cdnas{$id}{"GSS sequence"}){
		$reason_id = $reason_list{"GSS sequence"};
		$query_score = 'NULL';
		$reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{"See kill list"}){
		$reason_id = $reason_list{"See kill list"};
		$query_score = 'NULL';
		$reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{"Low coverage"}){
		$reason_id = $reason_list{"Low coverage"};
		$query_score = $cdnas{$id}{"Low coverage"};
		$reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{"Low percent_id"}){
		$reason_id = $reason_list{"Low percent_id"};
		$query_score = $cdnas{$id}{"Low percent_id"};
		$reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{"Processed pseudogene"}){
		$reason_id = $reason_list{"Processed pseudogene"};
		$query_score = 'NULL';
		$reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{"All long introns"}){
		$reason_id = $reason_list{"All long introns"};
		$query_score = 'NULL';
		$reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{"Low coverage with long intron"}){
		$reason_id = $reason_list{"Low coverage with long intron"};
		$query_score = $cdnas{$id}{"Low coverage with long intron"};
                $reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{"Low percent_id with long intron"}){
		$reason_id = $reason_list{"Low percent_id with long intron"};
		$query_score = $cdnas{$id}{"Low percent_id with long intron"};
		$reasons{$reason_id} = $query_score;
	}
	if (exists $cdnas{$id}{">10% N-strings"}){
		$reason_id = $reason_list{">10% N-strings"};
		$query_score = 'NULL';
		$reasons{$reason_id} = $query_score;

	}
	if (exists $cdnas{$id}{"No output from Exonerate"}){
		$reason_id = $reason_list{"No output from Exonerate"};
		$query_score = 'NULL';
		$reasons{$reason_id} = $query_score;
	}
	
	
	foreach my $reason_id (keys %reasons){
		print OUT "INSERT INTO unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id, query_score,".
		          "target_score, ensembl_id, ensembl_object_type) VALUES ('$type', $analysis_id, $external_db_id, '$id', $reason_id, ".
				  "$reasons{$reason_id}, $target_score, $ensembl_id, $ensembl_object_type);\n";
    }
}

close OUT;



		

