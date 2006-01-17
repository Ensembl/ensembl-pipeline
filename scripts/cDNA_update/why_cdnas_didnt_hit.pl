#!/usr/local/ensembl/bin/perl -w 

#script to gather evidence from the various analyses I have run to determine why cDNAs have no alignments in the db

#perl identify_hit_failings.pl 

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;



my ($kill_list, $gss, $seq_file, $user, $host);
my ($port, $dbname, $species, $vertrna, $refseq);
my ($vertrna_update, $infile, $outfile, $findN_prog);

&GetOptions(
            'kill_list=s'       => \$kill_list,
            'gss=s'      		=> \$gss,
            'seq_file=s'      	=> \$seq_file,
            'user=s'          	=> \$user,
            'host=s'    		=> \$host,
            'port=s'    		=> \$port,
            'dbname=s'    		=> \$dbname,
            'species=s'    		=> \$species,
            'vertrna=s'    		=> \$vertrna,
            'refseq=s'    		=> \$refseq,
            'vertrna_update=s' 	=> \$vertrna_update,
            'infile=s' 			=> \$infile,
            'outfile=s' 		=> \$outfile,
            'findN_prog=s' 		=> \$findN_prog
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
		$cdnas{$tmp[1]}{"some_look_like_pseudos"} = 1;

	}elsif($tmp[0]=~/only/){
		#rejected because only contains introns > 250000bp
		$cdnas{$tmp[3]}{"all_long_introns"} = 1;

	}elsif($tmp[0]=~/reject/){
		#rejected because contains an intron 250000-400000bp long
		#but coverage and pid are not both > 98%
		#row looks like: reject: intron 350849 coverage 99.09 %id 94.80 AY182732.1
		if ($tmp[4] < 98){
			if (!exists $cdnas{$tmp[7]}{"intron_length"}){
				$cdnas{$tmp[7]}{"long_intron_low_cov"} = $tmp[4];
				$cdnas{$tmp[7]}{"intron_length"} = $tmp[2];
				#print $tmp[7]." ".$cdnas{$tmp[7]}{"long_intron_low_cov"}." ".$cdnas{$tmp[7]}{"intron_length"}."\n";
			}
		}elsif ($tmp[6] < 98){
			#store the best match - ie highest pid - sorted in ExonerateTranscriptFilter so just want first hit
			if (!exists $cdnas{$tmp[7]}{"intron_length"}){
				$cdnas{$tmp[7]}{"long_intron_low_pid"} = $tmp[6];
				$cdnas{$tmp[7]}{"intron_length"} = $tmp[2];
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
			$cdnas{$tmp[6]}{"cov_is_too_low"} = $tmp[3];	
		}
		elsif ($tmp[5] < 97){ #do we only want to record highest pid hit?
			if ((!defined $cdnas{$tmp[6]}{"pid_is_too_low"}) || ($tmp[5] > $cdnas{$tmp[6]}{"pid_is_too_low"})){
				$cdnas{$tmp[6]}{"pid_is_too_low"} = $tmp[5];
			}	
		}

	}else{
		print "strange line $_\n";
	}
}



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
	close(IN);
}


#read RefSeq file - has diff format so have to do separately
open(IN, "<", $refseq) or die("can t read $refseq.\n");

while (my $entry = <IN>){
	my $header = "";
	if($entry =~ m/^>gi.+ref\|(NM_.+)\| $species.*/){
		$header = $1;
	}
	elsif($entry =~ m/^>gi.+ref\|(NR_.+)\| $species.*/){
		$header = $1;
	}
	else{
		next;
	}
	if($header){
	  #reduce header to accession number
	  $EMBL_ids{$header} = 1;
	}
}
close IN;   

#list those which did not give an output from ExonerateTranscriptFilter:
#this may be because they are in the database...!
foreach my $k (keys %EMBL_ids){
	if (!exists $cdnas{$k}){	
		$cdnas{$k}{"no_output"} = 1;
	}
}

#now run the findN.pl to see if have large strings of Ns:
`perl $findN_prog $seq_file > many_n.out`;

open IN, "many_n.out";
while (<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	$cdnas{$tmp[0]}{"N"} = 1;
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
	#if exists in cdnas - remove because it has an alignement in the db
	if ($cdnas{$hit}){
		#print "$hit\n";
		delete $cdnas{$hit}; 
	}
}	

#now have done other analyses need to check whether the reason we have no output is because
#cdnas are on the kill list

open(LIST, "<", $kill_list) or die("can't open kill list $kill_list");
my %kill_list;
while (<LIST>){
	my @tmp = split/\s+/, $_;
	$kill_list{$tmp[0]} = 1;
}
close LIST;
foreach my $k (keys %kill_list){
	if (exists $cdnas{$k}){
		delete $cdnas{$k};	
		$cdnas{$k}{"kill_list"} = 1;
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
		$cdnas{$sv}{"gss_list"} = 1;
	}
}	 



#print out details about these missing sequences - possible reasons they are missing:

#sort the hit_names to make the output file easier to browse:
my @sorted = sort keys %cdnas;

for my $k (@sorted){
	#don't use "other_hit_with_higher_cov"
	#only want to print out ones with maximum coverage because using 'best in genome'
	print OUT "$k ";
	if (exists $cdnas{$k}{"gss_list"}){
		print OUT ": gss sequence ";
	}
	if (exists $cdnas{$k}{"kill_list"}){
		print OUT ": see kill list ";
	}
	if (exists $cdnas{$k}{"cov_is_too_low"}){
		print OUT ": low coverage ".$cdnas{$k}{"cov_is_too_low"}." ";
	}
	if (exists $cdnas{$k}{"pid_is_too_low"}){
		print OUT ": low percent_id ".$cdnas{$k}{"pid_is_too_low"}." ";
	}
	if (exists $cdnas{$k}{"some_look_like_pseudos"}){
		print OUT ": looks like a processed pseudogene ";
	}
	if (exists $cdnas{$k}{"all_long_introns"}){
		print OUT ": all introns >250,000bp ";
	}
	
	if (exists $cdnas{$k}{"long_intron_low_cov"}){
		print OUT ": low coverage ".$cdnas{$k}{"long_intron_low_cov"}." with long intron ".$cdnas{$k}{"intron_length"}." ";
	}
	elsif (exists $cdnas{$k}{"long_intron_low_pid"}){
		print OUT ": low percent_id ".$cdnas{$k}{"long_intron_low_pid"}." with long intron ".$cdnas{$k}{"intron_length"}." ";
	}
	
	if (exists $cdnas{$k}{"N"}){
		print OUT ": >10% N-strings ";
	}
	
	if (exists $cdnas{$k}{"no_output"}){
		print OUT ": No output from Exonerate 0.9.0 ";
	}
	print OUT "\n";
}

close OUT;



		

