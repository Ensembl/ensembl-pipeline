#!/usr/local/ensembl/bin/perl
# Author: Emmanuel Mongin
# Creation: 03.19.2001


=head1 Run_protein_RunnableDB

=head2 Description

This script will submit run_protein_RunnableDB
    This is the script which is used to run the full protein analysis. The protein analysis is specific in a way that\'s the algorithms used have different speed. For example the difference of speed between hmmpfam and low complexity is huge. Moreover the analysis is done at the peptide level (much less material than for a whole).
    There is then 3 ways of running an annotation:
    + The classical way, protein by protein (used for hmmpfam for example)
    + The intermediate way, using chunks of the main peptide file. For example using chunks containing 100 peptides in fasta format (used for PRINTS for example).
    + All against all. Uses a file containing all the peptide (used for low complexity for example)

    For information read:
    ensembl-doc/protein-annotation.txt

=head2 Usage

  Simply: perl submit_prot_analysis.pl
  But the configuration file Bio/EnsEMBL/Pipeline/Prot_analysis_Conf.pl has to be properly filled in.

=cut

use strict;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/Prot_analysis_Conf.pl";
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;

my %db_conf =  %::db_conf;
my %scripts_conf = %::scripts_conf;

my %runnables;

#Get the DB variables
my $dbname = $::db_conf{'dbname'};
my $dbhost = $::db_conf{'dbhost'};
my $dbuser = $::db_conf{'dbuser'};
my $dbpass = $::db_conf{'dbpass'};

#Get the DB handler

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $::db_conf{'dbhost'},
					    -user   => $::db_conf{'dbuser'},
					    -dbname => $::db_conf{'dbname'},
					    -pass => $::db_conf{'dbpass'},
					    );

#Get the location of the peptide file
my $pep_file = $::scripts_conf{'pep_file'};

#Get the location of the scratch directory
my $scratchdir =  $::scripts_conf{'tmpdir'};

print STDERR "Using the following scratch directory: ".$scratchdir."\n";

#Give for each analysis its corresponding runnable

#prints, prosite, profile, pfam, scanprosite, tmhmm, coils, signalp, seg

$runnables{'Prints'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints";
$runnables{'Prosite'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ScanProsite";
$runnables{'Pfam'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Hmmpfam   ";
$runnables{'Tmhmm'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Tmhmm";
$runnables{'ncoils'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Coil";
$runnables{'Signalp'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Signalp";
$runnables{'Seg'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Seg";
$runnables{'Profile'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Profile";
$runnables{'ParacelHMM'} = "Bio::EnsEMBL::Pipeline::RunnableDB::Protein::ParacelHMM";

#Get all of the protein Ids for the whole genome
my @ids = &get_ids();

&make_directories();
&chunk_pepfile();

&run_jobs(@ids);

sub get_ids {
    my @id;
    
#Get all of the protein ids for the given database
    
    my $q = "select translation_id from transcript";
    
    my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
    my $res = $sth->execute || $db->throw("can't execute: $q");
    
    while( my ($prot_id) = $sth->fetchrow_array) {
	push(@id,$prot_id);
    }
    return @id;
}


sub make_directories {
#Find out which analysis should be run
    my $run = $::scripts_conf{'2berun'};
    #print STDERR "making directories\n";
    my @toberun = split /,/, $run;

    my $chunks = $scratchdir."/chunks";
    &makedir($chunks);

#First make the tmp directories for each analysis which is supposed to run
    foreach my $r(@toberun) {
	&makedir($scratchdir."/".$r);

	my $dir_err = $scratchdir."/$r/stderr";
	my $dir_out = $scratchdir."/$r/stdout";
	&makedir($dir_err);
	&makedir($dir_out);
    }
}

sub chunk_pepfile {
#Chunk the peptide file
    open (PEPFILE, "$pep_file");
    my $count = 0;
    my $chunk = 1;
    #print STDERR "chunking peptide file\n";
    my $size = $::scripts_conf{'chunk_size'};

    $/ = "\>";
    #print "have opened ".$pep_file."\n";
    while(<PEPFILE>){
	#print $_."\n";
	if ($_ ne "\>") {
	    if ($count == 0) {
		open (CHUNK,">".$scratchdir."/chunks/chunk.$chunk");
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

sub run_jobs {
    my (@id) = @_;
    
    #Find out which analysis should be run
    my $run = $::scripts_conf{'2berun'};
    #print STDERR "running jobs\n";
    my @toberun = split /,/, $run;

    my $queue = $::scripts_conf{'queue'};
    my $runner = $::scripts_conf{'runner'};
    
    my $pep_file = $::scripts_conf{'pep_file'};
    
     foreach my $r(@toberun) {
	 #First get the analysisId of the module which is supposed to run
	 #print STDERR "running ".$r."\n";
	 my $q = "select analysis_id from analysis where module = '$r'";
    
	 my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
	 my $res = $sth->execute || $db->throw("can't execute: $q");
    
	 my ($analysis_id) = $sth->fetchrow_array;

	 unless ($analysis_id) {die "AnalysisId not defined for ".$r."\n"};

	 #Try to know how this analysis should be run
	 my $chunk_name = $r."_chunk";

	 my $chunk = $::scripts_conf{$chunk_name};

	 unless ($chunk) {die "No chunk option defined for ".$r."\n"};

	 my $runnabledb = $runnables{$r};

	 unless ($runnabledb) {die "No runnableDB defined for ".$r."\n"};

	 my $check = "-E \"".$runner." -check -runnable ".$runnabledb."\"";
	 #print STDERR "using chunk ".$chunk."\n";
#Three diffenrent way of running the analysis 1,2,3 (see documentation in the configuration file and in the perldoc of this script)
	 if ($chunk == 1) {
	     my @chunk;
	     my $count = 0;
	     foreach my $i(@id) {
		 $count++;
		 push (@chunk,$i);
		 if ($count >= 25) {
		     my $subm = join(":",@chunk);
		     
		     print STDERR "ID: $subm\n";
		     my $command = "bsub -Ralpha -C 0 -J ".$r." -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$i.".out -e ".$scratchdir."/".$r."/stderr/".$i.".err ".$check." ".$runner." -runnable ".$runnabledb." -dbuser ".$dbuser." -dbpass ".$dbpass." -dbname ".$dbname." -host ".$dbhost." -input_id ".$subm." -analysis ".$analysis_id;
				 
		     print STDERR "RUNNING: $command\n";
		 
		     system($command)==0 || die "$0\Error running '$command' : $!";
		     @chunk = undef;
		     $count = 0;
		 }
	     }
	 }

	 
	 if ($chunk == 3) {
	     my $command = "bsub -Ralpha -C 0 -J ".$r." -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$r.".out -e ".$scratchdir."/".$r."/stderr/".$r.".err ".$check." ".$runner." -runnable ".$runnabledb." -dbuser ".$dbuser." -dbpass ".$dbpass." -dbname ".$dbname." -host ".$dbhost." -input_id ".$pep_file." -analysis ".$analysis_id;
	    

	     print STDERR "RUNNING: $command\n";
	     
	     system($command)==0 || die "$0\Error running '$command' : $!";
	     
	    
	 }
	 
	 if ($chunk == 2) {
	   print STDERR "have chunk ".$chunk."\n";
	     my $dir = $scratchdir."/chunks";     
	     opendir(DIR, $dir);
	     my @allfiles = readdir DIR;
	     closedir DIR;
	     
	     foreach my $f(@allfiles) {
	       print STDERR "have file ".$f."\n";
		 if (($f ne ".") && ($f ne "..")) {
		      my $command = "bsub -Ralpha -C 0 -J ".$r." -q ".$queue." -o ".$scratchdir."/".$r."/stdout/".$f.".out -e ".$scratchdir."/".$r."/stderr/".$f.".err ".$check." ".$runner." -runnable ".$runnabledb." -dbuser ".$dbuser." -dbpass ".$dbpass." -dbname ".$dbname." -host ".$dbhost." -input_id ".$dir."/".$f." -analysis ".$analysis_id;
		  

		      print STDERR "RUNNING: $command\n";
	     
		      system($command)==0 || die "$0\Error running '$command' : $!";

		  }
	     }
	 }
		     
     }
}


sub makedir{
  my ($dir) = @_;
  if(opendir(DIR, $dir)){ closedir(DIR); }
  else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}
