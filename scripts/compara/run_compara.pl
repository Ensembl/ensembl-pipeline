#!/usr/local/ensembl/bin/perl -w

# EXIT STATUS
# 0 all is fine
# 1 Problem in parsing input file

use strict;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer;
use Getopt::Long;

$| = 1;

my $usage = "
Usage: $0 [options] input_file

Line input in input_file should be 
sequence_source:species:dnafrag_type:dnafrag_name::sequence_source:species:dnafrag_type:dnafrag_name,

e.g.

ENSEMBL:Homo_sapiens:VirtualContig:1.106750001.107000000::ENSEMBL:Mus_musculus:VirtualContig:3.111750001.112000000

Get from a Ensembl core database a human DNA sequence, get it as a VirtualContig,
on chromosome 1 between position 106750001 and 107000000 (this sequence will be the reference sequence)
Get from a Ensembl core database a mouse DNA sequence, get it as a VirtualContig,
on chromosome 3 between position 111750001 and 112000000
And compare the human sequence to the mouse sequence

or

ENSEMBL:Homo_sapiens:VirtualContig:20.1.1000000::FASTA:Mus_musculus:VirtualContig:/path/to/fasta/file/mouse_chr2.fa

Get from a Ensembl core database a human DNA sequence, get it as a VirtualContig,
on chromosome 20 between position 1 and 1000000 (this sequence will be the reference sequence)
Get from a Fasta file mouse DNA sequences. Sequences are named as a VirtualContig,
chr_name.chr_start.chr_end
And compare the human sequemce to the mouse sequences present in the Fasta file.

 -help       show this menu
 -h          compara database hostname
 -d          compara database name
 -u          database user name
 -p          database password
 -alnprog    program used for alignment (default: bl2seq) see perldoc CrossComparer.pm for details
 -alntype    program used for alignment (default: blastn) see perldoc CrossComparer.pm for details
 -parameters parameters (options) added to the executable command line
             e.g. \"-g T -W 10 -G 1 -E 2\" (default
 -min_score  minimum score to take in account in results (default: 40)
 -masked     0 unmasked; 1 masked; 2 soft masked (default: 0)
 -filter     specify a perl module that will do hit filtering before writing them in compara database
             e.g. \"Bio::EnsEMBL::Compara::Filter::Greedy\" (default)
";


my $help = 0;
my ($host,$dbname,$dbuser,$pass);

if (scalar @ARGV == 0) {
  print $usage;
  exit 0;
}


my $alnprog = 'bl2seq';
my $alntype = '';
my $min_score = 40;
my $masked = 1;
my $filter = "Bio::EnsEMBL::Compara::Filter::Greedy";
my $parameters = "-g T -W 10 -G 1 -E 2";

unless (GetOptions('help' => \$help,
		   'h=s' => \$host,
		   'd=s' => \$dbname,
		   'u=s' => \$dbuser,
		   'p=s' => \$pass,
		   'alnprog=s' => \$alnprog,
		   'alntype=s' => \$alntype,
		   'parameters=s' => \$parameters,
		   'min_score=f' => \$min_score,
		   'masked=i' => \$masked,
		   'filter=s' => \$filter)) {
  die $usage;
}

if ($help) {
  print $usage;
  exit 0;
}

# Some checks before starting

unless (scalar @ARGV) {
  die $usage;
}

if ($alnprog eq "bl2seq" && $alntype eq "") {
  $alntype = "blastn";
}

print STDERR "
Unix start time           : ",time,"
Compara database used     : $dbname @ $host
Alignment program used    : alnprog \"$alnprog\"; alntype \"$alntype\"; parameters \"$parameters\";
Filter perl module used   : \"$filter\"
Type of masking           : $masked
";

my $input_file = $ARGV[0];

if (! defined $input_file) {
  die "Must call run_compara with an input file!";
}

# Connecting to compara database

my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor (-host => $host,
						      -user => $dbuser,
						      -pass => $pass,
						      -dbname => $dbname );

# Openning input file

open INPUT, "$input_file" || die "could not open $input_file, $!\n";

# Kepping traces of what was already run

my $log_file = $input_file.".log";
my %already_run;

if (-e $log_file) {
  open LOG, "$log_file";
  while (defined (my $line =<LOG>)) {
    chomp $line;
    $already_run{$line} = 1;
  }
  close LOG;
}

open LOG, ">> $log_file";

my $exit_jobs_file = "exit_jobs_file";

while (defined (my $line = <INPUT>)) {

# Clean exit if $exit_jobs_file exists :)

  if (-e $exit_jobs_file) {
    print STDERR "Job properly ended before all pairwise comparisons were completed\n";
    print STDERR "endtime: ",time,"\n";
    exit 0;
  }

  if ($line =~ /^\S+:\S+:\S+:\S+::\S+:\S+:\S+:\S+$/) {
    chomp $line;
  } else {
    warn "Line input in $input_file should be sequence_source:species:dnafrag_type:dnafrag_name::sequence_source:species:dnafrag_type:dnafrag_name. See -help.
exit 1";
    exit 1;
  }

  next if (defined $already_run{$line});
  
  # RunnableDB need an analysis to run properly
  my $analysis = new Bio::EnsEMBL::Analysis(-parameters => $parameters); 

  my $rundb = new Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer (-alnprog => $alnprog,
								     -alntype => $alntype,
								     -db => $db,
								     -input_id => $line,
								     -min_score => $min_score,
								     -masked => $masked,
								     -filter => $filter,
								     -analysis => $analysis);

  # run the runnabledb
  
  print STDERR "Input comparison          : $line\n";
  
  $rundb->fetch_input();
  $rundb->run();
  if ($rundb->write_output()) {
    $already_run{$line} = 1;
    print LOG "$line\n";
  } else {
    print STDERR "Error in write_output\n";
  }
}

close LOG;
close INPUT;

print STDERR "Unix end time             : ",time,"\n";
