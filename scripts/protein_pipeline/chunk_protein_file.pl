#!/usr/local/ensembl/bin/perl -w


use strict;
use Getopt::Long 'GetOptions';
use DBI;
use Bio::EnsEMBL::Pipeline::Config::Protein_Annotation::General qw(
								   PA_PEPTIDE_FILE
								   PA_CHUNKS_DIR
								   PA_CHUNK_SIZE
								  );



# check the last chunk file in a database: this is essential for 
# a db already has protein pipeline analysis. Id, the next chunk file number
# will be incremented which is used in the input_id_analysis table

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname, $help);

GetOptions('dbhost=s' => \$dbhost,
           'dbuser=s' => \$dbuser,
           'dbpass=s' => \$dbpass,
           'dbport=s' => \$dbport,
           'dbname=s' => \$dbname,
          );    # plus default options
				
exec('perldoc', $0) if !($dbhost && $dbuser && $dbpass && $dbport && $dbname);

my $dbh = DBI->connect("DBI:mysql:$dbname:$dbhost:$dbport", $dbuser, $dbpass, {RaiseError => 1, AutoCommit => 0})
        || die "cannot connect to $dbname, $DBI::errstr";

$dbh->debug();

my $last_chunk = get_last_chunk_file();

if ( $last_chunk =~ /\d+/ ){
  warn "\nLast chunk file in input_id_analysis table is $last_chunk: starting from ", $last_chunk+1, "\n";
}
else {
  warn "\nNo last chunk file in input_id_analysis table: starting from 1\n";
  $last_chunk=0;
}

&chunk_pepfile($PA_PEPTIDE_FILE, $PA_CHUNKS_DIR, $PA_CHUNK_SIZE, $last_chunk);

sub get_last_chunk_file {

  my $sql= $dbh->prepare(qq{ SELECT max(input_id) FROM input_id_analysis WHERE input_id_type='filename' });

  $sql->execute;

  my $file = $sql->fetchrow;

  if ( $file ) {
	$file =~ /(\d+)/;
	return $1;
  }
  else {
	return "new";
  }
}

sub chunk_pepfile {
  my ($pep_file, $scratchdir, $size) = @_;

  #Chunk the peptide file
  open (PEPFILE, "$pep_file") or die "couldn't open $pep_file $!";
  my $count = 0;

  #my $chunk = 1;
  my $chunk;
  $last_chunk ? $chunk = $last_chunk + 1 : $chunk = 1;

  #print STDERR "chunking peptide file\n";


  $/ = "\>";
  #print "have opened ".$pep_file."\n";
  while(<PEPFILE>){
    #print $_."\n";
    if ($_ ne "\>") {
      if ($count == 0) {
	open (CHUNK,">".$scratchdir."/chunk.$chunk") or die "couldn't open ".$scratchdir."/chunk.$chunk";
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
