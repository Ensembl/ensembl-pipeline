#!/usr/local/bin/perl -w
=head1 NAME

  delete_transcripts.pl

=head1 SYNOPSIS
 
  delete_transcripts.pl
  deletes transcripts from given database whose ids are passed in a file

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)

    -dbport      For RDBs, what port to connect to (port= in locator)

    -dbname      For RDBs, what name to connect to (dbname= in locator)

    -dbuser      For RDBs, what username to connect as (dbuser= in locator)

    -dbpass      For RDBs, what password to use (dbpass= in locator)

    -dnadbname, dnadbhost, dnadbuser, dnadbpass, dnadbport  - Options for using a reference DNA database

    -file        File containing list of transcript ids to be deleted
=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $dbhost;
my $dbport;
my $dbname;
my $dbuser;
my $dbpass;
my $dnadbhost;
my $dnadbuser;
my $dnadbpass;
my $dnadbport;
my $dnadbname;
my $file;

&GetOptions( 
	     'dbhost:s'      => \$dbhost,
	     'dbport:n'      => \$dbport,
	     'dbname:s'      => \$dbname,
	     'dbuser:s'      => \$dbuser,
	     'dbpass:s'      => \$dbpass,
	     'dnadbhost=s'   => \$dnadbhost,
	     'dnadbname=s'   => \$dnadbname,
	     'dnadbuser=s'   => \$dnadbuser,
	     'dnadbpass=s'   => \$dnadbpass,
	     'dnadbport=s'   => \$dnadbport,
	     'file=s'        => \$file,
	     ) or die("couldn't get options");

&check_options();

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbhost,
					    -user   => $dbuser,
					    -dbname => $dbname,
					    -pass   => $dbpass,
					   );
if($dnadbhost){
  my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor (
						  -host    => $dnadbhost,
						  -user    => $dnadbuser,
						  -dbname  => $dnadbname,
						  -pass    => $dnadbpass,
						  -port    => $dnadbport,
						 );
  $db->dnadb($dnadb);
}

open(FH, $file) or die("couldn't open ".$file);

while(<FH>){
  chomp;
  my $transcript_id= $_;
  # are we about to remove the last transcript from a gene? if so, shout!
  my $sth = $db->prepare("SELECT gene.transcript_count FROM gene, transcript WHERE gene.gene_id=transcript.gene_id AND transcript.transcript_id=$transcript_id");
  $sth->execute;
  my $row = $sth->fetchrow_hashref;
  if($row->{transcript_count} == 1){
    print "STOP! About to delete last transcript\n";
    last;
  }

  print STDERR "about to delete $transcript_id\n";
  my $transcript_adaptor = $db->get_TranscriptAdaptor;
  my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
  $transcript_adaptor->remove($transcript);

}

close FH or die ("Couldn't close [$file]\n");

### SUBROUTINES ###

sub check_options{
  &usage unless defined ($dbname && $dbhost && $dbuser);
  if(defined $dnadbname){
    &usage unless defined ($dnadbhost && $dnadbuser);
  }
  &usage unless defined $file;

}

sub usage{
  print "Usage: delete_transcripts.pl -dbname -dbuser -dbhost -dbpass [-dnadbname -dnadbuser -dnadbhost -dnadbpass] -file\n";
  exit(0);
}
