#!/usr/local/bin/perl -w
=head1 NAME

  set_stable_ids.pl

=head1 SYNOPSIS
 
  Script to create stable ids for the pre-release targetted gene set.
  Sets stable ids to be the name of the evidence used to buikld the structure, 
  plus a .version in case more than one gene is built from the same protein.

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)

    -dbport      For RDBs, what port to connect to (port= in locator)

    -dbname      For RDBs, what name to connect to (dbname= in locator)

    -dbuser      For RDBs, what username to connect as (dbuser= in locator)

    -dbpass      For RDBs, what password to use (dbpass= in locator)

=cut
use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host;
my $port;
my $name;
my $user;
my $pass;

&GetOptions( 
	     'dbhost:s'      => \$host,
	     'dbport:n'      => \$port,
	     'dbname:s'      => \$name,
	     'dbuser:s'      => \$user,
	     'dbpass:s'      => \$pass,
	     );

my %proteins;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $user,
					    '-dbname' => $name,
					    '-pass'   => $pass,
					    '-port'   => $port,
					   );

my $query = "insert into gene_stable_id values (?,?,?,now(),now())";
my $sth = $db->prepare($query);

GENE:


foreach my $gene_id(@{$db->get_GeneAdaptor->list_geneIds}) {
  print STDERR "gene id $gene_id\n";
  
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
  my @exons = @{$gene->get_all_Exons};
  my @sf = @{$exons[0]->get_all_supporting_features};

  my $protein_id;
  foreach my $f(@sf){
    last if defined $protein_id;
    eval{
      $protein_id = $f->hseqname;
    }
  }
  print "protein_id: $protein_id\n";
  next GENE unless defined $protein_id;
  $proteins{$protein_id}++;

  $sth->execute($gene_id,$protein_id, $proteins{$protein_id});

}
