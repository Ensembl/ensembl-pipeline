#!/usr/local/bin/perl -w
=head1 NAME

  delete_genes.pl

=head1 SYNOPSIS
 
  delete_genes.pl
  deletes genes from given database whose ids are passed in through STDIN

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)
    -help      summary of options


=head2 EXAMPLES

./delete_genes.pl -dbhost ecs2b -dbuser ensadmin -dbpass **** -dbname rat_Jun03_mk2 genes_to_delete

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host;
my $port;
my $dbname;
my $dbuser;
my $pass;
my $help = 0;

&GetOptions( 
	     'dbhost:s'      => \$host,
	     'dbport:n'      => \$port,
	     'dbname:s'    => \$dbname,
	     'dbuser:s'    => \$dbuser,
	     'dbpass:s'      => \$pass,
	     'help!'      => \$help, 
	     )or &useage();


if(!$host || !$dbuser || !$dbname ){
  $help = 1;
}

&useage() if($help);
 


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
					    -user => $dbuser,
					    -dbname => $dbname,
					    -pass  => $pass,
					    -port => $port,
					   );

while(<>){
  chomp;
  my $gene_id= $_;
  print STDERR "about to delete $gene_id\n";
  #my $sth = $db->prepare("delete from gene where gene_id = $gene_id");
  #$sth->execute;
  eval{
    my $gene_adaptor = $db->get_GeneAdaptor;
    my $gene = $gene_adaptor->fetch_by_dbID($gene_id);
    $gene_adaptor->remove($gene);
  };
  if($@){
    print STDERR "$@ couldn't remove gene\n";
    my $sth = $db->prepare("delete from gene where gene_id = $gene_id");
    $sth->execute;
  }
}



sub useage{
	exec('perldoc', $0);
	exit;
}
