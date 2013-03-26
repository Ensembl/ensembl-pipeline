#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/Pre-release-site/set_stable_ids.pl,v $
# $Revision: 1.12 $

=head1 NAME

  set_stable_ids.pl

=head1 SYNOPSIS

  Script to create stable IDs for the pre-release targetted gene set.
  Sets stable IDs to be the name of the evidence used to build the
  structure, plus a .version in case more than one gene is built from
  the same protein.

=head1 DESCRIPTION


=head1 OPTIONS

  --dbhost  Host name for database.
  --dbport  For RDBs, what port to connect to (optional).
  --dbname  For RDBs, what name to connect to.
  --dbuser  For RDBs, what username to connect as.
  --dbpass  For RDBs, what password to use.

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );

my $host;
my $port = '3306';
my $name;
my $user;
my $pass;

GetOptions( 'dbhost|h|host:s' => \$host,
            'dbport|P|port:n' => \$port,
            'dbname|D|db:s' => \$name,
            'dbuser|u|user:s' => \$user,
            'dbpass|p|pass:s' => \$pass, );

my %proteins;
my %transcript_proteins;
my %exons;
my $db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( '-host'   => $host,
                                      '-user'   => $user,
                                      '-dbname' => $name,
                                      '-pass'   => $pass,
                                      '-port'   => $port, );

my $query =
  'UPDATE gene SET stable_id = ?, version = ?, ' .
  'created_date = now(), modified_date = now() WHERE gene_id = ?;';

my $sth = $db->dbc->prepare($query);

my $trans_query =
  'UPDATE transcript SET stable_id = ?, version = ?, ' .
  'created_date = now(), modified_date = now() ' .
  'WHERE transcript_id = ?;';

my $trans_sth = $db->dbc->prepare($trans_query);

my $translation_query =
  'UPDATE translation SET stable_id = ?, version = ?, ' .
  'created_date = now(), modified_date = now() ' .
  'WHERE translation_id = ?;';

my $translation_sth = $db->dbc->prepare($translation_query);

my $exon_query =
  'UPDATE exon SET stable_id = ?, version = ?, ' .
  'created_date = now(), modified_date = now() WHERE exon_id = ?;';

my $exon_sth = $db->dbc->prepare($exon_query);

# Arrays to hold write execution so that write only happens when 
# there are no fails, i.e. no partial writes possible
my @trans_write = ();
my @translation_write = ();
my @exon_write = ();

GENE:
foreach my $gene_id ( @{ $db->get_GeneAdaptor->list_dbIDs() } ) {
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);

  if ( defined( $gene->stable_id() ) ) { next GENE }

  my $gene_protein_id = get_gene_protein_id($gene);

  $proteins{$gene_protein_id}++;

  eval {
    
    foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
      my $transcript_protein_id =
        get_transcript_protein_id($transcript);

      $transcript_proteins{$transcript_protein_id}++;

      push(@trans_write,$transcript_protein_id . "_" .$transcript_proteins{$transcript_protein_id}.
	             ":\n:".$transcript_proteins{$transcript_protein_id}.":\n:".$transcript->dbID());


      if($transcript->translation())
      {
        my $translation = $transcript->translation();


        push(@translation_write,$transcript_protein_id . "_" .$transcript_proteins{$transcript_protein_id}.
	                    ":\n:".$transcript_proteins{$transcript_protein_id}.":\n:".$translation->dbID());
      }

      my $exon_count = 1;

    EXON:
      foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
        my $exon_dbID = $exon->dbID();

        if ( exists( $exons{$exon_dbID} ) ) { next EXON }

        my $stable_id =
          $transcript_protein_id . "_" .
          $transcript_proteins{$transcript_protein_id} . "." .
          $exon_count++;


        push(@exon_write, $stable_id.":\n:".$transcript_proteins{$transcript_protein_id}.":\n:".$exon_dbID);

        $exons{$exon_dbID} = $stable_id;

      }


    } ## end foreach my $transcript ( @{...})
  };

  if ($@) {
    throw( "Stable ID insertion for gene " . $gene_id . " failed $@" );
  }

  else
  {
      my $gene_write_string = $gene_protein_id . "_" . $proteins{$gene_protein_id}.
	                       ":\n:". $proteins{$gene_protein_id}.":\n:".$gene_id;

      write_to_db($gene_id,$gene_write_string);
      clear_write_arrays();
  }

} ## end foreach my $gene_id ( @{ $db...})

sub get_gene_protein_id {
  my ($gene) = @_;

  foreach my $exon ( @{ $gene->get_all_Exons() } ) {
    foreach my $sf ( @{ $exon->get_all_supporting_features() } ) {
      my $protein_id = $sf->hseqname();
      if ( defined($protein_id) ) { return $protein_id }
    }
  }

  throw( "Found no protein id for " . $gene->dbID() );
}

sub get_transcript_protein_id {
  my ($transcript) = @_;

  foreach my $sf ( @{ $transcript->get_all_supporting_features() } ) {
    my $protein_id = $sf->hseqname();
    if ( defined($protein_id) ) { return $protein_id }
  }

  throw( "Found no protein id for " . $transcript->dbID() );
}

sub write_to_db
{
  my ($gene_id,$gene_write_string) = @_;

  my @tmp = split(':\n:',$gene_write_string);
  
  $sth->execute($tmp[0],$tmp[1],$tmp[2]);

  @tmp = ();

  foreach (@trans_write)
  {
      @tmp = split(':\n:',$_);
      $trans_sth->execute($tmp[0],$tmp[1],$tmp[2]);
  }

  @tmp = ();

  foreach (@translation_write)
  {
      @tmp = split(':\n:',$_);
      $translation_sth->execute($tmp[0],$tmp[1],$tmp[2]);
  }
          
  @tmp = ();

  foreach (@exon_write)
  {
      @tmp = split(':\n:',$_);
      $exon_sth->execute($tmp[0],$tmp[1],$tmp[2]);
  }

}

sub clear_write_arrays
{
    @trans_write = ();
    @translation_write = ();
    @exon_write = ();
}
