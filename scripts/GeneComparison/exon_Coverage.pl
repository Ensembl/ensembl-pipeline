#!/usr/local/ensembl/bin/perl


=head1 NAME

comapre_Exons
 
=head1 DESCRIPTION

reads the config options from Bi::Ensembl::Pipeline::GeneComparison::GeneCompConf
and reads as input an input_id in the style of other Runnables, i.e. -input_id chr_name.chr_start-chr_end

=head1 OPTIONS

    -input_id  The input id: chrname.chrstart-chrend

=cut

use strict;  
use diagnostics;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison;
use Getopt::Long;

## load all the parameters
use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf;

# annotation
my $host1   = $DBHOST1;
my $dbname1 = $DBNAME1;
my $path1   = $PATH1;
my $port1   = $PORT1 || 3306;
my $type1   = $GENETYPES1;
my $user1   = $DBUSER1;

# prediction
my $host2   = $DBHOST2;
my $dbname2 = $DBNAME2;
my $path2   = $PATH2;
my $port2   = $PORT2 || 3306;
my $type2   = $GENETYPES2;
my $user2   = $DBUSER2;

# reference db
my $ref_host1   = $REF_DBHOST1;
my $ref_dbname1 = $REF_DBNAME1;
my $ref_port1   = $REF_PORT1 || 3306;
my $ref_path1   = $REF_PATH1;
my $ref_user1   = $REF_DBUSER1;



my $ref_host2   = $REF_DBHOST2;
my $ref_dbname2 = $REF_DBNAME2;
my $ref_port2   = $REF_PORT2 || 3306;
my $ref_path2   = $REF_PATH2;
my $ref_user2   = $REF_DBUSER2;

my $input_id;
my $write  = 0;
my $check  = 0;
my $params;
my $pepfile;
my $lower_bound = 0;

# can override db options on command line
&GetOptions( 
	    'input_id:s'    => \$input_id,
	    'lower_bound:n' => \$lower_bound,
	   );
	
unless( $input_id){     
  print STDERR "Usage: exon_Coverage.pl -input_id <chrname>.<chrstart>-<chrend>\n";
  print STDERR "                        -lower_bound <lower_bound> (defaults to zero)\n";
  exit(0);
}
    
# get genomic region 
my $chr      = $input_id;
$chr         =~ s/\.(.*)-(.*)//;
my $chrstart = $1;
my $chrend   = $2;

unless ( $chr && $chrstart && $chrend ){
       print STDERR "bad input_id option, try something like 20.1-5000000\n";
}

# connect to the databases 
my $dna_db1= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $ref_host1,
					       -user  => $ref_user1,
					       -port  => $ref_port1,
					       -dbname=> $ref_dbname1);
$dna_db1->assembly_type($ref_path1); 

my $dna_db2= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $ref_host2,
						-user  => $ref_user2,
					        -port  => $ref_port2,
						-dbname=> $ref_dbname2);
$dna_db2->assembly_type($ref_path2);


my $db1= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host1,
					    -user  => $user1,
					    -port  => $port1,
					    -dbname=> $dbname1,
					    -dnadb => $dna_db1,
					   );
print STDERR "Connected to database $dbname1 : $host1 : $user1 \n";


my $db2= new Bio::EnsEMBL::DBSQL::DBAdaptor(-host  => $host2,
					    -user  => $user2,
					    -port  => $port2,
					    -dbname=> $dbname2,
					    -dnadb => $dna_db2);
print STDERR "Connected to database $dbname2 : $host2 : $user2 \n";



# use different golden paths
$db1->assembly_type($path1); 
$db2->assembly_type($path2); 

my $sa1 = $db1->get_SliceAdaptor;
my $sa2 = $db2->get_SliceAdaptor;




# get a virtual contig with a piece-of chromosome #
my ($slice1,$slice2);

print STDERR "Fetching region $chr, $chrstart - $chrend\n";
$slice1 = $sa1->fetch_by_chr_start_end($chr,$chrstart,$chrend);
$slice2 = $sa2->fetch_by_chr_start_end($chr,$chrstart,$chrend);

# get the genes of type @type1 and @type2 from $slice1 and $slice2, respectively #
my (@genes1,@genes2);

foreach my $type ( @{ $type1 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = @{$slice1->get_all_Genes_by_type($type)};
  my @more_trans = ();
  foreach my $gene ( @more_genes ){
    push ( @more_trans, @{$gene->get_all_Transcripts} );
  }
  push ( @genes1, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
  print STDERR "with ".scalar(@more_trans)." transcripts\n";
}


foreach my $type ( @{ $type2 } ){
  print STDERR "Fetching genes of type $type\n";
  my @more_genes = @{$slice2->get_all_Genes_by_type($type)};
  my @more_trans = ();
  foreach my $gene ( @more_genes ){
    push ( @more_trans, @{$gene->get_all_Transcripts} );
  }
  push ( @genes2, @more_genes ); 
  print STDERR scalar(@more_genes)." genes found\n";
  print STDERR "with ".scalar(@more_trans)." transcripts\n";
}


# get a GeneComparison object 
my $gene_comparison = 
  Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison->new(					     
							      '-annotation_genes' => \@genes1,
							      '-prediction_genes' => \@genes2,
							      '-input_id'         => $input_id,
							     );


#########################################################

$gene_comparison->exon_Coverage(\@genes1, \@genes2, $lower_bound);


