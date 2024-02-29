#!/usr/bin/env perl 

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 NAME

 load_refseq_synonyms.pl

=head1 SYNOPSIS

 load_refseq_synonyms.pl [arguments]


=head1 DESCRIPTION

 NCBI RefSeq gff files use RefSeq's own identifiers instead of INSDC
 identifiers for contig/scaffold/chromosome since they do their own
 fixes to the sequence. Currently we do not store the RefSeq
 accession in the seq_region_synonym table therefore we have no
 'hook' onto a slice in the database. 

 This script will populate the seq_region_synonym table with a mapping
 between INSDC and RefSeq genomic identifiers.

 Usage:

        perl load_refseq_synonyms.pl -workdir /dir/ \   
                                     -host hostname \
                                     -dbname dbname \
                                     -port port     \
                                     -user username \
                                     -pass password \
                                     -verbose       \
                                     -species scientific_name \
                                     -cleanup     


  -workdir    location used for downloading and writing logfile to

  -verbose    turn on more verbose logging                   

  -species    run on a single species formatted thus: anolis_carolinensis

  -cleanup    clean up files after use


=head1 AUTHOR

 Daniel Barrell <dbarrell@ebi.ac.uk>, Ensembl genebuild team

=head1 CONTACT

 Please post comments/questions to the Ensembl development list
 <http://lists.ensembl.org/mailman/listinfo/dev>

=cut 

use strict;
use warnings;          
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Net::FTP;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::UnmappedObject;

my $workdir = '';
my $verbose = '';     
my $species = '';
my $dbhost  = '';
my $dbname  = '';
my $dbport  = '3306';
my $dbuser  = 'ensro';
my $dbpass  = '';
my $help    = 0;
my $cleanup = 0;
my $ass_checked = 0;
my $cs      = 'toplevel';

&GetOptions(
  'workdir:s' => \$workdir,
  'verbose!'  => \$verbose,
  'species:s' => \$species,
  'dbhost:s'  => \$dbhost,
  'dbuser:s'  => \$dbuser,
  'dbname:s'  => \$dbname,
  'dbpass:s'  => \$dbpass,
  'dbport:n'  => \$dbport,   
  'help'      => \$help,
  'cleanup'   => \$cleanup,
  'coord_system|cs:s' => \$cs, # undocumented, only used for debugging
);  

if (!$workdir or !$species or !$dbhost or !$dbname or !$dbuser or !$dbpass or $help) {
  &usage; exit(1);
} 

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dbport,
  -user    => $dbuser,
  -host    => $dbhost,
  -dbname  => $dbname,
  -pass    => $dbpass
);

my $sa = $db->get_SliceAdaptor;     

my $refseq_external_db_id = get_refseq_external_db_id();
printf STDERR "# External DB ID for RefSeq is: %s\n", $refseq_external_db_id if $verbose;

my $ensembl_assembly_version = get_assembly_name_from_ensembl();
printf STDERR "# Assembly in Ensembl is: %s\n", $ensembl_assembly_version if $verbose;

my $ftphost = "ftp.ncbi.nlm.nih.gov";
my $ftpuser = "anonymous";
my $ftppassword = "";

my $f = Net::FTP->new($ftphost) or die "Can't open $ftphost\n";
$f->login($ftpuser, $ftppassword) or die "Can't log $ftpuser in\n";

# hash lookup of RefSeq identifiers given an INSDC identifier.
my %INSDC2RefSeq;

$species = ucfirst($species);
unless(-e $workdir) {
  my $cmd = "mkdir -p ".$workdir."/".$species;
  my $return = system($cmd);
  if($return) {
    die("Could not create dir for $species\n.Commandline used:\n".$cmd);
  }
}

printf STDERR "Looking for %s accession mapping files\n", $species;

# Scaffold_names file. Always look for this file whether we have
# assembled chromosomes or not. We always look here first to check that we
# are using the same assembly, if we do not then there is no point in continuing
my $fdir = "/genomes/".$species;
if ($f->cwd($fdir)) {
  local $, = ',';
  my @files = $f->ls;
   
  print STDERR "  Found files:\n    ", @files, "\n" if $verbose;

  foreach my $file_to_get ( @files ) {
    next unless $file_to_get =~ /^scaffold_names$/;

    print STDERR "  Downloading $file_to_get...\n";                  
    my $local_file_location = $workdir.'/'.$species.'/'.$file_to_get;
    $f->get($file_to_get, $local_file_location) or throw("Can't get $file_to_get from $fdir\n");
    print STDERR "  Download complete.\n" if $verbose;
    #TODO change: at the moment you need to be in the workdir

    open(ACC,"<",$local_file_location) or throw("Could not open $local_file_location");
    while(<ACC>) {
      my $line = $_;
      next if $line =~ /^#/;
      my ($ass, $insdc, $refseq) = ('','');
      if ($line =~ /^(.*)\t.*\t(.*)\t(.*)\t.*$/) {
        ($ass, $refseq, $insdc) = ($1, $2, $3);
        unless ($ass_checked) {
          $ass_checked = "true";
          $ass =~ s/[^1-9]//g; # to ensure matching of RNOR5 and Rnor_5.0  
          if ($ass ne $ensembl_assembly_version) {
            throw("Assembly in Ensembl ($ensembl_assembly_version) is different to what RefSeq has annotated ($ass)")
          }
        }
        $INSDC2RefSeq{$insdc} = $refseq
        unless exists $INSDC2RefSeq{$insdc};
      } else {
        throw ("line malformed:\n",$line)
      }
    }
    close ACC;
  }
} else {
  printf STDERR "%s: No scaffold_name file found - This species not on RefSeq ftp site?\n",$species;
}                                                                                     

# Does this organism have assembled chromosomes?
$fdir = "/genomes/".$species."/Assembled_chromosomes";
if ($f->cwd($fdir)) {
  # there are assembled chromosomes available for this species
  local $, = ',';
  my @files = $f->ls;
  
  if ($verbose) {
    print STDERR "  Found files:\n    ", @files, "\n";
  }

  foreach my $file_to_get ( @files ) {
    next unless $file_to_get =~ /^.*_accessions_.*$/;

    print STDERR "  Downloading $file_to_get...\n";                  
    my $local_file_location = $workdir.'/'.$species.'/'.$file_to_get;
    $f->get($file_to_get, $local_file_location)
      or die "Can't get $file_to_get from $fdir\n";
    print STDERR "  Download complete.\n" if $verbose;

    open(ACC,"<",$local_file_location) or throw("Could not open $local_file_location");
    while(<ACC>) {
      my $line = $_;
      next if $line =~ /^#/;
      my ($insdc, $refseq) = ('','');
      if ($line =~ /^\S+\t(\S+)\t\S+\t(\S+)\t\S+$/) {
        ($refseq, $insdc) = ($1, $2);
        $INSDC2RefSeq{$insdc} = $refseq
        unless exists $INSDC2RefSeq{$insdc};
      } else {
        throw ("line malformed:\n",$line)
      }
    }       
    close ACC;
  }
} else {
  printf STDERR "%s: No Asssembled_chromosome folder on refseq ftp site\n",$species;
}                                                                                     

if ($cleanup) {
  `rm -rf $species`;
  print STDERR "Removed mapping directory\n";
}

$f->quit;


# Now we have a mapping, update the database:
foreach my $insdc (keys %INSDC2RefSeq) {

  my $refseq = $INSDC2RefSeq{$insdc};

  # For anolis everything had a version except AAWZ.* strip the version for these
  $insdc =~ s/\.\d+$// if $insdc =~ /^AAWZ.*/;
  next if $insdc =~ /^EU.*/; # basically, skip it

  # TODO we may not want to skip MT for some species
  
  my $slice = $sa->fetch_by_region($cs, $insdc);

  # some accesions found are not on the assembly we have loaded, if so, report and move on
  if (!defined $slice) {
    print "Warning: " . $insdc . " not found in the assembly that is loaded in " . $dbname . "\n";
    next;
  }

  # TODO may be the case that some species do not have chromosome accessions.

  # write the synonyms directly to the database
  $slice->add_synonym($refseq, $refseq_external_db_id);
  $sa->update($slice); 
}

exit 0;


# SUBS

# not necessarily needed - external db id should be constant 
sub get_refseq_external_db_id {
  my $sth_refseq = $db->dbc->prepare('SELECT external_db_id 
    FROM external_db
    WHERE db_name = "RefSeq_genomic"');
  $sth_refseq->execute();
  my ($refseq_db_id) = $sth_refseq->fetchrow_array;
  throw("No RefSeq DB in Ensembl") unless defined $refseq_db_id;
  return $refseq_db_id;
}

# Sometimes NCBI has a different version of the assembly - we 
# need to check and to skip these, since it would be a waste of time
sub get_assembly_name_from_ensembl {
  my $sth_assembly = $db->dbc->prepare('SELECT meta_value
    FROM meta
    WHERE meta_key = "assembly.name"');
  $sth_assembly->execute();
  my ($assembly_version) = $sth_assembly->fetchrow_array;
  throw("No Assembly version in Ensembl meta table") unless defined $assembly_version;
  $assembly_version =~ s/[^1-9]//g; # to ensure matching of RNOR5 and Rnor_5.0
  return $assembly_version;
}

sub usage {
    print <<EOF

  Usage:

    $0 -workdir <workdir> -sqlfile <sqlfile> -species <scientific_name> -dbhost <dbhost> [-dbport <dbport>] -dbname <dbname> -dbuser <dbuser> -dbpass <dbpass> [-verbose] [-help]

    -workdir  		Local path where work will take place

    -species      which species to run e.g. 'homo_sapiens'

    -dbhost    		host name where the reference database will be created

    -dbport    		what port to connect (default 3306)

    -dbname    		name for the new reference db that will be created 

    -dbuser    		what username to connect as

    -dbpass    		what password to use

    -verbose      Use this option to get more print statements to follow
                  the script. Set to 0 (not verbose) by default to get
                  only the final summary.

    -help		      Show usage.

    -cleanup      Remove files downloaded and species directory after creation of mapping

  Example:

  bsub -M 1000 -R 'select[mem>1000] rusage[mem=1000]' -o refseq_synonym.out -e refseq_synonym.err "perl load_refseq_synonyms.pl -workdir ./refseq_synonyms -species anolis_carolinensis -dbhost host -dbname anolis_carolinensis_core_75_2 -dbuser *** -dbpass *** -verbose"

EOF
}

1; 
