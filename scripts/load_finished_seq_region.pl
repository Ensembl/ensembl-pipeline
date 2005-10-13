#!/usr/local/ensembl/bin/perl -w

=head1 NAME

load_finished_seq_region.pl

=head1 SYNOPSIS

  load_finished_seq_region.pl 

=head1 DESCRIPTION

This script uses the entries in an agp file to load the chromosome level data into the seq_region, seq_region_attrib, assembly, dna tables of the pipeline database. Appropriate options like the coord_system version, rank for the specific chromosome coord_system level, has to be given from the coord_system table of the corresponding database.
a set of seq_regions into the seq_region object. The sequence that is loaded into the dna table is Pfetched. 

format of an agp file is shown here,

http://www.sanger.ac.uk/Projects/C_elegans/DOCS/agp_files.shtml

here is an example commandline

this will load the assembled pieces from the agp file into the seq_region
and other sequence region related tables.

./load_finished_seq_region.pl -dataset_name human -chromosome_cs_version Otter -chromosome_name chr22 -chromosome_cs_rank 2 -agp_file  chr22.agp -verbose -dbhost ecs4 -dbport 3352 -dbname human22 -dbuser pipuser -dbpass *****




=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -chromosome_name the name of the coordinate system being stored
    -chromosome_cs_version the version of the coordinate system being stored
    -chromosome_cs_rank the rank of the coordinate system.  
    -agp_file the name of the agp file to be parsed
    -verbose, prints the name which is going to be used can be switched 
              off with -noverbose
    -help      displays this documentation with PERLDOC

=head1 CONTACT

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk

=head1 APPENDIX

to be added :
-check for type of input (list vs agp)
-check if chromosome is already present
-check if contig is already present
-supercontigs
-transactions
-add loading to otter database parallely

=cut

use strict;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SequenceAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Attribute;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host   = '';
my $port   = '';
my $dbname = '';
my $dbuser = '';
my $dbpass = '';
my $help;
my $ds_name;
my $chromosome_name;
my $chromosome_cs_version;
my $chromosome_cs_rank;
my $default = 0;
my $agp;
my $verbose = 0;
my $INPUT_TYPE;

&GetOptions(
            'dbhost:s'   => \$host,
            'dbport:n'   => \$port,
            'dbname:s'   => \$dbname,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'dataset_name:s' => \$ds_name,
            'chromosome_cs_version:s' => \$chromosome_cs_version,
	    'chromosome_name:s' => \$chromosome_name,
            'chromosome_cs_rank:i' => \$chromosome_cs_rank,
            'agp_file:s' => \$agp,
            'verbose!' => \$verbose,
            'h|help'     => \$help,
           ) or ($help = 1);

if(!$host || !$dbuser || !$dbname || !$dbpass){
  print STDERR "Can't store sequence without database details\n";
  print STDERR "-dbhost $host -dbuser $dbuser -dbname $dbname ".
    " -dbpass $dbpass\n";
  $help = 1;
}

if(!$ds_name || !$chromosome_cs_version || !$chromosome_name || !$chromosome_cs_rank || !$agp){
  print STDERR "Need dataset_name and chromosome_name and chromosome_cs_version and chromosome_cs_rank and agp_file to be able to run\n";
  print STDERR "-dataset_name $ds_name -chromosome_name $chromosome_name -chromosome_cs_version $chromosome_cs_version -chromosome_cs_rank $chromosome_cs_rank -agp_file $agp\n";
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}

#change chromosome name from chr22 to 22 
$chromosome_name =~ s/^chr//i;


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbname => $dbname,
    -host   => $host,
    -user   => $dbuser,
    -port   => $port,
    -pass   => $dbpass
);
my $dbc= $db->dbc();
my $aa = $db->get_AttributeAdaptor();

my $csa = $db->get_CoordSystemAdaptor();
my $cs;
my $cs1;
eval{
  $cs = $csa->fetch_by_name("chromosome", $chromosome_cs_version, $chromosome_cs_rank);
};
if(!$cs){
    print STDERR "A coord_system matching the arguments does not exsist in the cord_system table, please ensure you have the right coord_system entry in the database";
}
my $sa  = $db->get_SliceAdaptor();
my $seqa = $db->get_SequenceAdaptor();
my %agp_fileparams;
my $href=\%agp_fileparams;
my $asm_seq_reg_id;
$href= &insert_chromosome_seq_region($agp,$sa,$cs,\%agp_fileparams,\$asm_seq_reg_id);
$csa = $db->get_CoordSystemAdaptor();
eval{
  $cs = $csa->fetch_by_name("clone", '', 3);
};
if(!$cs){
    print STDERR "A coord_system matching the arguments does not exsist in the cord_system table, please ensure you have the right coord_system entry in the database";
}
eval{
  $cs1 = $csa->fetch_by_name("contig", '', 4);
};
if(!$cs1){
    print STDERR "A coord_system matching the arguments does not exsist in the cord_system table, please ensure you have the right coord_system entry in the database";
}
$sa  = $db->get_SliceAdaptor();

&insert_clone_contig_seq_region_attrib_dna_assembly($agp,$sa,$cs,$cs1,$href,$seqa,$asm_seq_reg_id,$dbc,$aa);


################start of subroutines

sub insert_chromosome_seq_region{

    my ($agp_file,$sa,$cs,$hashref,$asm_seq_reg_id_ref) = @_;
    my %end_value;
    open(FH, $agp_file) or throw("Can't open ".$agp_file." ".$!);
  LINE:while(<FH>){   
      chomp;
      next if /^\#/;

      #chr22   14430001        14467693        4       F       AP000522.1      1       37693   +
      #chr22   14467694        14506727        5       F       AP000523.1      2273    41306   +
      my @values = split  /\s+/, $_, 10;
      if($values[4] eq 'N'){
	  next LINE; 
      }
      if ($values[4] eq 'F') {
	  $INPUT_TYPE = "agp";
      }
      else {
	  die "please check file format";
      }
      $values[0] =~ s/^chr//i;
      if ($values[0] != $chromosome_name) {
	  die " the supplied chromosome name does not match with the agp file \n\n";
      }
      if(!$end_value{$values[0]}){
	  $end_value{$values[0]} = $values[2];
      }else{
      
	  if($values[2] > $end_value{$values[0]}){
	      $end_value{$values[0]} = $values[2];
	  }
      }
      $hashref->{$values[3]} = \@values;

  }//while
      
      my $endv;
      foreach my $name (keys(%end_value)){
	  $endv = $end_value{$name};
	  my $slice = &make_slice($name, 1, $endv, $endv, 1, $cs);
	  $$asm_seq_reg_id_ref=$sa->store($slice);
      }
  
      close(FH) or throw("Can't close ".$agp_file);
      return $hashref;
}

sub insert_clone_contig_seq_region_attrib_dna_assembly{

    my ($agp_file,$sa,$cs,$cs1,$hashref,$seqa,$asm_seq_reg_id,$dbc,$aa) = @_;
    my $insert_query = qq {INSERT INTO assembly (asm_seq_region_id, cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori) values (?,?,?,?,?,?,?)};
    my $sth=$dbc->prepare($insert_query);

    while( my ($k, $v) = each %$hashref ) {
	my @values = @$v;
	my $acc_ver=$values[5];
	##fetch the dna sequence from pfetch server with acc_ver id
	my $seqobj ||= &pfetch_acc_sv($acc_ver);
	my $seq   =$seqobj->seq;
	my $seqlen= $seqobj->length;
	my $clone;
	my $clone_seq_reg_id;
        #split into accesion number and version number 
        my ($acc, $sv) = $acc_ver =~ /^(.+)\.(\d+)$/;
        die "Can't parse '$acc_ver' into accession and sv\n" unless $acc and $sv;
        eval { $clone = $sa->fetch_by_region('clone', $acc_ver); };
	if($clone){
	    warn "clone <".$clone->name."> is already in the pipeline database\n";
	    $clone_seq_reg_id = $clone->get_seq_region_id;
### code to be added to grab contigs related to the clone and alter the contig insertion section by moving it to the else part of this.

	}
	else {
	    ##make clone and insert clone to seq_region table
	    my $slice = &make_slice($acc_ver,1,$seqlen,$seqlen,1,$cs);
	    $clone_seq_reg_id = $sa->store($slice);
	    if (!defined ($clone_seq_reg_id)) {
		print "clone seq_region_id has not been returned for the accession $acc_ver";
		exit;
	    }

	    ##make clone attributes
	    my @attrib;
	    my $attrib = &make_attribute('htgs_phase','HTGS Phase','High Throughput Genome Sequencing Phase','3');
	    push @attrib, $attrib;
	    push @attrib, &make_attribute('intl_clone_name','International Clone Name','',$acc_ver);
	    push @attrib, &make_attribute('embl_accession','EMBL Accession','',$acc);
	    push @attrib, &make_attribute('embl_version','EMBL Version','',$sv);
	    $aa->store_on_Slice($slice,\@attrib);
	}
	my $chr_name = $values[0];
	my $chr_start= $values[1];
	my $chr_end  = $values[2];
	my $ctg_start= $values[6];
	my $ctg_end  = $values[7];
	my $contig = $acc_ver."."."1".".".$seqlen;
	##make contig and insert contig, and associated dna sequence to seq_region & dna table
	my $slice = &make_slice($contig,1,$seqlen,$seqlen,1,$cs1);
	my $ctg_seq_reg_id = $sa->store($slice,\$seq);
	if (!defined ($ctg_seq_reg_id)) {
	    print "contig seq_region_id has not been returned for the contig $contig";
	    exit;
	}
       	#translate orientation to integer
	my $ctg_ori;
        if ($values[8] eq '-') {
            $ctg_ori = -1;
        }
        elsif ($values[8] eq '+') {
            $ctg_ori = 1;
        }
	else {
            die "Invalid orientation $ctg_ori\n";
        }	
	##insert chromosome to contig assembly data into assembly table
	$sth->execute($asm_seq_reg_id,$ctg_seq_reg_id,$chr_start,$chr_end,$ctg_start,$ctg_end,$ctg_ori);
	##insert clone to contig assembly data into assembly table
	$sth->execute($clone_seq_reg_id,$ctg_seq_reg_id,1,$seqlen,1,$seqlen,$ctg_ori);
	
    }##while

}##sub

sub make_attribute{
    my ($code,$name,$description,$value) = @_;
    my $attrib = Bio::EnsEMBL::Attribute->new
	(
	 -CODE => $code,
	 -NAME => $name,
	 -DESCRIPTION => $description,
	 -VALUE => $value
	 );
    return $attrib;
}
	 

sub make_slice{

  my ($name, $start, $end, $length, $strand, $coordinate_system) = @_;
  my $slice = Bio::EnsEMBL::Slice->new
      (
       -seq_region_name   => $name,
       -start             => $start,
       -end               => $end,
       -seq_region_length => $length,
       -strand            => $strand,
       -coord_system      => $coordinate_system,
      );
  return $slice;
}

#
# Pfetch the sequences
#-------------------------------------------------------------------------------
#  If the sequence isn't available from the default pfetch
#  the archive pfetch server is tried.
#

{
    my( $pfetch, $pfetch_archive );

    sub pfetch_acc_sv {
        my( $acc_ver ) = @_;
        
        $pfetch         ||= Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
        $pfetch_archive ||= Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new(

            -PFETCH_PORT => 23100,
            );
        my $seq = $pfetch->get_Seq_by_acc($acc_ver);
        unless ($seq) {
            warn "Fetching '$acc_ver' from archive\n";
            $seq = $pfetch_archive->get_Seq_by_acc($acc_ver);
        }
	unless ($seq){
            my $seq_file = "$acc_ver.seq";
	    warn "Attempting to read fasta file <$acc_ver.seq> in current dir.\n";
	    my $in = Bio::SeqIO->new(
                -file   => $seq_file,
                -format => 'FASTA',
                );
            $seq = $in->next_seq;
            my $name = $seq->display_id;
            unless ($name eq $acc_ver) {
                die "Sequence in '$seq_file' is called '$name' not '$acc_ver'";
            }
	}

        return $seq;
	
    }
}
1;




