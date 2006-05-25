#!/usr/local/ensembl/bin/perl -w

=head1 NAME

load_medaka.pl 

=head1 DESCRIPTION

Checks that all CDS have length divisible by 3 
-> shortens the end of the first exon by 1 if this is the case and if the exon is on the (-) strand.

=cut

use strict;
use MedakaConf;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning);
my %opt; 

# options submitted with commandline override MedakaConf.pm 
GetOptions(
           \%opt ,
           '-h|help' ,  
           'gff_file=s' ,
	   'outfile=s' ,
           ) ; 

$MED_GFF_FILE     = $opt{gff_file} if $opt{gff_file} ; 
$MED_FIXED_GFF = $opt{outfile} if $opt{outfile} ;


# end of set-up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Parse GFF file
print STDERR "Fixing $MED_GFF_FILE...\n";
fix_gff(); 


print STDERR "Done!\n";



# start of subroutines
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Parse GFF file. Use a hash of hashes.
# OuterKey = gene_stable_id
# InnerKey = exon_number
sub fix_gff {
  # example of GFF file.
  #	scaffold1       UTOLAPRE05100100001     initial-exon    129489  129606  .       +       0
  #	scaffold1       UTOLAPRE05100100001     internal-exon   129920  130027  .       +       2
  #	scaffold1       UTOLAPRE05100100001     internal-exon   130753  130839  .       +       2
  #	scaffold1       UTOLAPRE05100100001     final-exon      131859  132262  .       +       2
  #	scaffold6469    UTOLAPRE05100120178     single-exon     1604    2746    .       -       0

  # read in the file. 
  my @gff;
  open (GFF,$MED_GFF_FILE) or die "Cannot open gff file $MED_GFF_FILE\n";
  open (FIXED,">>$MED_FIXED_GFF") or die "Cannot open gff file $MED_FIXED_GFF\n";
  my $line;
  while (<GFF>){
    chomp;
    $line = $_;
    next if ($line =~ m/^\#/);
    my @fields = split/\s+/, $line;
    # [scaffold, gene_id, feature, start, end, ".", strand, frame]
    push @gff, [@fields]; 
    # make the Gene obj and load into $output_db
    if (($fields[2] eq 'final-exon' && $fields[6] eq '+') ||
        ($fields[2] eq 'initial-exon' && $fields[6] eq '-') ||
	($fields[2] eq 'single-exon')) {
      # check if CDS is a multiple of 3  
      my @fixed = @{fix_cds(\@gff)};
      # print to outfile
      for (my $i=0; $i<scalar(@fixed); $i++) {
        for (my $j=0; $j<scalar(@{$fixed[0]}); $j++) {
          print FIXED $fixed[$i][$j]."\t";
	  print STDOUT $fixed[$i][$j]."\t";
	}
	print FIXED "\n";
	print STDOUT "\n";
      } 
      @gff = ();
      @fixed = ();
    }
    $line = '';
  }    
  close GFF; 
  close FIXED;
}

sub fix_cds {
  my ($gff) = @_;
  my @gene = @{$gff};
  my $cds_length = 0;
  my $exon_length;  
  for (my $i=0; $i<scalar(@gene); $i++){
    # [scaffold, gene_id, feature, start, end, ".", strand, frame]
    $exon_length = $gene[$i][4] - $gene[$i][3] +1;
    $cds_length += $exon_length;
  }
  if ($cds_length%3 != 0) { 
    print STDERR "CDS problem: ".$gene[0][1]."\n" ;
    if ($gene[-1][6] eq '-') {
      my $saved = $gene[-1][3];
      $gene[-1][3] = $saved + 1;
      print STDOUT "Fixed $saved\n";
    } else {
      print STDERR "  Have not fixed ".$gene[0][1]."\n";
    }
  }
  return \@gene;    
}    
  

