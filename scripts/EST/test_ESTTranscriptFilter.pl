#!/usr/local/ensembl/bin/perl -w
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::Runnable::ESTTranscriptFilter;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $dbname;
my $dbhost;
my $input_id;

&GetOptions( 
	     'input_id:s'    => \$input_id,
	     'dbname:s'      => \$dbname,
             'dbhost:s'      => \$dbhost,
	     );

unless ( $input_id && $dbname && $dbhost ){
    print STDERR "script to read ESTs and pass them through\n";
    print STDERR "ESTTranscriptFilter.pm as a test\n";
    print STDERR "$0 -input_id -dbname -dbhost\n";
    exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host             => $dbhost,
					    -user             => 'ensro',
					    -dbname           => $dbname,
					    );

$input_id =~ /(\S+)\.(\d+)-(\d+)/;
my $chrname  = $1;
my $chrstart = $2;
my $chrend   = $3;
my $slice    = $db->get_SliceAdaptor->fetch_by_chr_start_end($chrname,$chrstart,$chrend);

print STDERR "fetching ests in region $chrname.$chrstart-$chrend\n";
my @genes  = @{$slice->get_all_Genes_by_type('exonerate')};

my @forward_trans;
my @reverse_trans;
my @unknown;

foreach my $gene ( @genes ){
    foreach my $tran ( @{$gene->get_all_Transcripts} ){
	my @exons1 = @{$tran->get_all_Exons};
	#print STDERR "transcript ".$tran->dbID." exons ".scalar(@exons1)."\n";

	my $strand = $exons1[0]->strand;
	#print STDERR "strand = $strand\n";
	if ( $strand == 1 ){
	    push ( @forward_trans, $tran );
	}
	elsif( $strand == -1 ){
	    push ( @reverse_trans, $tran );
	}
	else{
	    push ( @unknown, $tran );
	}
    }
}
if ( @unknown ){
    print STDERR scalar(@unknown)." ests with no strand ??\n";
}


my @all_trans;
if ( @forward_trans ){
    push ( @all_trans, @forward_trans );
}
if ( @reverse_trans ){
    push ( @all_trans, @reverse_trans );
}

print_GFF( "input_est.gff", \@all_trans, "input");

@all_trans = ();

my $est_filter = Bio::EnsEMBL::Pipeline::Runnable::ESTTranscriptFilter
    ->new( -coverage => 95,
	   -perc_id  => 99,
	   -depth    => 10,
	   );

my @forward_filtered_ests;
my @reverse_filtered_ests;

if ( @forward_trans ){
    print STDERR "forward: ".scalar( @forward_trans )."\n";
    @forward_filtered_ests = $est_filter->filter(\@forward_trans);
    push ( @all_trans,  @forward_filtered_ests);
}

if ( @reverse_trans ){
    print STDERR "reverse: ".scalar( @reverse_trans )."\n";
    @reverse_filtered_ests = $est_filter->filter(\@reverse_trans);
    push ( @all_trans, @reverse_filtered_ests );
}

print_GFF( "output_est.gff", \@all_trans, "output");

############################################################

sub print_GFF{
    my ( $file_name, $trans, $type ) = @_;
    open( OUT, ">$file_name") or die("cannot open file $file_name");  
    
    my $exon_count = 0;
    my $trans_count = 0;
    foreach my $tran ( @$trans ){
	$trans_count++;
	my ($evidence_id,$score) = get_evidence_score($tran);
	my @exons = @{$tran->get_all_Exons};
	#print STDERR "transcript ".$tran->dbID." exons ".scalar(@exons)."\n";
	foreach my $exon ( @exons ){
	    $exon_count++;
	    my $strand = "+";
	    if ($exon->strand == -1) {
		$strand = "-";
	    }
	    my $phase = ".";
	    my $g_type = 'input';
	    print OUT $exon_count."\t".
		$type."\t".
		    "exon\t". 
			($exon->start) . "\t" . 
			    ($exon->end) . "\t" .
				$score. "\t" . 
				    $strand . "\t" . 
					$phase . "\t" . 
					    $evidence_id. "\n";
	}
    }
    close (OUT);
    return 1;
}

############################################################

sub get_evidence_score{
    my $t = shift;
    my @exons = @{$t->get_all_Exons};
    my $evidence;
    my $score;
  EXON:
    foreach my $exon ( @exons ){
	foreach my $evi (  @{$exon->get_all_supporting_features} ){
	    $evidence = $evi->hseqname;
	    $score    = $evi->score;
	    last EXON if ( defined $evidence && defined $score );
	}
    }
    return ($evidence, $score);
}

############################################################






