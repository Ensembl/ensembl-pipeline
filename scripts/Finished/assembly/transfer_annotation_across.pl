#!/usr/bin/env perl

=head1 NAME

transfer_annotation_across.pl - transfer annotation from one sequence set on the reference database 
to the same sequence set (same name and coordinate system) on the alternative database. This script 
doesn't need any mapping, it just loops through the list of genes in the reference database and 
save them into the alternative database. Print a warning and skip any gene that is already in 
the alternative database (this is based on the gene's locus name, start, end, strand, biotype and status).

=head1 SYNOPSIS

transfer_annotation_across.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  --assembly=ASSEMBLY                 assembly version ASSEMBLY
  --chromosomes, --chr=LIST           only process LIST chromosomes
  --altdbname                         alternative database name NAME
  --althost                           alternative database host HOST
  --altport                           alternative database port PORT
  --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY
  --altchromosomes, --altchr=LIST     supply alternative chromosome names (the two lists must agree)

Optional arguments:
  -n, --dry_run, --dry=0|1            don't write results to database (default: true)
  --start                             start coordinate on reference and alternative chromosomes
  --end                               end coordinate on reference and alternative chromosomes
  --author                            author login used to lock the assemblies (must be set to write the changes)
  --prefix                            add a prefix to the transferred gene name
  --skip_stable_id                    Gene stable id to skip  

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION


=head1 AUTHOR

ml6@sanger.ac.uk

=head1 CONTACT

ml6@sanger.ac.uk

=cut

use strict;
use warnings;
use List::Util qw[min max];
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);
BEGIN {
    $SERVERROOT = "$Bin/../../../..";
    unshift(@INC, "$Bin");
}

use Getopt::Long;
use Pod::Usage;
use Sys::Hostname;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::Vega::DBSQL::DBAdaptor;
use Bio::Vega::ContigLockBroker;
use Bio::Vega::Author;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

$support->param('verbose', 1);
$support->param('interactive', 0);
$support->param('dry_run', 1);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'althost=s', 
    'altport=s',
    'altdbname=s',
    'altassembly=s',
    'altchromosomes|altchr=s@',
    'assembly=s',         
    'chromosomes|chr=s@',
    'start=i', 
    'end=i',
    'author=s',
    'email=s',
    'skip_stable_id=s@',
    'prefix=s',
);
$support->allowed_params( $support->get_common_params, 'assembly',
    'altassembly', 'chromosomes', 'altchromosomes', );

if ( $support->param('help') or $support->error ) {
    print STDOUT $support->error if $support->error;
    pod2usage(1);
}



$support->comma_to_list( 'chromosomes', 'altchromosomes' );

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params( 'assembly', 'altassembly' );

for my $prm (qw(host port user pass dbname)) {
    $support->param( "alt$prm", $support->param($prm) )
      unless ( $support->param("alt$prm") );
}

# database connection
my $R_dba = $support->get_database('ensembl');
my $A_dba = $support->get_database( 'ensembl', 'alt' );
bless $R_dba, "Bio::Vega::DBSQL::DBAdaptor";
bless $A_dba, "Bio::Vega::DBSQL::DBAdaptor";
my $R_dbh = $R_dba->dbc->db_handle;
my $A_dbh = $A_dba->dbc->db_handle;

my $ref_dbname    = $support->param('dbname');
my $alt_dbname = $support->param('altdbname');
my $assembly    = $support->param('assembly');
my $altassembly = $support->param('altassembly');
my $write_db    = not $support->param('dry_run');
my $start = $support->param('start') || undef;
my $end = $support->param('end') || undef;
my $author      = $support->param('author');
my $email       = $support->param('email') || $author;
my $prefix       = $support->param('prefix');

my $skip_gene = join(',',$support->param('skip_stable_id'));
$skip_gene ||= "NOGENETOSKIP";

throw("must set author name to lock the assemblies if you want to write the changes") if (!$author && $write_db);


my $geneAd   = $R_dba->get_GeneAdaptor;
my $sliceAd  = $R_dba->get_SliceAdaptor;
my $sfeat_Ad = $R_dba->get_SimpleFeatureAdaptor();

my $geneAd_alt   = $A_dba->get_GeneAdaptor;
my $sliceAd_alt  = $A_dba->get_SliceAdaptor;
my $sfeat_Ad_alt = $A_dba->get_SimpleFeatureAdaptor();

my @R_chr_list = $support->param('chromosomes');
die "Reference chromosome list must be defined!" unless scalar(@R_chr_list);

my @A_chr_list = $support->param('altchromosomes');
if ( scalar(@R_chr_list) != scalar(@A_chr_list) ) {
    die "Chromosome lists do not match by length";
}

SET: for my $i ( 0 .. scalar(@R_chr_list) - 1 ) {
    my $R_chr = $R_chr_list[$i];
    my $A_chr = $A_chr_list[$i];


    # for stats
    my $total_genes      = 0;
    my $transfered_genes = 0;
    my $total_sf         = 0;
    my $total_created_genes  = 0;
    my $transfered_sf    = 0;
    my $skipped_sf       = 0;
    my $skipped_gene     = 0;
    my $missed_g        = 0;
    my $locked = 0;

    my $ref_sl =
        $sliceAd->fetch_by_region( 'chromosome', $R_chr, $start, $end, undef,$assembly );
    my $alt_sl =
        $sliceAd_alt->fetch_by_region( 'chromosome', $A_chr, $start, $end, undef,$altassembly );


    # Lock the reference and alternative assemblies
    my ($cb,$author_obj);
    if($write_db){
        eval {
            $cb = Bio::Vega::ContigLockBroker->new(-hostname => hostname);
            $support->log_verbose("Locking $ref_dbname:$R_chr and $alt_dbname:$A_chr\n");
            $author_obj = Bio::Vega::Author->new(-name => $author, -email => $email);
            $cb->lock_clones_by_slice($ref_sl,$author_obj,$R_dba);
            $locked = 1;
            $author_obj = Bio::Vega::Author->new(-name => $author, -email => $email);
            $cb->lock_clones_by_slice($alt_sl,$author_obj,$A_dba);
        };
        if($@){
            warning("Cannot lock assemblies $ref_dbname:$R_chr and $alt_dbname:$A_chr with author name $author\n$@\n");
            $cb->remove_by_slice($ref_sl,$author_obj,$R_dba) if $locked;
            next SET;
        }
    }

    $A_dbh->begin_work;

    eval {
        $support->log_verbose("Annotation transfer $ref_dbname:$R_chr => $alt_dbname:$A_chr...\n");


        my @genes;
        @genes       = @{ $ref_sl->get_all_Genes };
        $total_genes = scalar @genes;
        @genes = sort { $a->start <=> $b->start } @genes;

        # transfer simple features (PolyA_site/_signal)
        my @proj_feat;
        my @simple_features = @{ $sfeat_Ad->fetch_all_by_Slice($ref_sl) };
        $total_sf = scalar @simple_features;
        SF: foreach my $f (@simple_features) {
            $f->dbID(undef);
            $f->adaptor(undef);
            $f->slice->{'coord_system'} = $alt_sl->coord_system();
            $f->slice->adaptor($sliceAd_alt);
            my $existing_sf = $sfeat_Ad_alt->fetch_all_by_Slice( 
                $f->feature_Slice, 
                $f->analysis->logic_name );
            if ( ! @$existing_sf ) {
                push @proj_feat, $f;
                $transfered_sf++;
            } else {
            	$skipped_sf++;
            }
        }
        
      GENE: foreach my $g (@genes) {
            next if $g->stable_id =~ /$skip_gene/;
            my ($gene_name_attrib) = @{ $g->get_all_Attributes('name') };
            my $gene_name = $gene_name_attrib ? $gene_name_attrib->value : "UNDEF";
            my $stable_id = $g->stable_id;
            &remove_all_db_ids($g);
            $g->slice->{'coord_system'} = $alt_sl->coord_system();
            $g->slice->adaptor($sliceAd_alt);
            
            my $existing_gene;
            
            eval {
                $existing_gene = $geneAd_alt->fetch_all_by_Slice( 
                    $g->feature_Slice,
                    $g->analysis->logic_name );
			};
			if($@) {
			    $support->log_verbose($@);
			    next GENE;
			}
			my $exist = 0;
			my $tg_key = join(":",$gene_name,$g->seq_region_start,$g->seq_region_end,$g->biotype,$g->status,$g->source);
			foreach (@$existing_gene) {
				my ($existing_name_att) = @{ $_->get_all_Attributes('name') };
                my $existing_gene_name = $existing_name_att ? $existing_name_att->value : "UNDEF";
			    my $existing_gkey = join(":",$existing_gene_name,$_->seq_region_start,$_->seq_region_end,$_->biotype,$_->status,$_->source);
			    $exist = 1 if $tg_key eq $existing_gkey;
			}
			
			if($exist) {
			    $support->log_verbose(
			        sprintf(
			            "WARNING: SKIP GENE %s %s ($gene_name) already transfered\n",
			            $stable_id, $g->biotype, $R_chr,
			        ));
			    $skipped_gene++; 
			    next GENE;
			} 

            $support->log_verbose(
                sprintf(
                    "Copying gene %s %s %d %d\n",
                    $stable_id, $g->biotype, $g->start, $g->end
                )
            );


            if ( $geneAd_alt->store($g) ) {
                $transfered_genes++;
                $support->log_verbose(
                sprintf(
                    "GENE %s %s successfully COPIED (%s:%d-%d)\n",
                    $g->stable_id, $g->biotype, $R_chr,
                    $g->start,      $g->end)
                );
                
            } else {
                throw(
                    sprintf(
                        "GENE %s %s cannot be copied across (%s:%d-%d)",
                        $g->stable_id, $g->biotype, $R_chr,
                        $g->start,      $g->end
                     )
                );
            }
        }
                
        # save polyA features if gene transfered
        $sfeat_Ad_alt->store(@proj_feat) if ( @proj_feat && $transfered_genes );

        $write_db ? $A_dbh->commit : $A_dbh->rollback;
    };

    if ($@) {
        $A_dbh->rollback;
        $support->log_verbose(
                "UNABLE TO TRANSFER ANNOTATION FROM $ref_dbname:$R_chr:$assembly to $alt_dbname:$A_chr:$altassembly\n[$@]\n" );
    } else {
        # print annotation transfer stats
        $support->log_verbose(
            sprintf(
"INFO: Annotation transfer %s:%s:%s => %s:%s:%s
INFO: transfered genes: %d/%d
INFO: skipped genes: %d/%d
INFO: transfered PolyA features: %d/%d
INFO: skipped PolyA features: %d/%d\n",
                $ref_dbname, $R_chr, $assembly,
                $alt_dbname, $A_chr, $altassembly,
                ($transfered_genes - $total_created_genes), $total_genes,
                $skipped_gene,$total_genes,
                $transfered_sf,    $total_sf,
                $skipped_sf,       $total_sf
            )
        );
    }
    
    $ref_sl =
        $sliceAd->fetch_by_region( 'chromosome', $R_chr, $start, $end, undef,$assembly );
    $alt_sl =
        $sliceAd_alt->fetch_by_region( 'chromosome', $A_chr, $start, $end, undef,$altassembly );
    
    if($write_db){
    # remove the assemblies locks
        eval {
            $support->log_verbose("Removing $ref_dbname:$R_chr and $alt_dbname:$A_chr Locks\n");
            $cb->remove_by_slice($ref_sl,$author_obj,$R_dba);
            $cb->remove_by_slice($alt_sl,$author_obj,$A_dba);
        };
        if($@){
            warning("Cannot remove locks from assemblies $ref_dbname:$R_chr and $alt_dbname:$A_chr with author name $author\n$@\n");
        }
    }
}


sub remove_all_db_ids {
    my $g = shift;
    
    $g->dbID(undef);
    $g->adaptor(undef);
    $g->stable_id(undef);
    
    foreach my $t (@{$g->get_all_Transcripts;}){
	    $t->dbID(undef);
	    $t->stable_id(undef);
	    $t->adaptor(undef);
	    if ( $t->translation ) {
	        $t->translation->dbID(undef);
	        $t->translation->adaptor(undef);
	        $t->translation->stable_id(undef);
	    }
	    foreach my $e ( @{ $t->get_all_Exons } ) {
	        $e->dbID(undef);
	        $e->adaptor(undef);
	        $e->stable_id(undef);
	    }
    }
}
