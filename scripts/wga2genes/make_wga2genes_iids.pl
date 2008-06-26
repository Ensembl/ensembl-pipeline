#!/usr/local/ensembl/bin/perl

# Generation of input ids for WGA2Genes. Segments the
# genome into regions that (1) contain at least one gene,
# and (2) do not interrupt a chain.

use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::Analysis;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $query_dbname,
    $query_dbhost,
    $query_dbport,
    $query_dbuser,
    $query_dbpass,
    $target_dbname,
    $target_dbhost,
    $target_dbport,
    $target_dbuser,
    $target_dbpass,
    $compara_dbname,
    $compara_dbhost,
    $compara_dbport,
    $compara_dbuser,
    $compara_dbpass,
    $seq_region_name,
    $source_type,
    $source_align_type,
    $kill_list,
    $logic_name,
    $write,
    $verbose,
);

$dbuser = 'ensro';
$dbport = 3306;
$query_dbuser = 'ensro';
$query_dbport = 3306;
$target_dbuser = 'ensro';
$target_dbport = 3306;
$compara_dbuser = "ensro";
$compara_dbport = 3306;

$source_type = 'protein_coding';
$source_align_type = 'BLASTZ_NET';

&GetOptions(
            
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'querydbname=s' => \$query_dbname,
            'querydbhost=s' => \$query_dbhost,
            'querydbport=s' => \$query_dbport,
            'querydbuser=s' => \$query_dbuser,
            'querybpass=s' => \$query_dbpass,
            'targetdbname=s' => \$target_dbname,
            'targetdbhost=s' => \$target_dbhost,
            'targetdbport=s' => \$target_dbport,
            'targetdbuser=s' => \$target_dbuser,
            'targetbpass=s'  => \$target_dbpass,
            'comparadbname=s' => \$compara_dbname,
            'comparadbhost=s' => \$compara_dbhost,
            'comparadbport=s' => \$compara_dbport,
            'comparadbuser=s' => \$compara_dbuser,
            'comparadbpass=s' => \$compara_dbpass,
            'logic=s' => \$logic_name,
            'aligntype=s'   => \$source_align_type,
            'genetype=s'    => \$source_type,
            'kill=s'        => \$kill_list,
            'seq_region_name=s' => \$seq_region_name,
            'write' => \$write,
            'verbose' => \$verbose,
            );


my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
        '-pass' => $dbpass,
);

my $ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
if (not defined $ana) {
  warn "Could not find analysis with logic name '$logic_name' in pipe db; creating\n";
  $ana = Bio::EnsEMBL::Pipeline::Analysis->new;
  $ana->input_id_type("SLICE");
}



my ($query_db, $target_db, $compara_db);

$query_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
	'-dbname' => $query_dbname,
	'-host' => $query_dbhost,
	'-port' => $query_dbport,
	'-user' => $query_dbuser,
        '-pass' => $query_dbpass,
        );

if (not defined $target_dbname) {
  $target_db = $db;
} else {
 $target_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
	'-dbname' => $target_dbname,
	'-host' => $target_dbhost,
	'-port' => $target_dbport,
	'-user' => $target_dbuser,
        '-pass' => $target_dbpass,
        );

}

if (defined $compara_dbname) {
  $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->
      new(
          '-dbname' => $compara_dbname,
          '-host' => $compara_dbhost,
          '-port' => $compara_dbport,
          '-user' => $compara_dbuser,
          '-pass' => $compara_dbpass,
          );
}

my $q_species =
    $query_db->get_MetaContainerAdaptor->get_Species->binomial;
my $t_species =
    $target_db->get_MetaContainerAdaptor->get_Species->binomial;

if (defined $kill_list) {
  my $list = {};

  open LIST, $kill_list or die "Could not open '$kill_list' for reading\n";
  while(<LIST>) {
    /^\#/ and next;
    /^(\S+)/ and $list->{$1} = 1;
  }

  $kill_list = $list;

  $verbose and printf "Read KILL list (%d items)\n", scalar(keys %$kill_list);
} else {
  $kill_list = {};
}

my @tl_slices = @{$query_db->get_SliceAdaptor->fetch_all('toplevel')};


my @slices;

if (defined $seq_region_name) {
  my $sl = $query_db->get_SliceAdaptor->fetch_by_name("toplevel::${seq_region_name}:::");
  foreach my $tl_sl (@tl_slices) {
    if ($sl->coord_system->name eq $tl_sl->coord_system->name and
        $sl->seq_region_name eq $tl_sl->seq_region_name) {
      push @slices, $tl_sl
    }
  }
} else {
  @slices = @tl_slices;
}


my ($gdb_adap, $q_gdb, $t_gdb, $mlss, $gaba);

if (defined $compara_db) {
  $gdb_adap = $compara_db->get_GenomeDBAdaptor;
  $q_gdb = $gdb_adap->fetch_by_name_assembly($q_species);
  $t_gdb = $gdb_adap->fetch_by_name_assembly($t_species);

  $mlss = $compara_db->get_MethodLinkSpeciesSetAdaptor
      ->fetch_by_method_link_type_GenomeDBs($source_align_type,
                                          [$q_gdb, $t_gdb]);
  $gaba = $compara_db->get_GenomicAlignBlockAdaptor;
}

foreach my $sl (@slices) { 
  my @genes = @{$sl->get_all_Genes};
  @genes = grep { $_->biotype eq $source_type } @genes;

  next if not @genes;
  
  my @gene_regions;

  foreach my $g (@genes) {
    next if exists $kill_list->{$g->stable_id};
    
    my ($min_start, $max_end);
    foreach my $t (@{$g->get_all_Transcripts}) {
      next if exists $kill_list->{$t->stable_id};

      foreach my $e (@{$t->get_all_translateable_Exons}) {
        if (not defined $min_start or $e->start < $min_start) {
          $min_start = $e->start;
        }
        if (not defined $max_end or $e->end > $max_end) {
          $max_end = $e->end;
        }
      }
    }

    if (defined $min_start and defined $max_end) {
      push @gene_regions, {
        qstart => $min_start + $sl->start - 1,
        qend   => $max_end   + $sl->start - 1,
      };
    }
  }

  @gene_regions = sort { $a->{qstart} <=> $b->{qstart} } @gene_regions;


  my (%chains, @chains, @ov_chains);

  if (defined $compara_db) {
    my $dnafrag = $compara_db->get_DnaFragAdaptor->
        fetch_by_GenomeDB_and_name($q_gdb,
                                   $sl->seq_region_name);
    
    for(my $gen_start = $sl->start; $gen_start < $sl->end; $gen_start += 1000000) {
      my $gen_end = $gen_start + 1000000 - 1;
      $gen_end = $sl->end if $gen_end > $sl->end;
      
      $verbose and printf(STDERR "Fetching blocks for $source_align_type %s %d %d...\n", 
                          $dnafrag->name, 
                          $gen_start,
                          $gen_end);
      
      my $gen_al_blocks =
          $gaba->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,
                                                           $dnafrag,
                                                           $gen_start,
                                                           $gen_end);
      
      $verbose and printf(STDERR "Extracting chain information...\n");
      
      foreach my $block (@$gen_al_blocks) {
        my $qga = $block->reference_genomic_align;
        my ($tga) = @{$block->get_all_non_reference_genomic_aligns};
        
#        my $chain_id =  $qga->genomic_align_group_id_by_type("chain");
        my $chain_id =  $block->group_id;
        print "Group id = " . $block->group_id . "\n";
        next if !$block->group_id;
        
        if ($block->reference_genomic_align->dnafrag_strand < 0) {
          $block->reverse_complement;
        }
        
        if (not exists $chains{$chain_id}) {
          $chains{$chain_id} = {
            qstart => $qga->dnafrag_start,
            qend   => $qga->dnafrag_end,
          };
        } else {
          if ($qga->dnafrag_start < $chains{$chain_id}->{qstart}) {
            $chains{$chain_id}->{qstart} = $qga->dnafrag_start;
          }
          if ($qga->dnafrag_end > $chains{$chain_id}->{qend}) {
            $chains{$chain_id}->{qend} = $qga->dnafrag_end;
          }
        }
      }
    }
  }

  map { push @chains, $chains{$_} } keys %chains;
  @chains = sort { $a->{qstart} <=> $b->{qstart} } @chains;
  
  # remove chains that do not overlap with genes
  my @gchains;
  foreach my $ch (@chains) {
    my $ov = 0;

    foreach my $gr (@gene_regions) {
      if ($gr->{qstart} > $ch->{qend}) {
        last;
      } elsif ($gr->{qend} >= $ch->{qstart}){ 
        $ov = 1;
        last;
      }
    }
    
    if ($ov) {
      push @gchains, $ch;
    }
  }

  my @sl_regs;

  # merge ov_chains and gene clusters
  foreach my $cm (sort { $a->{qstart} <=> $b->{qstart} } (@gchains, @gene_regions)) {
    if (not @sl_regs or $sl_regs[-1]->{qend} < $cm->{qstart}) {
      push @sl_regs, $cm;
    } elsif ($cm->{qend} > $sl_regs[-1]->{qend}) {
      $sl_regs[-1]->{qend} = $cm->{qend};
    }
  }
  
  foreach my $cl (@sl_regs) {
    #printf("REGION %s %d %d %d\n",
    #      $sl->seq_region_name,
    #      $cl->{qstart},
    #      $cl->{qend},
    #      $cl->{qend} - $cl->{qstart} + 1);
    #next;
    
    my $iid = sprintf("%s::%s:%d:%d:", 
                      $sl->coord_system->name,
                      $sl->seq_region_name,
                      $cl->{qstart},
                      $cl->{qend});
    
    if ($write) {
      my $s_inf_cont = $db->get_StateInfoContainer;
      
      eval {
         print $iid . "\n";
        $db->get_StateInfoContainer->
            store_input_id_analysis( $iid, $ana, '' );
      };
      if ($@) {
        print STDERR "Input id $iid already present\n";
      } else {
        print STDERR "Stored input id $iid\n";
      }       
    } else {
      printf("INSERT into input_id_analysis(input_id, input_id_type, analysis_id) " . 
             " values('%s', '%s', %d);\n", 
             $iid, 
             $ana->input_id_type,
             $ana->dbID);
      
    }
  }
}


