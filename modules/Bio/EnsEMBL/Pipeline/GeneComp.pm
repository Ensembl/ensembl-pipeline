
#
# BioPerl module for Bio::EnsEMBL::Pipeline::GeneComp
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneComp - Comparison of old and new Gene/Transcript/Exons

=head1 SYNOPSIS

    $vc is a Bio::EnsEMBL::Virtual::Contig (it has to be attached to a crossmatch
    database, in order to support the get_old_Exons and get_old_Genes calls)
    
    $arcdb is the archive database to which all dead genes, transcripts and
    translations will be written

    $finaldb is a new database object to which all new, mapped genes will be
    written (with all their transcripts, exons, etc.)
    
    $log is the logfile filehandle to which all mappings are logged

    my $gc =  Bio::EnsEMBL::Pipeline::GeneComp->new(-vc => $vc,
						    -archive => $arcdb,
						    -finaldb => $finaldb,
						    -log => \*LOG);

    $gc->map();

=head1 DESCRIPTION

This object deals with mapping exons, transcript and genes from old 
versions through to new versions. To do the mapping we need to get 
the old exons out in the new coordinates (remapping). This currently 
is hidden behind the method call get_old_Exons() on a contig object. 
This call returns old exon objects in the new coordinates, by getting
crossmatch feature pairs from a crossmatch database.

All in all this object holds on to 5 databases. It holds on to:

1-Temp (new) database ($vc->dbobj)
2-Old database ($vc->dbobj->_crossdb->old_dbobj->)
3-Crossmatch database ($vc->dbobj->_crossdb)
4-Archive database ($self->archive)
5-Final database ($self->finaldb)

MAPPING LOGIC:

This module is complex. I dont think there is anyway around
this. There are two basic pieces of logic - rules for exon migration
and rules for gene/transcript migration.

For exon migration, if the start/end/strand in new coordinates of the
exons are the same then it gets the old exon id. If the dna sequence
has changed, this increments the version number.  If not it stays the
same.

For gene/transcript migration, the following things happen. Old and
New genes are clustered into 4 sets on the basis of shared exons (this
occurs after exon mapping, done outside of this module)

   Simple - one old gene, one new gene
   Split  - one old gene, >1 new genes
   Merges - >1 new genes, one old gene
   Unassigned new genes 

There is the possibility of >1 old gene with >1 new gene. Depending on
the order of discovery, this will be classified as a split or
merge. This is a known bug/feature.

For each cluster, old transcripts are sorted by length and then fitted
to new transcripts, with the best fit taking a win (fit on the number
of co-linear identical id''d exons). Perfect matches (all exons the
same id) trigger a direct assignment.

Versioning for transcripts is that any addition/removal of an exon, or
any update in sequence of an exon rolls up the transcript version. The
gene version clicks up on any transcript version or any transcript
addition/deletion.



=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::GeneComp;
use vars ('@ISA');
use Bio::Root::RootI;
use strict;
use Carp;

@ISA = ('Bio::Root::RootI');


=head2 new

 Title   : new
 Usage   : $genecomp = Bio::EnsEMBL::Pipeline::GeneComp->new(-vc => $vc, -archive => $archive,-finaldb => $db,-log => $log);
 Function: Builds a new genecomp object, based around a virtual contig, 
 Example :
 Returns : a new genecomp object
 Args    : 


=cut

sub new {
  my($class,@args) = @_;
  my $self = {};
  bless ($self,$class);
  my ($vc,$archive,$log,$map) = 
      $self->_rearrange([qw(VC
			    ARCHIVE
			    LOG
			    MAP
			    )],@args);

  #REmoved finaldb, now writing....

  if( !defined $vc || !ref $vc || !$vc->isa('Bio::EnsEMBL::Virtual::Contig') ) {
      $self->throw("must have a virtual contig, got a [$vc]");
  }

  if( !defined $archive || !ref $archive || !$archive->isa('Bio::EnsEMBL::DBArchive::Obj') ) {
      $self->throw("must have a archive, got a [$archive]");
  }

  #if( !defined $finaldb || !ref $finaldb || !$finaldb->isa('Bio::EnsEMBL::DBSQL::Obj') ) {
  #    $self->throw("must have a final database, got a [$finaldb]");
  #}

  if( !defined $log ) {
      $self->throw("Must have a log file");
  }

  $self->vc($vc);
  $self->archive($archive);
  #$self->finaldb($finaldb);
  $self->log($log);
  $self->mapfile($map);
  $self->{_fitted_trans_hash}=();

  # this usually dies if we can't write
  print $log "Built GeneComp object\n";

  return $self;
}


=head2 map

 Title   : map
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub map{
    my ($self) = @_;
    my $log = $self->log();
    my $mapfile = $self->mapfile();
    my $vc = $self->vc();
    
    #Map exons
    print $log "Mapping exons:\n";
    my ($before_e,$after_e)=$self->map_temp_Exons_to_real_Exons();
    if ($before_e != $after_e) {
	print $log "EXON MAPPING BUG: Size of temp exons and final exons not equal, bad bug! Temp exons: $before_e, final: $after_e";
    }
    #Map genes
    print $log "Mapping genes:\n";
    my ($before_g,$finalgenes) = $self->map_temp_Genes_to_real_Genes();
    my $after_g = scalar (@$finalgenes);
    if ($before_g != $after_g) {
	print $log "GENE MAPPING BUG: Size of temp genes and final genes not equal, bad bug! Temp genes: $before_g, final: $after_g";
    }

    #Archiving is dealt with within the mapping methods
    #We archive all dead gene ids, and infor on dead transcripts
    #and proteins

    #convert back to raw coordinates
    #my @rawgenes;
    #print $log "Converting genes to raw coordinates...\n";
    foreach my $gene (@$finalgenes) {
	#print $log "Dumping structure of genes\n";
	$gene->_dump(\*$mapfile);
	#foreach my $exon ($gene->all_Exon_objects) {
	#    $exon->contig_id($vc->id);
	#}
	#my $newgene;
	#eval {
	#    $newgene = $vc->convert_Gene_to_raw_contig($gene);
	#};
	#if ($@) {
	#    print $log "CONVERSION BUG: Could not convert gene ".$gene->id." to raw contigs coordinates because of $@\n";
	#}
	#else {
	#    push (@rawgenes,$newgene);
	#}
    }
    #print $log "Writing genes to final database...\n";
    #write final set to final db
    #foreach my $rgene (@rawgenes) {
	#eval {
	#    $self->finaldb->gene_Obj->write($rgene);
	#};
	#if ($@) {
	#    print $log "GENE WRITING BUG: Could not write back gene ".$rgene->id." to database,error:\n$@";
	#}
    #}
    return 1;
}

=head2 map_temp_Exons_to_real_Exons

 Title   : map_temp_Exons_to_real_Exons
 Usage   : $gc->map_temp_Exons_to_real_Exons
 Function: mapping of temp exons to old exons
 Example : 
 Returns : nothing
 Args    : none

=cut

sub map_temp_Exons_to_real_Exons{
   my ($self) = @_;

   my $vc = $self->vc();
   my $log = $self->log;
   print $log "Getting all new exons...";
   my @tempexons=$vc->get_all_Exons();
   my $t_size=scalar(@tempexons);
   print $log " got $t_size new exons\n";
   
   if( scalar(@tempexons) == 0 ) {
       print $log "No new exons found on this contig. Nothing to map - returning an empty list";
       return ();
   }

   # ok - we are ready to rock and roll...

   # we need some internal data structures during the mapping process.

   my %contig;   # hash of contig objects
   my %oldexons; # a hash of arrays of old exon objects
   my %moved;    # shows which old exons we have moved (or not).
   my %oldexonhash; # direct hash of old exons for final dead call
   my %ismapped; # tells us whether temporary exons have been mapped or not
   my %temp_old;
   my @new;      # new exons
   my @mapped;   #mapped exons
  
   # initialise ismapped to 0 for all temp exons and
   # build a hash of contig objects we want to see and exons.
   # better to do this is a separate loop as we have many exons to one contig

   foreach my $tempexon ( @tempexons ) {
       if( exists $ismapped{$tempexon->id} ) {
	   $self->throw("you have given me temp exons with identical ids. bad...");
       }

       $ismapped{$tempexon->id} = 0;
   }
   print $log "Getting all old exons...\n";
   my @oldexons=$vc->get_old_Exons($log);
   my $size=scalar(@oldexons);
   print $log "got $size old exons mapped\n";
   # get out one date for this mapping....

   my $time = time();

   @tempexons  = sort { $a->start <=> $b->start } @tempexons;
   @oldexons   = sort { $a->end   <=> $b->end } @oldexons;
   # go over each exon and map old->new...
   my $s1=scalar(@tempexons);
   my $s2=scalar(@oldexons);
   TEMPEXON :

   foreach my $tempexon ( @tempexons ) {
       foreach my $oldexon ( @oldexons ) {
	   if ($tempexon->start > $oldexon->end) {
	       #Skip this old exon, it ends before temp starts
	       next;
	   }
	   if( $oldexon->start > $tempexon->end ) {
	       #Exiting comparison loop because this old exon 
	       #starts after end of temp
	       last;
	   }

	   # if start/end/strand is identical, it is definitely up for moving.
	   
	   if( $tempexon->start == $oldexon->start &&
	       $tempexon->end   == $oldexon->end   &&
	       $tempexon->strand == $oldexon->strand ) {

	       # ok - if we have already moved this - ERROR

	       if( $moved{$oldexon->id} == 1 ) {
		   print $log "attempting to move old exon twice ".$oldexon->id." with identical start/end/strand. ".$tempexon->start.":".$tempexon->end.":".$tempexon->strand." Not clever!\n";
		   next;
	       }

	       # set ismapped to 1 for this tempexon

	       $ismapped{$tempexon->id} = 1;
	       $temp_old{$tempexon->id}=$oldexon->id;
	       print $log "EXON MAP IDENTICAL: ",$tempexon->id," ",$oldexon->id," "; # will add version
	       # we move the id. Do we move the version?
	       $tempexon->id($oldexon->id);
	       push (@mapped,$tempexon);
	       
	       if( $oldexon->seq eq $tempexon->seq) {
		   # version and id
		   $tempexon->version($oldexon->version);
		   print $log "Identical\n";
		   # don't update modified
	       } else {
		   my $v=$oldexon->version()+1;
		   $tempexon->version($oldexon->version()+1);
		   print $log $v."\n";
		   $tempexon->modified($time);
	       }
	       
	       $moved{$oldexon->id} = 1;

	       # remove old exons we never need to look at again
	       # because they will never overlap with tempexons
	       
	       while( my $tempoldexon = shift @oldexons ) {
		   if( $tempoldexon->end >= $tempexon->start ) {
		       unshift(@oldexons,$tempoldexon);
		       last;
		   }
	       }

	       next TEMPEXON;
	   }
	   
       }
   



      # we can't map this exon directly. Second scan over oldexons to see whether
      # there is a significant overlap

      my $biggestoverlap = undef;
      my $overlapsize = 0;
      foreach my $oldexon ( @oldexons ) {
	  if ($tempexon->start > $oldexon->end) {
	      #Skip this old exon, it ends before temp starts
	      next;
	  }
	  if( $oldexon->start > $tempexon->end ) {
	      #Go out of loop, reached old exon after temp
	      last;
	  }

	  
	  if( $oldexon->overlaps($tempexon) && $moved{$oldexon->id} == 0 ) {
	      my ($tstart,$tend,$tstrand) = $oldexon->intersection($tempexon);
	      if( !defined $biggestoverlap ) {
		  $biggestoverlap = $oldexon;
		  $overlapsize = ($tend - $tstart +1); 
              } else {
		  if( ($tend - $tstart +1) > $overlapsize ) {
		      $biggestoverlap = $oldexon;
		      $overlapsize = ($tend - $tstart +1);
		  }
	      }
	  }
      }


       # remove old exons we never need to look at again
       # because they will never overlap with tempexons
       
       #while( my $tempoldexon = shift @oldexons ) {
#	   print STDERR "Shifting out ".$tempoldexon->id."\n";
	#   if( $tempoldexon->end >= $tempexon->start ) {
	#       print STDERR "Putting back in ".$tempoldexon->id."\n";
	#       unshift(@oldexons,$tempoldexon);
	#       last;
	#   }
       #}
       

       # if I have got a biggest overlap - map across...
       if( defined $biggestoverlap ) {
	   # set ismapped to 1 for this tempexon
	   
	   $ismapped{$tempexon->id} = 1;
	   
	   # we move the id. Do we move the version?
	   print $log "EXON MAP BEST FIT: ",$tempexon->id," ",$biggestoverlap->id," ",$biggestoverlap->version()+1,"\n";

	   $temp_old{$tempexon->id}=$biggestoverlap->id;
	   $tempexon->id($biggestoverlap->id);
	   $tempexon->version($biggestoverlap->version()+1);
	   $tempexon->modified($time);
	   push(@mapped,$tempexon);
	   $moved{$biggestoverlap->id} = 1;
	   my $size=scalar (@oldexons);
	   while( my $tempoldexon = shift @oldexons ) {
	       if( $tempoldexon->end >= $tempexon->start ) {
		   unshift(@oldexons,$tempoldexon);
		   last;
	       }
	   }

	   next TEMPEXON;
       } else {
	   # ok - new Exon
	   my $tempid = $tempexon->id();
	   $tempexon->id($self->archive->get_new_stable_ids('exon',1));
	   $tempexon->created($time);
	   $tempexon->modified($time);
	   $tempexon->version(1);
	   $temp_old{$tempid}=$tempexon->id;
	   print $log "EXON MAP NEW: ",$tempid," ",$tempexon->id," 1\n"; 
	   push(@new,$tempexon);
	   next TEMPEXON;
       }
    
       $self->throw("Error - should never reach here!");
   }
	      
   $self->{'_exon_map_hash'} = \%temp_old;
   $self->{'_mapped_exons'} = \@mapped;
   $self->{'_new_exons'} = \@new;
   my $final_size= scalar (@mapped)+scalar(@new);
   return ($t_size,$final_size);
}


=head2 map_temp_Genes_to_real_Genes

 Title   : map_temp_Genes_to_real_Genes
 Usage   : $gc->map_temp_Genes_to_real_Genes
 Function: mapping of temp genes to old genes
 Example :
 Returns : nothing
 Args    : none

=cut

sub map_temp_Genes_to_real_Genes{
    my $self = shift;

    my $vc= $self->vc();
    my $log = $self->log();
    my $exon_map = $self->{'_exon_map_hash'};
    
    if( !defined $vc) {
	die("Not passing the Virtual Contig to do gene mapping");
    }
    print $log "Getting all new genes... ";
    my @tempgenes = $vc->get_Genes_by_Type('merged');
    my $temp_size = scalar @tempgenes;
    print $log "got $temp_size new genes\n";
    foreach my $gene (@tempgenes) {
	foreach my $exon ($gene->all_Exon_objects) {
	    if ($exon_map->{$exon->id}) {
		$exon->id($exon_map->{$exon->id});
	    }
	}
    }
    print $log "Getting all old genes... ";
    my @oldgenes  = $vc->get_old_Genes();
    my $size=scalar(@oldgenes);
    print $log "got $size old genes\n";

    my @final_genes; #final set of genes with mapped and new ids
    my %olde2t; # hash of exon id to array of transcript id
    my %olde2g; # hash of exon id to gene id
    my %oldt;   # hash of transcript on transcript id (effectively t->e mapping)
    my %oldg;   # hash of gene on gene id (effectively g->e mapping)
    my %newg;   # hash of new genes on geneid
    
    my %has_done_new;  # 1 or 0 depending on whether this gene has be moved or not
    my %has_moved_old; # 1 or 0 depending on whether this gene has mapped forward or not
    my @dead_gene_ids; #Array of dead genes
    my %newe2g; # hash of new exon to new gene objects
    
    # map old exons to gene and transcripts
    
    foreach my $og ( @oldgenes ) {
	
	$oldg{$og->id} = $og;
	$has_moved_old{$og->id} = 0;
	
	foreach my $ot ( $og->each_Transcript ) {
	    $oldt{$ot->id} = $ot;
	    foreach my $oe ( $ot->each_Exon ) {
		if( ! exists $olde2t{$oe->id} ) {
		    $olde2t{$oe->id} = []; 
		}
		push(@{$oldt{$oe->id}},$ot->id);
		$olde2g{$oe->id} = $og->id;
	    }
	}
    }
    
    # build the newe2g hash and newg
    
    foreach my $ng ( @tempgenes ) {
	$newg{$ng->id} = $ng;
	$has_done_new{$ng->id} = 0;
	
	foreach my $ne ( $ng->each_unique_Exon ) {
	    $newe2g{$ne->id} = $ng;
	}
   }
    
    
    # build hashes for the classes of new genes that we have.
    # simple, split, merge
    
    my %simple;  # key is old gene id, value is the new gene id to map to
    my %reversed_simple; # key is the new gene id, value is the old gene id
    
    my %split;   # key is old gene id, value an array of new gene id to map to
    my %merge;   # hash on the new gene id -> value is an array of old gene ids
    
    foreach my $og ( @oldgenes ) {
	
	my $currentgeneid = undef;

	foreach my $oe ( $og->each_unique_Exon ) {
	    if( ! exists $newe2g{$oe->id} ) {
		# this exon does not exist in the new transcripts!
		next;
	    }
	    
	    
	    my $tgeneid = $newe2g{$oe->id}->id;
	    
	    # if the tgeneid is the same as current then we have already 
	    # looked at this case
	    
	    if( defined $currentgeneid && $tgeneid eq $currentgeneid ) {
		next;
	    }
	    
	    $currentgeneid = $tgeneid;
 	    
            # this could be more of a merge or split.
	    if( exists $split{$og->id} ) {
		# then this is more of a split
		push(@{$split{$og->id}},$tgeneid);
		next;
	    }
	    
	    if( exists $merge{$tgeneid} ) {
		# then this is more of a merge
		push(@{$merge{$tgeneid}},$og->id);
		next;
	    }
	    
	    
	    # ok - if reversed_simple has this id, then it is not simple anymore
	    # - becomes a merge.
	    
	    if( exists $reversed_simple{$tgeneid} ) {
	       # remove both sides from simple
	       my $oldid = $reversed_simple{$tgeneid};

	       delete $reversed_simple{$tgeneid};
	       delete $simple{$oldid};

	       # start a merge with oldid and og->id

	       $merge{$tgeneid} = [];
	       push(@{$merge{$tgeneid}},$og->id,$oldid);
	       next;
	   }

	   # if simple has this oldid then it is not simple anymore - 
	   # becomes a split

	   if( exists $simple{$og->id} ) {
	       # remove both sides from simple/reverse
	       my $previous_new = $simple{$og->id};
	       delete $simple{$og->id};
	       delete $reversed_simple{$previous_new};
	       
	       # start a split, on og->id with previous and tgeneid

	       $split{$og->id} = [];
	       push(@{$split{$og->id}},$previous_new,$tgeneid);
	       next;
	   }

	   # it is simple - hurray!
	
	   $simple{$og->id} = $tgeneid;
	   $reversed_simple{$tgeneid} = $og->id;
       }
   }


   # now we have a list of potential mappings. We step over these, mapping transcripts
   # from old to new.

   # simple is easy ;).
    my $now = time();
    foreach my $oldgeneid ( keys %simple ) {
       my $newgeneid = $simple{$oldgeneid};

       # flag that we have done this move before the ids change ;)

       $has_done_new{$newgeneid} = 1;
       $has_moved_old{$oldgeneid} = 1;


       my @newtrans = $newg{$newgeneid}->each_Transcript;
       my @oldtrans = $oldg{$oldgeneid}->each_Transcript;
       my %old_tg;
       foreach my $trans ($oldg{$oldgeneid}->each_Transcript) {
	   $old_tg{$trans->id}=$oldg{$oldgeneid};
       }
       
       my $should_increment = $self->map_temp_Transcripts_to_real_Transcripts($oldg{$oldgeneid},\@newtrans,\@oldtrans,%old_tg);


       # deal with the mapping of ids
       print $log "GENE MAP: MOVED ",$newgeneid," ",$oldgeneid,"\n";
       $self->{_done_gene_hash}->{$newgeneid}=1;
       $self->{_done_gene_hash}->{$oldgeneid}=1;
       $newg{$newgeneid}->id($oldgeneid);
       
       if( $should_increment ) {
	   $newg{$newgeneid}->version($oldg{$oldgeneid}->version+1);
	   $newg{$newgeneid}->modified($now);
       }
       #print $log "About to dump after move!\n";
       #$newg{$newgeneid}->_dump(\*STDERR);
       push (@final_genes,$newg{$newgeneid});
   }

   # merges are also quite easy.

   foreach my $newgeneid ( keys %merge ) {
       my @newtrans = $newg{$newgeneid}->each_Transcript;
       my @oldtrans;
       my $largest;
       my $size = 0;


       foreach my $oldgeneid ( @{$merge{$newgeneid}} ) {
	   $has_moved_old{$oldgeneid} = 1;
	   if( $oldgeneid ne $largest ) {
	       push (@dead_gene_ids,$oldgeneid);
	   }
       }

       #Into merge code....
       my %old_tg;
       foreach my $oldgeneid ( @{$merge{$newgeneid}} ) {
	   if ($self->{_done_gene_hash}->{$oldgeneid}) {
	       next;
	   }
	   my $tsize = scalar ( $oldg{$oldgeneid}->each_unique_Exon );

	   if( $tsize > $size ) {
	       $largest = $oldgeneid;
	   }
	   foreach my $trans ($oldg{$oldgeneid}->each_Transcript) {
	       $old_tg{$trans->id}=$oldg{$oldgeneid};
	   }
	   push(@oldtrans,$oldg{$oldgeneid}->each_Transcript);
       }
       my $should_increment = $self->map_temp_Transcripts_to_real_Transcripts($newg{$newgeneid},\@newtrans,\@oldtrans,%old_tg);
       
       if ($oldg{$largest}) {
	   $has_done_new{$newgeneid} = 1;
	   # deal with the mapping of ids
	   print $log "GENE MAP: MERGED ",$newgeneid," ",$oldg{$largest}->id,"\n";
	   $self->{_done_gene_hash}->{$newgeneid}=1;
	   $self->{_done_gene_hash}->{$oldg{$largest}->id}=1;
	   $newg{$newgeneid}->id($oldg{$largest}->id);
	   
	   # we increment irregardless of anything else
	   $newg{$newgeneid}->version($oldg{$largest}->version+1);
	   $newg{$newgeneid}->modified($now);
	   push (@final_genes,$newg{$newgeneid});
       }
   }

   # splits **suck** big time
    foreach my $oldgeneid ( keys %split ) {
	my @newgeneid = @{$split{$oldgeneid}};

       # we take old transcripts one at a time, and fit to all possible
       # new transcripts. We take the best. The first case wins the gene id 
       #and the rest get assigned new geneids

	my $assigned = 0;
	
	foreach my $trans ( $oldg{$oldgeneid}->each_Transcript ) {
	    print $log "Checking if ".$trans->id." has been used\n";
	    if (exists $self->{_fitted_trans_hash}->{$trans->id}) {
		print $log "Skipping ".$trans->id."\n";
		next;
	    }
	    my $score = 0;
	    my $current_fit = undef;
	    
	    # needs to up here so we can assign it later on.
	    

	   # flag that we have moved these
	   $has_moved_old{$oldgeneid} =1;
	   
	   
	   foreach my $newgeneid ( @newgeneid ) {
	       $has_done_new{$newgeneid} = 1;
	   }

	   my $lastgeneid;
	   GENE: foreach my $newgeneid ( @newgeneid ) {
	       foreach my $newtrans ( $newg{$newgeneid}->each_Transcript ) { 
		   if( exists $self->{_fitted_trans_hash}->{$newtrans->id} ) {
		       next;
		   }

		   my ($tscore,$perfect) = Bio::EnsEMBL::Pipeline::GeneComp::overlap_Transcript($trans,$newtrans);
		   if( $tscore > $score || $perfect == 1) {
		       $current_fit = $newtrans;
		       $lastgeneid = $newgeneid;
		   }

		   if( $perfect == 1 ) {
		       $lastgeneid=$newgeneid;
		       # just in case
		       $current_fit = $newtrans;
		       last GENE;
		   }
	       }
	   }
	   if( !defined $current_fit ) {
	       #Dead, archive!
	       #my $fe = $trans->start_exon;
	       #my $clone = $vc->dbobj->_crossdb->old_dbobj->get_Clone($fe->clone_id);
	       #my $old_trans = $vc->dbobj->_crossdb->old_dbobj->gene_Obj->get_Transcript($trans->id);
	       #$self->archive->write_seq($old_trans->dna_seq,$old_trans->version,'transcript',$oldgeneid, $oldg{$oldgeneid}->version,$clone->id,$clone->embl_version);
	       #$self->archive->write_seq($old_trans->translate,$old_trans->translation->version,'protein',$oldgeneid,$oldg{$oldgeneid}->version,$clone->id,$clone->embl_version);
	       print $log "TRANSCRIPT MAP: KILLED ",$trans->id,"\n";
	       print $log "TRANSLATION MAP: KILLED ".$trans->translation->id."\n";
	       next;
	   }

	   my $tempid=$current_fit->id;
	   # current_fit is the best fit of this old transcript
	   $current_fit->id($trans->id);
	   if( Bio::EnsEMBL::Pipeline::GeneComp::increment_Transcript($trans,$current_fit) == 1 ) {
	       $current_fit->version($trans->version()+1);
	   } else {
	       $current_fit->version($trans->version());
	   }
	   print $log "TRANSCRIPT MAP: SPLIT $tempid ".$current_fit->id."\n";
	   $self->{_fitted_trans_hash}->{$current_fit->id} = 1;
	   $self->{_fitted_trans_hash}->{$tempid} = 1;
	   #Deal with translation
	   my $tempid=$current_fit->translation->id;
	   my $id=$current_fit->id;
	   $id =~ s/ENST/ENSP/;
	   $current_fit->translation->id($id);
	   my $oldid = $current_fit->translation->start_exon_id();
	   if ($exon_map->{$oldid}) {
	       $current_fit->translation->start_exon_id($exon_map->{$oldid});
	   }
	   else {
	       print $log "START EXON BUG: Could not map start exon ".$oldid." for translation $tempid\n";
	   } 
	   $oldid = $current_fit->translation->end_exon_id();
	   if ($exon_map->{$oldid}){
	       $current_fit->translation->end_exon_id($exon_map->{$oldid});
	   }
	   else {
	       print $log "END EXON BUG: Could not map end exon ".$oldid." for translation $tempid\n";
	   }  
	   print $log "TRANSLATION MAP: ".$tempid." ".$current_fit->translation->id."\n";
	   
	   if($assigned == 0) {
	       # this gene id wins. Hurray!
	       print $log "GENE MAP: SPLIT ",$lastgeneid," ",$oldgeneid,"\n";
	       $self->{_done_gene_hash}->{$lastgeneid}=1;
	       $newg{$lastgeneid}->id($oldgeneid);
	       $newg{$lastgeneid}->version($oldg{$oldgeneid}->version+1);
	       $newg{$lastgeneid}->modified($now);
	       $assigned=1;
	       push (@final_genes,$newg{$lastgeneid});
	   }

       }
   
       # we need to create new transcripts for the remainder of the new transcripts not
       # fitted
       foreach my $newgeneid ( @newgeneid ) {
	   foreach my $newtrans ( $newg{$newgeneid}->each_Transcript ) { 
	       if(exists $self->{_fitted_trans_hash}->{$newtrans->id}) {
		   next;
	       }
	       my $tempid = $newtrans->id; 
	       $newtrans->id($self->archive->get_new_stable_ids('transcript',1));
	       $self->{_fitted_trans_hash}->{$newtrans->id}=1;
	       print $log "TRANSCRIPT MAP: SPLIT ",$tempid," ",$newtrans->id,"\n";
	       
	       $newtrans->version(1);
	       $newtrans->created($now);
	       $newtrans->modified($now);
	       
	       #Deal with the translation mapping
	       my $tempid=$newtrans->translation->id();
	       $newtrans->translation->id($self->archive->get_new_stable_ids('translation',1));
	       my $oldid = $newtrans->translation->start_exon_id();
	       if ($exon_map->{$oldid}) {
	          $newtrans->translation->start_exon_id($exon_map->{$oldid});
	       }
	       else {
		   print $log "START EXON BUG: Could not map start exon ".$oldid." for translation $tempid-".$newtrans->translation->id." on transcript ".$newtrans->id."\n";
	       }
	       $oldid = $newtrans->translation->end_exon_id();
	       if ($exon_map->{$oldid}) {
		   $newtrans->translation->end_exon_id($exon_map->{$oldid});
	       }
	       else {
		   print $log "END EXON BUG: Could not map end exon ".$oldid." for translation $tempid-".$newtrans->translation->id." on transcript ".$newtrans->id."\n";
	       } 
		   print $log "TRANSLATION MAP: SPLIT ".$newtrans->translation->id()." ".$tempid."\n";
	   }
	   if(($newg{$newgeneid}->id ne $oldgeneid) && ($newgeneid !~ /ENSG/)) {
	       if ($self->{_done_gene_hash}->{$newgeneid}) {
		   next;
	       }
	       my $newgene = $newg{$newgeneid};
	       my $tempid = $newgene->id();
	       # it is an unassigned gene...
	       $newgene->id($self->archive->get_new_stable_ids('gene',1));
	       print $log "GENE MAP: SPLIT NEW ",$tempid," ",$newgene->id,"\n";
	       $self->{_done_gene_hash}->{$tempid}=1;
	       $newgene->created($now);
	       $newgene->modified($now);
	       $newgene->version(1);
	       $self->{_done_gene_hash}->{$newgeneid}=1;
	       push (@final_genes,$newgene);
	   }
       }


   }


   # Now - handle all the cases which have not been handled already, and assign
   # new everything to them!

   foreach my $newgene_id ( keys %newg ) {
       if( $has_done_new{$newgene_id} ) {
	   next;
       }

       my $newgene = $newg{$newgene_id};

       $newgene->id($self->archive->get_new_stable_ids('gene',1));
       $newgene->created($now);
       $newgene->modified($now);
       $newgene->version(1);

       foreach my $t ( $newgene->each_Transcript ) {
	   if ($self->{_fitted_trans_hash}->{$t->id}) {
	       next;
	   }
	   my $tempid=$t->id;
	   $t->id($self->archive->get_new_stable_ids('transcript',1));
	   $t->created($now);
	   $t->modified($now);
	   $t->version(1);
	   $self->{_fitted_trans_hash}->{$t->id}=1;
	   print $log "TRANSCRIPT MAP: LEFT NEW ".$tempid." ".$t->id."\n";
	   my $tempid=$t->translation->id;
	   $t->translation->id($self->archive->get_new_stable_ids('translation',1));
	   my $oldid = $t->translation->start_exon_id();
	   if ($exon_map->{$oldid}) {
	       $t->translation->start_exon_id($exon_map->{$oldid});
	   }
	   else {
	       print $log "START EXON BUG: Could not map start exon ".$oldid." for translation $tempid\n";
	   } 
	   $oldid = $t->translation->end_exon_id();
	   if ($exon_map->{$oldid}) {
	       $t->translation->end_exon_id($exon_map->{$oldid});
	   }
	   else {
	       print $log "END EXON BUG: Could not map end exon ".$oldid." for translation $tempid\n";
	   }  
	   print $log "TRANSLATION MAP: NEW ".$tempid." ".$t->translation->id."\n";
       }
       print $log "GENE MAP: NEW ".$newgene_id." ".$newgene->id."\n";
       $self->{_done_gene_hash}->{$newgene_id}=1;
       push (@final_genes,$newgene);
   }
    
    foreach my $gene ( @oldgenes ) {
	if( $has_moved_old{$gene->id} ) {
	    next;
	}
	push(@dead_gene_ids,$gene->id);
    }
    #Making array of genes unique, is this normal? hmm.,...
    my %seen = ();
    my @unique = grep { ! $seen{$_} ++ } @dead_gene_ids;

    #Archiving dead genes
    #foreach my $gene (@unique) {
	#$self->archive->write_dead_geneid($gene);
    #}
    my $final_size=scalar(@final_genes);
    return ($temp_size,\@final_genes);
}


=head2 map_temp_Transcripts_to_real_Transcripts

 Title   : map_temp_Transcripts_to_real_Transcripts
 Usage   : $self->map_temp_Transcripts_to_real_Transcripts
 Function: mapping of temp transcripts to old ones
 Example : 
 Returns : nothing
 Args    : Old gene object, array of new transcripts, 
           array of old transcripts

=cut

sub map_temp_Transcripts_to_real_Transcripts{
   my ($self,$oldgene,$new,$old,%old_tg) = @_;

   my $vc = $self->vc();
   my $log = $self->log();
   my @newt = @$new;
   my @oldt = @$old;
   my $should_change = 0;
   my @dead;
   my $exon_map = $self->{'_exon_map_hash'};
   my %newt;
   
   foreach my $t ( @newt ) {
       $newt{$t->id} = $t;
   }
   

   # sort old by number of exons - largest first
   @oldt = sort { $b->number <=> $a->number } @oldt;
   my $now = time();

   # best fit to new transcripts...
   foreach my $oldt ( @oldt ) {
       if ($self->{_fitted_trans_hash}->{$oldt->id}) {
	   next;
       }
       my $score = 0;
       my $fitted_trans = undef;

       foreach my $newt ( @newt ) {
	   if( exists $self->{_fitted_trans_hash}->{$newt->id} ) {
	       next;
	   }
	   my ($tscore,$perfect) = Bio::EnsEMBL::Pipeline::GeneComp::overlap_Transcript($oldt,$newt);
	   if( $perfect ) {
	       $fitted_trans = $newt;
	       last;
	   }
	   if( $tscore > $score ) {
	       $fitted_trans = $newt;
	   }
       }

       if ( defined $fitted_trans ) {
	   print $log "TRANSCRIPT MAP: ",$fitted_trans->id," ",$oldt->id,"\n";

	   $fitted_trans->id($oldt->id);
	   $self->{_fitted_trans_hash}->{$fitted_trans->id}=1;
	   #print $log "Setting ".$oldt->id." to 1\n";
	   $self->{_fitted_trans_hash}->{$oldt->id}=1;
	   if( Bio::EnsEMBL::Pipeline::GeneComp::increment_Transcript($oldt,$fitted_trans) == 1 ) {
	       $fitted_trans->version($oldt->version()+1);
	       $fitted_trans->modified($now);
	       $should_change = 1;
	   } else {
	       $fitted_trans->version($oldt->version());
	   }

	   #Deal with translation
	   my $tempid=$fitted_trans->translation->id;
	   my $id=$fitted_trans->id;
	   $id =~ s/ENST/ENSP/;
	   $fitted_trans->translation->id($id);
	   my $oldid = $fitted_trans->translation->start_exon_id();
	   if ($exon_map->{$oldid}) {
	       $fitted_trans->translation->start_exon_id($exon_map->{$oldid});
	   }
	   else {
	       print $log "START EXON BUG: Could not map start exon ".$oldid." for translation $tempid\n";
	   } 
	   $oldid = $fitted_trans->translation->end_exon_id();
	   if ($exon_map->{$oldid}){
	       $fitted_trans->translation->end_exon_id($exon_map->{$oldid});
	   }
	   else {
	        print $log "END EXON BUG: Could not map end exon ".$oldid." for translation $tempid\n";
	   }  
	   print $log "TRANSLATION MAP: ".$tempid." ".$fitted_trans->translation->id."\n";
       } else {
	   # it is dead ;)
	   $should_change = 1;
	   #my $gene=$old_tg{$oldt->id};
	   #my $fe = $oldt->start_exon;
	   #my $clone = $vc->dbobj->_crossdb->old_dbobj->get_Clone($fe->clone_id);
	   #my $old_trans = $vc->dbobj->_crossdb->old_dbobj->gene_Obj->get_Transcript($oldt->id);
	   #$self->archive->write_seq($old_trans->dna_seq,$old_trans->version,'transcript',$gene->id,$gene->version,$clone->id,$clone->version);
	   #$self->archive->write_seq($old_trans->translate,$old_trans->translation->version,'protein',$gene->id,$gene->version,$clone->id,$clone->version);
	   print $log "TRANSCRIPT MAP: KILLED ",$oldt->id,"\n";
	   print $log "TRANSLATION MAP: KILLED ".$oldt->translation->id."\n";
       }
   }

   foreach my $newt ( @newt ) {
       if( exists $self->{_fitted_trans_hash}->{$newt->id} ) {
	   next;
       }
       my $tempid = $newt->id;
       $should_change = 1;
       $newt->id($self->archive->get_new_stable_ids('transcript',1));
       $newt->version(1);
       $newt->created($now);
       $newt->modified($now);
       print $log "TRANSCRIPT MAP: NEW ",$tempid," ",$newt->id,"\n";
       $self->{_fitted_trans_hash}->{$newt->id}=1;
       $tempid = $newt->translation->id;
       $newt->translation->id($self->archive->get_new_stable_ids('translation',1));
       my $oldid = $newt->translation->start_exon_id();
       if ($exon_map->{$oldid}) {
	   $newt->translation->start_exon_id($exon_map->{$oldid});
       }
       else {
	    print $log "START EXON BUG: Could not map start exon ".$oldid." for translation $tempid\n";
       }
       $oldid = $newt->translation->end_exon_id();
       if ($exon_map->{$oldid}) {
	   $newt->translation->end_exon_id($exon_map->{$oldid});
       }
       else {
	    print $log "END EXON BUG: Could not map end exon ".$oldid." for translation $tempid\n";
       }
       print $log "TRANSLATION MAP: NEW ".$tempid." ".$newt->translation->id."\n";
   }
   return $should_change;
}


=head2 overlap_Transcript

 Title   : overlap_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub overlap_Transcript{
   my ($old,$new) = @_;
   my ($score,$perfect);

   if( !defined $old || !defined $new || ! ref $old || !$old->isa('Bio::EnsEMBL::Transcript') || !$new->isa('Bio::EnsEMBL::Transcript')) {
       croak ('Did not give me both old and new transcripts in overlap Transcript');
   }

   
   $perfect = 1;

   my ($i,$j);
   my @newe = $new->each_Exon();
   my @olde = $old->each_Exon();
   my $jj;

   MAIN_LOOP:

   for($i=0,$j=0;$i<= $#olde && $j <= $#newe ;) {
       if( $olde[$i]->id eq $newe[$j]->id ) {
	   $score++;
	   $i++; $j++; next;
       }
       # can't be perfect if we get into here...
       $perfect = 0;
       # see whether this exon is anywhere here...
       for($jj=$j+1;$jj <= $#newe;$jj++) {
	   if( $olde[$i]->id eq $newe[$jj]->id ) {
	       $j=$jj;
	       next MAIN_LOOP;
	   }
       }
       # move i along
       $i++;

   }

   return ($score,$perfect);
	   
}


=head2 increment_Transcript

 Title   : increment_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub increment_Transcript{
   my ($old,$new) = @_;

   if( !defined $old || !defined $new || ! ref $old || !$old->isa('Bio::EnsEMBL::Transcript') || !$new->isa('Bio::EnsEMBL::Transcript')) {
       croak ('Did not give me both old and new transcripts in increment Transcript');
   }

   my ($i,$j);
   my @newe = $new->each_Exon();
   my @olde = $old->each_Exon();

   if( $#newe != $#olde ) {
       return 1;
   }

   my $jj;

 MAIN_LOOP:

   for($i=0,$j=0;$i<= $#olde && $j <= $#newe ;) {
       #if( $olde[$i]->id eq $newe[$j]->id && $olde[$i]->has_changed_version == 0 && $newe[$j]->has_changed_version == 0 ) {
       if( $olde[$i]->id eq $newe[$j]->id && ($olde[$i]->version != $newe[$j]->version) ) {
	   $i++; $j++; next;
       }
       last;
   }

   if( $i <= $#olde ) {
       return 1;
   }

   return 0;
}


=head1 Get/Set functions


=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)


=cut

sub vc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'vc'} = $value;
    }
    return $obj->{'vc'};

}

=head2 archive

 Title   : archive
 Usage   : $obj->archive($newval)
 Function: 
 Returns : value of archive
 Args    : newvalue (optional)


=cut

sub archive{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'archive'} = $value;
    }
    return $obj->{'archive'};

}

#At the moment, no final db writing...
#=head2 finaldb
#
# Title   : finaldb
# Usage   : $obj->finaldb($newval)
# Function: 
# Returns : value of finaldb
# Args    : newvalue (optional)
#
#
#=cut
#
#sub finaldb{
#   my $obj = shift;
#   if( @_ ) {
#      my $value = shift;
#      $obj->{'finaldb'} = $value;
#    }
#    return $obj->{'finaldb'};
#
#}

=head2 log

 Title   : log
 Usage   : $obj->log($newval)
 Function: 
 Returns : value of log
 Args    : newvalue (optional)


=cut

sub log{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'log'} = $value;
    }
    return $obj->{'log'};

}
=head2 mapfile

 Title   : mapfile
 Usage   : $obj->mapfile($newval)
 Function: 
 Returns : value of mapfile
 Args    : newvalue (optional)


=cut

sub mapfile{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'mapfile'} = $value;
    }
    return $obj->{'mapfile'};

}
1;
