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

    $vc is a Bio::EnsEMBL::Virtual::Contig (it has to be attached to a
    crossmatch database, in order to support the get_old_Exons and
    get_old_Genes calls).

 
    $arcdb is the archive database to which all dead genes,
    transcripts and translations will be written.

    $finaldb is a new database object to which all new, mapped genes will be
    written (with all their transcripts, exons, etc.)
    
    $log is the logfile filehandle to which all mappings are logged

    my $gc =  Bio::EnsEMBL::Pipeline::GeneComp->new('-vc' => $vc,
						    '-archive' => $arcdb,
						    '-finaldb' => $finaldb,
						    '-log' => \*LOG);

    $gc->map();

=head1 DESCRIPTION

This object deals with mapping exons, transcript and genes from old 
versions through to new versions. To do the mapping we need to get 
the old exons out in the new coordinates (remapping). This currently 
is hidden behind the method call get_old_Exons() on a contig object. 
This call returns old exon objects in the new coordinates, by getting
crossmatch feature pairs from a crossmatch database.

All in all this object holds on to 4 databases. It holds on to:

1-Temp (new) database ($vc->dbobj)
2-Old database ($vc->dbobj->_crossdb->old_dbobj->)
3-Crossmatch database ($vc->dbobj->_crossdb)
4-Archive database ($self->archive)

And it also hold onto two files, a log file and a map file. The map file is a tab delimited temp-mapped text file.

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
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Root;
use Carp;

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

 Title   : new
 Usage   : $genecomp = Bio::EnsEMBL::Pipeline::GeneComp->new(-vc => $vc,
                      -archive => $archive,-finaldb => $db,-log => $log);
 Function: Builds a new genecomp object, based around a virtual contig, 
 Example :
 Returns : a new genecomp object
 Args    : 


=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($vc,$archive,$log,$map,$maphref) = 
      $self->_rearrange([qw(VC
			    ARCHIVE
			    LOG
			    MAP
			    HASHREF
			    )],@args);

  #REmoved finaldb, now writing....

  if( !defined $vc || !ref $vc || 
      !$vc->isa('Bio::EnsEMBL::Virtual::Contig') ) {
      $self->throw("must have a virtual contig, got a [$vc]");
  }

  if( !defined $archive || !ref $archive || 
      !$archive->isa('Bio::EnsEMBL::DBArchive::Obj') ) {
      $self->throw("must have an archive database, got a [$archive]");
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
  $self->maphref($maphref);
  $self->{'_fitted_trans_hash'} = ();

  # this usually dies if we can't write
  print $log "Built GeneComp object\n";

  return $self;
}


=head2 maphref

 Title   : maphref
 Usage   : $obj->maphref($newval)
 Function: Getset for maphref value
 Returns : value of maphref
 Args    : newvalue (optional)


=cut

sub maphref{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'maphref'} = $value;
    }
    return $obj->{'maphref'};

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
    my $vc = $self->vc();
    
    #Map exons
    print $log "Mapping exons:\n";
    $self->map_temp_Exons_to_real_Exons();
    #Map genes
    print $log "Mapping genes:\n";
    my ($before_g,$finalgenes) = $self->map_temp_Genes_to_real_Genes();
    my $after_g = scalar (@$finalgenes);
    print $log "Temp genes: $before_g Final genes: $after_g\n";
    if ($before_g != $after_g) {
	print $log "GENE MAPPING BUG: Size of temp genes and final genes not equal, bad bug! Temp genes: $before_g, final: $after_g";
    }
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
   my $arc = $self->archive;
   my $vc = $self->vc();
   my $log = $self->log;
   my $map = $self->mapfile;
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

   my %moved;    # shows which old exons we have moved (or not).
   my %ismapped; # tells us whether temporary exons have been mapped or not
   my %temp_old; # Hash of temp to old mapping to be used later in the gene code
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
   print $log "Getting all old exons...";
   my @oldexons=$vc->get_old_Exons($log,$self->maphref);

   my $size=scalar(@oldexons);
   print $log " got $size old exons\n";
   # get out one date for this mapping....

   my $time = time();

   @tempexons  = sort { $a->start <=> $b->start } @tempexons;

   @oldexons   = sort { $a->end   <=> $b->end } @oldexons;

   #if(my $dd == 1){
   #    my @all = @tempexons;
   #    push(@all,@oldexons);
   #    @all   = sort { $a->end   <=> $b->end } @all;
       #foreach my $exon (@all){
	   #print"  ". $exon->id." ".$exon->start." ".$exon->end."\n";
       #}
       #exit 0;
   #}
   my %e_v;
   foreach my $oldexon ( @oldexons ) {
       $e_v{$oldexon->id}=$oldexon->version;
   }
   

   my @tempexons2;
   my %bestfit;
   # go over each exon and map old->new...

   #First sort out identical exons...
   print $log "Checking for identical exons...\n";
   IDENTICAL:foreach my $tempexon ( @tempexons ) {
       #print $log "Trying to find identical exon to ".$tempexon->id." (".$tempexon->start."-".$tempexon->end.")\n";
       foreach my $oldexon ( @oldexons ) {
	   if ($tempexon->start > $oldexon->end) {
	       #Skip this old exon, it ends before temp starts
	       next;
	   }
	   if( $oldexon->start > $tempexon->end ) {
	       #Exiting comparison loop because this old exon 
	       #starts after end of temp
	       if (! $bestfit{$tempexon->id}) {
		   push (@tempexons2,$tempexon);
		   $bestfit{$tempexon->id}=1;
	       }
	       next IDENTICAL;
	   }
	   #print $log "Checking it against ".$oldexon->id." (".$oldexon->start."-".$oldexon->end.")\n";
	   # if start/end/strand is identical, it is definitely up for moving.
	   
	   if( $tempexon->start == $oldexon->start &&
	       $tempexon->end   == $oldexon->end   &&
	       $tempexon->strand == $oldexon->strand ) {

	       # ok - if we have already moved this - ERROR
	       if( $moved{$oldexon->id} == 1 ) {
		   print $log "Duplicate! ".$tempexon->id." maps to ".$oldexon->id." with identical start/end/strand. ".$tempexon->start.":".$tempexon->end.":".$tempexon->strand."\n";
		   if (!$bestfit{$tempexon->id} ) {
		       push (@tempexons2,$tempexon);
		       $bestfit{$tempexon->id}=1;
		   }
		   next IDENTICAL;
	       }

	       # set ismapped to 1 for this tempexon

	       $ismapped{$tempexon->id} = 1;
	       $temp_old{$tempexon->id}=$oldexon->id;
	       
	       print $log "EXON MAP IDENTICAL: ",$tempexon->id," ",$oldexon->id," "; 
	       print $map $tempexon->id,"\t",$oldexon->id."\t";
               
               # will add version
	       # we move the id. Do we move the version?
	       $tempexon->id($oldexon->id);
	       
	       if( $oldexon->seq eq $tempexon->seq) {
		   # version and id
		   $tempexon->version($oldexon->version);
		   # don't update modified
	       } else {
		   my $v=$oldexon->version()+1;
		   $tempexon->version($v);
		   $tempexon->modified($time);
	       }
	       push (@mapped,$tempexon);
	       print $log $tempexon->version."\n";
	       print $map $tempexon->version."\n";
	       $moved{$oldexon->id} = 1;
	       next IDENTICAL;
	   }
       }
       if (!$bestfit{$tempexon->id} ) {
	   push (@tempexons2,$tempexon);
	   $bestfit{$tempexon->id}=1;
       }
   }

   print $log "Looking for BEST fits for remaining exons...\n";
   my $size = scalar (@tempexons2);
   print $log "Going to do bestfit/new on $size remaining exons...\n";
   $size = scalar (@oldexons);

   # we can't map some exons directly. Second scan over all exons to see whether
   # there is significant overlap
   #Hash to hold all alternative bestfits, 
   #to archive those not used (in split cases)
   my %bf;
   
 BESTFIT:foreach my $tempexon2 (@tempexons2) {
     my $biggestoverlap = undef;
     my $overlapsize = 0;
     #print $log "BESTFIT Trying to find overlapping exon to ".$tempexon2->id." (".$tempexon2->start."-".$tempexon2->end.")\n";
     if ($ismapped{$tempexon2->id} == 1) {
	 next;
     }
     foreach my $oldexon ( @oldexons ) {
	 
	 if ($tempexon2->start > $oldexon->end) {
	     #Skip this old exon, it ends before temp starts
	     next;
	 }
	 if( $oldexon->start > $tempexon2->end ) {
	     #Go out of loop, reached old exon after temp
	     last;
	 }
	 #print $log "Checking it against ".$oldexon->id." (".$oldexon->start."-".$oldexon->end.")\n";
	 if( (my $intexon = $oldexon->intersection($tempexon2)) && ($moved{$oldexon->id} == 0) ) {
	     $bf{$oldexon->id}=$tempexon2->id;
	     my $tovsize = ($intexon->end - $intexon->start +1); 
	     if( !defined $biggestoverlap ) {
		 $biggestoverlap = $oldexon;
		 $overlapsize = $tovsize; 
	     } else {
		 if( $tovsize > $overlapsize ) {
		     $biggestoverlap = $oldexon;
		     $overlapsize = $tovsize;
		 }
	     }
	 }
     }
     # if I have got a biggest overlap - map across...
     if( defined $biggestoverlap ) {
	 # set ismapped to 1 for this tempexon2
	 
	 $ismapped{$tempexon2->id} = 1;
	 delete $bf{$biggestoverlap->id};
	 # we move the id. Do we move the version?
	 print $log "EXON MAP BEST FIT: ",$tempexon2->id," ",$biggestoverlap->id." ";
	 print $map $tempexon2->id,"\t",$biggestoverlap->id."\t";
	 
	 $temp_old{$tempexon2->id}=$biggestoverlap->id;
	 $tempexon2->id($biggestoverlap->id);
	 $tempexon2->version($biggestoverlap->version()+1);
	 print $log $tempexon2->version."\n";
	 print $map $tempexon2->version."\n";
	 
	 $tempexon2->modified($time);
	 push(@mapped,$tempexon2);
	 $moved{$biggestoverlap->id} = 1;
	 my $size=scalar (@oldexons);
	 while( my $tempoldexon = shift @oldexons ) {
	     if( $tempoldexon->end >= $tempexon2->start ) {
		 unshift(@oldexons,$tempoldexon);
		 last;
	     }
	 }
	 
	 next BESTFIT;
     } else {
	 # ok - new Exon
	 my $tempid = $tempexon2->id();
	 $tempexon2->id($self->archive->get_new_stable_ids('exon',1));
	 $tempexon2->created($time);
	 $tempexon2->modified($time);
	 $tempexon2->version(1);
	 $temp_old{$tempid}=$tempexon2->id;
	 print $log "EXON MAP NEW: ",$tempid," ",$tempexon2->id," 1\n"; 
	 print $map $tempid,"\t",$tempexon2->id,"\t1\n"; 
	 push(@new,$tempexon2);
	 next BESTFIT;
     }
    $self->throw("Error - should never reach here!");
}
   foreach my $alternative (keys (%bf)) {
       my $tempid = $bf{$alternative};
       my $v = $e_v{$alternative};
       $arc->write_deleted_id('exon',$alternative,$v,$temp_old{$tempid});
       print $log "EXON MAP KILLED $alternative (ALTERNATIVE BESTFIT ".$temp_old{$tempid}.")\n";
   }
   $self->{'_exon_map_hash'} = \%temp_old;
   $self->{'_mapped_exons'} = \@mapped;
   $self->{'_new_exons'} = \@new;

   my $ns = scalar @new;
   my $nm = scalar @mapped;
   
   my $f_size= $ns+$nm;
   print STDERR "Temp size: $t_size, final size: $f_size (new:$ns mapped:$nm)\n";
   if ($t_size != $f_size) {
       print $log "EXON MAPPING BUG: Size of temp exons and final exons not equal, bad bug! Temp exons: $t_size, final: $f_size";
   }
   return 1;
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

    my $map = $self->mapfile;
    my $vc= $self->vc();
    my $log = $self->log();
    my $exon_map = $self->{'_exon_map_hash'};
    
    if( !defined $vc) {
	die("Not passing the Virtual Contig to do gene mapping");
    }
    print $log "Getting all new genes... ";
    my @tempgenes = $vc->get_all_Genes();
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
    my @oldgenes  = $vc->get_old_Genes($self->maphref);
    my $size=scalar(@oldgenes);
    print $log "got $size old genes\n";
    if (!$temp_size && $size) {
	print $log "No genes to map in the new database!\n";
    }
	
    my @final_genes; #final set of genes with mapped and new ids
    my %olde2t; # hash of exon id to array of transcript id
    my %olde2g; # hash of exon id to gene id
    my %oldt;   # hash of transcript on transcript id (effectively t->e mapping)
    my %oldg;   # hash of gene on gene id (effectively g->e mapping)
    my %newg;   # hash of new genes on geneid
    
    my %has_done_new;  # 1 or 0 depending on whether this gene has be moved or not
    my %has_moved_old; # 1 or 0 depending on whether this gene has mapped forward or not
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
    my %killed; #Holds a record of genes archived in merge code
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
	print $log "GENE MAP: MOVED $newgeneid $oldgeneid ";
	print $map $newgeneid,"\t",$oldgeneid."\t";
	$self->{'_done_gene_hash'}->{$newgeneid}=1;
	$self->{'_done_gene_hash'}->{$oldgeneid}=1;
	$newg{$newgeneid}->id($oldgeneid);
	
	if( $should_increment ) {
	    $newg{$newgeneid}->version($oldg{$oldgeneid}->version+1);
	    $newg{$newgeneid}->modified($now);
	}
	print $log $newg{$newgeneid}->version."\n";
	print $map $newg{$newgeneid}->version."\n";
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

	#Into merge code....
	my %old_tg;
	foreach my $oldgeneid ( @{$merge{$newgeneid}} ) {
	    if ($self->{'_done_gene_hash'}->{$oldgeneid}) {
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
	
	if ($oldg{$largest}) {
	    $has_done_new{$newgeneid} = 1;
	    # deal with the mapping of ids
	    print $log "GENE MAP: MERGED $newgeneid ",$oldg{$largest}->id." ";
	    print $map $newgeneid,"\t",$oldg{$largest}->id."\t";
	    $self->{'_done_gene_hash'}->{$newgeneid}=1;
	    $self->{'_done_gene_hash'}->{$oldg{$largest}->id}=1;
	    $newg{$newgeneid}->id($oldg{$largest}->id);
	    
	    # we increment irregardless of anything else
	    $newg{$newgeneid}->version($oldg{$largest}->version+1);
	    print $log $newg{$newgeneid}->version()."\n";
	    print $map $newg{$newgeneid}->version()."\n";
	    $newg{$newgeneid}->modified($now);
	    push (@final_genes,$newg{$newgeneid});
	}

	foreach my $oldgeneid ( @{$merge{$newgeneid}} ) {
	    $has_moved_old{$oldgeneid} = 1;
	    if( $oldgeneid ne $largest ) {
		#Double check this is right...
		#Killing all genes in merge that are not largest
		print $log  "GENE MAP KILLED: $oldgeneid merged into $largest\n";
		my $oldv = $oldg{$oldgeneid}->version;
		$self->archive->write_deleted_id('gene',$oldgeneid,$oldv,$largest);
		$killed{$oldgeneid}=1;
	    }
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
	    if (exists $self->{'_fitted_trans_hash'}->{$trans->id}) {
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
		if( exists $self->{'_fitted_trans_hash'}->{$newtrans->id} ) {
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
		my $fe = $trans->start_exon;
		eval {
		    my $clone = $vc->dbobj->_crossdb->old_dbobj->get_Clone($fe->clone_id);
		    my $old_trans = $vc->dbobj->_crossdb->old_dbobj->gene_Obj->get_Transcript($trans->id);
		    $self->archive->write_seq($old_trans->dna_seq,$old_trans->version,'transcript',$oldgeneid, $oldg{$oldgeneid}->version,$clone->id,$clone->embl_version);
		    $self->archive->write_seq($old_trans->translate,$old_trans->translation->version,'protein',$oldgeneid,$oldg{$oldgeneid}->version,$clone->id,$clone->embl_version);
		};
		$self->warn("Could not archive ".$trans->id." because of $@\n");
		print $log "TRANSCRIPT MAP: KILLED ",$trans->id,"\n";
		print $log "TRANSLATION MAP: KILLED ".$trans->translation->id."\n";
		$self->archive->write_deleted_id('transcript',$trans->id,$trans->version);
		$self->archive->write_deleted_id('translation',$trans->id,$trans->version);
		next;
	    }
	    
	    my $tempid=$current_fit->id;
	    # current_fit is the best fit of this old transcript
	    $current_fit->id($trans->id);
	    my $should_increment = $self->increment_Transcript($trans,$current_fit);
	    my $v;
	    if( $should_increment == 1 ) {
		$v=$trans->version()+1;
	    } else {
		$v=$trans->version();
	    }
	    $current_fit->version($v);
	    print $log "TRANSCRIPT MAP: SPLIT $tempid ".$current_fit->id." ".$current_fit->version()."\n";
	    print $map "$tempid\t".$current_fit->id."\t".$current_fit->version()."\n";
	    $self->{'_fitted_trans_hash'}->{$current_fit->id} = 1;
	    $self->{'_fitted_trans_hash'}->{$tempid} = 1;
	    #Deal with translation
	    my $tempid=$current_fit->translation->id;
	    my $id=$current_fit->id;
	    $id =~ s/ENST/ENSP/;
	    $current_fit->translation->id($id);
	    $current_fit->translation->version($v);
	    print $log "TRANSLATION MAP: ".$tempid." ".$current_fit->translation->id." $v\n";
	    print $map "$tempid\t".$current_fit->translation->id."\t$v\n";
	    
	    if($assigned == 0) {
		# this gene id wins. Hurray!
		print $log "GENE MAP: SPLIT $lastgeneid $oldgeneid ";
		print $map $lastgeneid,"\t",$oldgeneid."\t";
		$self->{'_done_gene_hash'}->{$lastgeneid}=1;
		$newg{$lastgeneid}->id($oldgeneid);
		$newg{$lastgeneid}->version($oldg{$oldgeneid}->version+1);
		$newg{$lastgeneid}->modified($now);
		$assigned=1;
		print $log $newg{$lastgeneid}->version."\n";
		print $map $newg{$lastgeneid}->version."\n";
		push (@final_genes,$newg{$lastgeneid});
	    }
	    
	}
	
	# we need to create new transcripts for the remainder of the 
	# new transcripts not fitted
	foreach my $newgeneid ( @newgeneid ) {
	    foreach my $newtrans ( $newg{$newgeneid}->each_Transcript ) { 
		if(exists $self->{'_fitted_trans_hash'}->{$newtrans->id}) {
		    next;
		}
		my $tempid = $newtrans->id; 
		$newtrans->id($self->archive->get_new_stable_ids('transcript',1));
		$self->{'_fitted_trans_hash'}->{$newtrans->id}=1;
		print $log "TRANSCRIPT MAP: SPLIT NEW ",$tempid," ",$newtrans->id," 1\n";
		print $map $tempid,"\t",$newtrans->id,"\t1\n";
		$newtrans->version(1);
		$newtrans->created($now);
		$newtrans->modified($now);
		
		#Deal with the translation mapping
		my $tempid=$newtrans->translation->id();
		my $id=$newtrans->id;
		$id =~ s/ENST/ENSP/;
		$newtrans->translation->id($id);
		$newtrans->translation->version(1);
		
		print $log "TRANSLATION MAP: SPLIT $tempid ".$newtrans->translation->id()." 1\n";
		print $map "$tempid\t".$newtrans->translation->id()."\t1\n";
	    }
	    if(($newg{$newgeneid}->id ne $oldgeneid) && ($newgeneid !~ /ENSG/)) {
		if ($self->{'_done_gene_hash'}->{$newgeneid}) {
		    next;
		}
		my $newgene = $newg{$newgeneid};
		my $tempid = $newgene->id();
		# it is an unassigned gene...
		$newgene->id($self->archive->get_new_stable_ids('gene',1));
		print $log "GENE MAP: SPLIT NEW ",$tempid," ",$newgene->id," 1\n";
		print $map $tempid,"\t",$newgene->id,"\t1\n";
		$self->{'_done_gene_hash'}->{$tempid}=1;
		$newgene->created($now);
		$newgene->modified($now);
		$newgene->version(1);
		$self->{'_done_gene_hash'}->{$newgeneid}=1;
		push (@final_genes,$newgene);
	    }
	}
	
	
    }
    
    
    # Now - handle all the cases which have not been handled already, 
    # and assign new everything to them!
    
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
	    if ($self->{'_fitted_trans_hash'}->{$t->id}) {
		next;
	    }
	    my $tempid=$t->id;
	    $t->id($self->archive->get_new_stable_ids('transcript',1));
	    $t->created($now);
	    $t->modified($now);
	    $t->version(1);
	    $self->{'_fitted_trans_hash'}->{$t->id}=1;
	    print $log "TRANSCRIPT MAP: NEW ".$tempid." ".$t->id." 1\n";
	    print $map "$tempid\t".$t->id."\t1\n";
	    my $tempid=$t->translation->id;
	    my $id=$t->id;
	    $id =~ s/ENST/ENSP/;
	    $t->translation->id($id);
	    print $log "TRANSLATION MAP: NEW ".$tempid." ".$t->translation->id." 1\n";
	    print $map $tempid."\t".$t->translation->id."\t1\n"; 
	}
	print $log "GENE MAP: NEW ".$newgene_id." ".$newgene->id." 1\n";
	print $map $newgene_id."\t".$newgene->id."\t1\n";
	$self->{'_done_gene_hash'}->{$newgene_id}=1;
	push (@final_genes,$newgene);
    }

    my @dead_gene_ids;
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
    foreach my $gene (@unique) {
	print $log "GENE MAP KILLED: $gene\n";
	if (!$killed{$gene}) {
	    my $oldv = $oldg{$gene}->version;
	    $self->archive->write_deleted_id('gene',$gene,$oldv);
	}
    }
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

   my $map = $self->mapfile;
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
       if ($self->{'_fitted_trans_hash'}->{$oldt->id}) {
	   next;
       }
       my $score = 0;
       my $fitted_trans = undef;

       foreach my $newt ( @newt ) {
	   if( exists $self->{'_fitted_trans_hash'}->{$newt->id} ) {
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
	   print $log "TRANSCRIPT MAP: ",$fitted_trans->id," ",$oldt->id." ";
	   print $map $fitted_trans->id,"\t",$oldt->id."\t";
	   
	   $fitted_trans->id($oldt->id);
	   $self->{'_fitted_trans_hash'}->{$fitted_trans->id}=1;
	   #print $log "Setting ".$oldt->id." to 1\n";
	   $self->{'_fitted_trans_hash'}->{$oldt->id}=1;
	   my $should_increment = $self->increment_Transcript($oldt,$fitted_trans);
	   my $v;
	   if( $should_increment == 1 ) {
	       $v=$oldt->version()+1;
	       $fitted_trans->modified($now);
	       $should_change = 1;
	   } else {
	       $v=$oldt->version();
	   }
	   $fitted_trans->version($v);
	   print $log "$v\n";
	   print $map "$v\n";
	   #Deal with translation
	   my $tempid=$fitted_trans->translation->id;
	   my $id=$fitted_trans->id;
	   $id =~ s/ENST/ENSP/;
	   $fitted_trans->translation->id($id);
	   $fitted_trans->translation->version($v);
	   print $log "TRANSLATION MAP: ".$tempid." ".$fitted_trans->translation->id." $v\n";
	   print $map $tempid."\t".$fitted_trans->translation->id."\t$v\n";
       } else {
	   # it is dead ;)
	   $should_change = 1;
	   my $gene=$old_tg{$oldt->id};
	   my $fe = $oldt->start_exon;

	   eval {
	       my $clone = $vc->dbobj->_crossdb->old_dbobj->get_Clone($fe->clone_id);
	       my $old_trans = $vc->dbobj->_crossdb->old_dbobj->gene_Obj->get_Transcript($oldt->id);
	       $self->archive->write_seq($old_trans->dna_seq,$old_trans->version,'transcript',$gene->id,$gene->version,$clone->id,$clone->version);
	       $self->archive->write_seq($old_trans->translate,$old_trans->translation->version,'protein',$gene->id,$gene->version,$clone->id,$clone->version);
               $self->archive->write_deleted_id('transcript',$old_trans->id,$old_trans->version);
               $self->archive->write_deleted_id('translation',$old_trans->translation->id,$old_trans->translation->version);
	   };
	   if ($@) {
	       $self->warn("Could not archive ".$oldt->id." because of $@");
	   }
	   
	   print $log "TRANSCRIPT MAP: KILLED ",$oldt->id,"\n";
	   print $log "TRANSLATION MAP: KILLED ".$oldt->translation->id."\n";
       }
   }

   foreach my $newt ( @newt ) {
       if( exists $self->{'_fitted_trans_hash'}->{$newt->id} ) {
	   next;
       }
       my $tempid = $newt->id;
       $should_change = 1;
       $newt->id($self->archive->get_new_stable_ids('transcript',1));
       $newt->version(1);
       $newt->created($now);
       $newt->modified($now);
       print $log "TRANSCRIPT MAP: NEW ",$tempid," ",$newt->id," 1\n";
       print $map $tempid,"\t",$newt->id,"\t1\n";
       $self->{'_fitted_trans_hash'}->{$newt->id}=1;
       $tempid = $newt->translation->id;
       my $id=$newt->id;
       $id =~ s/ENST/ENSP/;
       $newt->translation->id($id);
       print $log "TRANSLATION MAP: NEW ".$tempid." ".$newt->translation->id." 1\n";
       print $map $tempid."\t".$newt->translation->id."\t1\n";
       
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

   if( !defined $old || !defined $new || ! ref $old || 
       !$old->isa('Bio::EnsEMBL::Transcript') || 
       !$new->isa('Bio::EnsEMBL::Transcript')) {
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
   my ($self,$old,$new) = @_;

   if( !defined $old || !defined $new || ! ref $old || !$old->isa('Bio::EnsEMBL::Transcript') || !$new->isa('Bio::EnsEMBL::Transcript')) {
       $self->throw("Did not give me both old and new transcripts in increment Transcript, got: $old,$new");
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
