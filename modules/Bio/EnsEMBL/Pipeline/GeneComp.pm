
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

    # $dbobj is the analysis database
    # $timdb is the tim database implementing get_old_Exons on the clone
    # (in the future $dbobj and $timdb could be the same object)
    # @newexons is a set of exons with temporary ids assigned
    # @mappedexons is the set of exons with olds ids, version numbers and new ids etc...

    # this module *does not* deal with writing the new, mapped, exons into the database.


    ($mapped,$new,$untransfered) = 
             Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Exons_to_real_Exons($dbobj,
									    $timdb,
									    @newexons);


    # $mapped - reference to array of tempexons mapped to their new ids, with versions
    # and modified stamps correctly placed

    # $new - reference to array of new exons, with new ids, versions set to 1 and
    # modified/created time stamps.

    # $untransfered - reference to array of old exons, with old ids which although were 
    # remapped did not have exons in the new database to map to.


=head1 DESCRIPTION

This is a methods bag, not a real object. It deals with mapping exons,
transcript and genes from old versions through to new versions. This
is where calls to get_new_ExonID etc are actually made, and where the
version logic happens. To do the mapping we need to get the old exons
out in the new coordinates (remapping). This currently is hidden
behind the method call get_old_Exons() on a contig object. This call
returns old exon objects in the new coordinates, with the method
->identical_dna set to true or not.

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

This is thick code. It is less thick than tims code. There are comments,
but probably not enough. I cant see many other areas for subroutines...


=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::GeneComp;

use strict;
use Carp;


=head2 map_temp_Exons_to_real_Exons

 Title   : map_temp_Exons_to_real_Exons
 Usage   : ($mapped,$new,$untransfered) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Exons_to_real_Exons($newvc,$logfile,@tempexons);
 Function:
 Example :
 Returns : exon objects with valid ids 
 Args    : New Virtual Contig, logfile filehandle


=cut

sub map_temp_Exons_to_real_Exons{
   my ($vc,$logfile) = @_;


   if( !ref $vc || !$vc->isa('Bio::EnsEMBL::Virtual::Contig') ) {
       die("No virtual contig provided - remember you need a virtual contig from the old !");
   }
   
   if( !defined $logfile ) {
       $vc->throw("Must pass logfile in for logging");
   }
   print STDERR "Getting all new exons...";
   my @tempexons=$vc->get_all_Exons();
   my $size=scalar(@tempexons);
   print STDERR " got $size new exons\n";

  

       
   if( scalar(@tempexons) == 0 ) {
       $vc->warn("No new exons found on this contig. Nothing to map - returning an empty list");
       return ();
   }

   # ok - we are ready to rock and roll...

   # we need some internal data structures during the mapping process.

   my %contig;   # hash of contig objects
   my %oldexons; # a hash of arrays of old exon objects
   my %moved;    # shows which old exons we have moved (or not).
   my %oldexonhash; # direct hash of old exons for final untransfered call
   my %ismapped; # tells us whether temporary exons have been mapped or not
   my %temp_old;
   my @new;      # new exons
   my @mapped;   #mapped exons
  
   # initialise ismapped to 0 for all temp exons and
   # build a hash of contig objects we want to see and exons.
   # better to do this is a separate loop as we have many exons to one contig

   foreach my $tempexon ( @tempexons ) {
       if( exists $ismapped{$tempexon->id} ) {
	   $vc->throw("you have given me temp exons with identical ids. bad...");
       }

       $ismapped{$tempexon->id} = 0;
   }
   print STDERR "Getting all old exons... ";
   my @oldexons=$vc->get_old_Exons($logfile);
   my $size=scalar(@oldexons);
   print STDERR "got $size old exons mapped\n";
   # get out one date for this mapping....

   my $time = time();

   
   # go over each exon and map old->new...

   TEMPEXON :

   foreach my $tempexon ( @tempexons ) {
       
       foreach my $oldexon ( @oldexons ) {
       
	   # if start/end/strand is identical, it is definitely up for moving.
	   
	   if( $tempexon->start == $oldexon->start &&
	       $tempexon->end   == $oldexon->end   &&
	       $tempexon->strand == $oldexon->strand ) {

	       # ok - if we have already moved this - ERROR

	       if( $moved{$oldexon->id} == 1 ) {
		   print $logfile "attempting to move old exon twice ".$oldexon->id." with identical start/end/strand. ".$tempexon->start.":".$tempexon->end.":".$tempexon->strand." Not clever!\n";
		   next;
	       }

	       # set ismapped to 1 for this tempexon

	       $ismapped{$tempexon->id} = 1;

	       $temp_old{$tempexon->id}=$oldexon->id;
	       # we move the id. Do we move the version?
	       $tempexon->id($oldexon->id);
	       print $logfile "EXON MAP IDENTICAL: ",$tempexon->id," ",$oldexon->id," "; # will add version
	       push (@mapped,$tempexon);
	       	       
	       if( $oldexon->seq eq $tempexon->seq) {
		   # version and id
		   $tempexon->version($oldexon->version);
		   print $logfile "Identical\n";
		   # don't update modified
	       } else {
		   $tempexon->version($oldexon->version()+1);
		   $tempexon->modified($time);
		   print $logfile "Version increment to ",$oldexon->version()+1,"\n";
	       }
	       
	       $moved{$oldexon->id} = 1;
	       
	       next TEMPEXON;
	   }
	   
       }
   



      # we can't map this exon directly. Second scan over oldexons to see whether
      # there is a significant overlap

      my $biggestoverlap = undef;
      my $overlapsize = 0;
      foreach my $oldexon ( @{$oldexons{$tempexon->contig_id}} ) {
	  
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

       # if I have got a biggest overlap - map across...
       if( defined $biggestoverlap ) {
	   # set ismapped to 1 for this tempexon
	   
	   $ismapped{$tempexon->id} = 1;
	   
	   # we move the id. Do we move the version?
	   print $logfile "EXON MAP BEST FIT: ",$tempexon->id," ",$biggestoverlap->id," ",$biggestoverlap->version()+1,"\n";

	   $temp_old{$tempexon->id}=$biggestoverlap->id;
	   $tempexon->id($biggestoverlap->id);
	   $tempexon->version($biggestoverlap->version()+1);
	   $tempexon->modified($time);
	   push(@mapped,$tempexon);
	   
	   $moved{$biggestoverlap->id} = 1;
	   next TEMPEXON;
       } else {
	   # ok - new Exon
	   my $tempid = $tempexon->id();
	   $tempexon->id($vc->dbobj->gene_Obj->get_new_ExonID);
	   $tempexon->created($time);
	   $tempexon->modified($time);
	   $tempexon->version(1);
	   print $logfile "EXON MAP NEW: ",$tempid," ",$tempexon->id," 1\n"; 

	   push(@new,$tempexon);
	   next TEMPEXON;
       }
    
       $vc->throw("Error - should never reach here!");
   }
	      
   # find exons which have not moved, and push on to untransfered array
   
   my @untransfered;

   foreach my $oldexonid ( keys %moved ) {
       if( $moved{$oldexonid} == 0 ) {
	   push(@untransfered,$oldexonhash{$oldexonid});
       }
   }

   return (%temp_old,\@mapped,\@new,\@untransfered);
}


=head2 map_temp_Genes_to_real_Genes

 Title   : map_temp_Genes_to_real_Genes
 Usage   : ($deadgeneid,$deadtranscriptid) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Genes_to_real_Genes($vc,$logfile,\@tempgenes,\@oldgenes);
 Function:
 Example :
 Returns : 
 Args    :

 NB you must have mapped exons first.

=cut

sub map_temp_Genes_to_real_Genes{
    my ($vc,$logfile,%temp_old) = @_;
    
    print STDERR "In genecomp->map_genes, got $vc, $logfile, %temp_old\n";
    
    if( !defined $vc) {
	die("Not passing the Virtual Contig to do gene mapping");
    }
    print STDERR "Getting all new genes... ";
    my @tempgenes = $vc->get_all_Genes();
    my $size=scalar @tempgenes;
    print STDERR "got $size new genes\n";

    if (%temp_old) {
	print STDERR "In gene mapping, started exon renaming\n";
	foreach my $gene (@tempgenes) {
	    foreach my $exon ($gene->each_unique_Exon) {
		if ($temp_old{$exon->id}) {
		    #print STDERR "In gene mapping, mapping ".$exon->id." to ".$temp_old{$exon->id}."\n";
		    $exon->id($temp_old{$exon->id});
		}
	    }
	}
    }

    print STDERR "Getting all old genes... ";
    my @oldgenes  = $vc->get_old_Genes();
    my $size=scalar(@oldgenes);
    print STDERR "got $size old genes\n";
    
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
    
    my %simple;  # key is old gene id, value is the new gene id to map to
    my %reversed_simple; # key is the new gene id, value is the old gene id
    
    my %split;   # key is old gene id, value an array of new gene id to map to
    my %merge;   # hash on the new gene id -> value is an array of old gene ids
    
    foreach my $og ( @oldgenes ) {
	
	my $currentgeneid = undef;

	foreach my $oe ( $og->each_unique_Exon ) {
	    print $logfile "Looking at ",$oe->id,"which is [",$newe2g{$oe->id},"]\n";
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
	    
	    print $logfile "Looking at $tgeneid for ".$og->id."\n";
	    
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
	   #print STDERR "Making it a simple!\n";

	   $simple{$og->id} = $tgeneid;
	   $reversed_simple{$tgeneid} = $og->id;
       }
   }


   # now we have a list of potential mappings. We step over these, mapping transcripts
   # from old to new.

   # simple is easy ;).

   my @dead_transcript_ids;
   my @dead_gene_ids;
   my $now = time();

   foreach my $oldgeneid ( keys %simple ) {

       my $newgeneid = $simple{$oldgeneid};

       #print STDERR "Doing $oldgeneid for $newgeneid\n";

       # flag that we have done this move before the ids change ;)

       $has_done_new{$newgeneid} = 1;
       $has_moved_old{$oldgeneid} = 1;


       my @newtrans = $newg{$newgeneid}->each_Transcript;
       my @oldtrans = $oldg{$oldgeneid}->each_Transcript;

       my ($should_increment,$additional_dead_transcript_ids) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Transcripts_to_real_Transcripts($vc,$logfile,\@newtrans,\@oldtrans);

       push(@dead_transcript_ids,@$additional_dead_transcript_ids);


       # deal with the mapping of ids

       print $logfile "GENE MAP: MOVED ",$newgeneid," ",$oldgeneid,"\n";

       $newg{$newgeneid}->id($oldgeneid);

       if( $should_increment ) {
	   $newg{$newgeneid}->version($oldg{$oldgeneid}->version+1);
	   $newg{$newgeneid}->modified($now);
       }
       #print STDERR "About to dump after move!\n";

       $newg{$newgeneid}->_dump(\*STDERR);
   }

   # merges are also quite easy.

   foreach my $newgeneid ( keys %merge ) {
       my @newtrans = $newg{$newgeneid}->each_Transcript;
       my @oldtrans;
       my $largest;
       my $size = 0;


       # flag that we have moved these
       $has_done_new{$newgeneid} = 1;

       foreach my $oldgeneid ( @{$merge{$newgeneid}} ) {
	   $has_moved_old{$oldgeneid} = 1;
	   if( $oldgeneid ne $largest ) {
	       push(@dead_gene_ids,$oldgeneid);
	   }
       }

       #print STDERR "Into merge code....\n";

       foreach my $oldgeneid ( @{$merge{$newgeneid}} ) {
	   my $tsize = scalar ( $oldg{$oldgeneid}->each_unique_Exon );

	   if( $tsize > $size ) {
	       $largest = $oldgeneid;
	   }

	   push(@oldtrans,$oldg{$oldgeneid}->each_Transcript);
       }
       
       my ($should_increment,$additional_dead_transcript_ids) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Transcripts_to_real_Transcripts($vc,$logfile,\@newtrans,\@oldtrans);

       push(@dead_transcript_ids,@$additional_dead_transcript_ids);


       # deal with the mapping of ids
       print $logfile "GENE MAP: MERGED ",$newgeneid," ",$oldg{$largest}->id,"\n";

       $newg{$newgeneid}->id($oldg{$largest}->id);

       # we increment irregardless of anything else
       $newg{$newgeneid}->version($oldg{$largest}->version+1);
       $newg{$newgeneid}->modified($now);
       
       
   }

   # splits **suck** big time

   foreach my $oldgeneid ( keys %split ) {
       my @newgeneid = @{$split{$oldgeneid}};

       # we take old transcripts one at a time, and fit to all possible
       # new transcripts. We take the best. The first case wins the gene id 
       #and the rest get assigned new geneids

       my $assigned = undef;
       

       my %fitted;
       foreach my $trans ( $oldg{$oldgeneid}->each_Transcript ) {
	   my $score = 0;
	   my $current_fit = undef;

	   # needs to up here so we can assign it later on.


	   # flag that we have moved these
	   $has_moved_old{$oldgeneid} =1;
	   
	   
	   foreach my $newgeneid ( @newgeneid ) {
	       $has_done_new{$newgeneid} = 1;
	   }



	   my $lastgeneid;
#	   print STDERR "Looking at splits\n";
	   GENE: foreach my $newgeneid ( @newgeneid ) {
	       #print STDERR "Going trough $newgeneid\n";
	       foreach my $newtrans ( $newg{$newgeneid}->each_Transcript ) { 

		   if( exists $fitted{$newtrans->id} ) {
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
	       # noone wants this transcript!
	       push(@dead_transcript_ids,$trans);
	       next;
	   }

	   
	   # current_fit is the best fit of this old transcript
	   
	   $current_fit->id($trans->id);
	   if( Bio::EnsEMBL::Pipeline::GeneComp::increment_Transcript($trans,$current_fit) == 1 ) {
	       $current_fit->version($trans->version()+1);
	   } else {
	       $current_fit->version($trans->version());
	   }
	   $fitted{$current_fit->id} = 1;
	   if( $assigned == 0 ) {
	       # this gene id wins. Hurray!
	       print $logfile "GENE MAP: SPLIT ",$lastgeneid," ",$oldgeneid,"\n";
	       #print STDERR "\n\nTrying to assign $oldgeneid to $lastgeneid\n\n";
	       $newg{$lastgeneid}->id($oldgeneid);
	       $newg{$lastgeneid}->version($oldg{$oldgeneid}->version+1);
	   }

       }
   
       # we need to create new transcripts for the remainder of the new transcripts not
       # fitted

       foreach my $newgeneid ( @newgeneid ) {
	   foreach my $newtrans ( $newg{$newgeneid}->each_Transcript ) { 
	       if( exists $fitted{$newtrans->id} ) {
		   next;
	       }
	       my $tempid = $newtrans->id; 
	       $newtrans->id($vc->dbobj->gene_Obj->get_new_TranscriptID);
	       print $logfile "TRANSCRIPT MAP: SPLIT ",$newgeneid," ",$oldgeneid,"\n";

	       $newtrans->version(1);
	       $newtrans->created($now);
	       $newtrans->modified($now);
	   }

	   if( $newg{$newgeneid}->id ne $oldgeneid ) {
	       my $newgene = $newg{$newgeneid};
	       my $tempid = $newgene->id();
	       # it is an unassigned gene...
	       $newgene->id($vc->dbobj->gene_Obj->get_new_GeneID());

	       print $logfile "GENE MAP: SPLIT NEW ",$tempid," ",$newgene->id,"\n";

	       $newgene->created($now);
	       $newgene->modified($now);
	       $newgene->version(1);
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

       $newgene->id($vc->dbobj->gene_Obj->get_new_GeneID());
       $newgene->created($now);
       $newgene->modified($now);
       $newgene->version(1);

       foreach my $t ( $newgene->each_Transcript ) {
	   $t->id($vc->dbobj->gene_Obj->get_new_GeneID());
	   $t->created($now);
	   $t->modified($now);
	   $t->version(1);
       }

   }

   foreach my $gene ( @oldgenes ) {
       if( $has_moved_old{$gene->id} ) {
	   next;
       }
       push(@dead_gene_ids,$gene);
   }

   return (\@dead_transcript_ids,\@dead_gene_ids);
}


=head2 map_temp_Transcripts_to_real_Transcripts

 Title   : map_temp_Transcripts_to_real_Transcripts
 Usage   : ($should_increment,$dead_transcript_array) = map_temp_Transcripts_to_real_Transcripts($dbobj,\@new,\@old)
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub map_temp_Transcripts_to_real_Transcripts{
   my ($vc,$logfile,$new,$old) = @_;

   my @newt = @$new;
   my @oldt = @$old;
   my $should_change = 0;
   my @dead;

   my %newt;

   foreach my $t ( @newt ) {
       $newt{$t->id} = $t;
   }
   

   # sort old by number of exons - largest first
   @oldt = sort { $b->number <=> $a->number } @oldt;

   my %fitted;
   my $now = time();

   # best fit to new transcripts...
   foreach my $oldt ( @oldt ) {
       my $score = 0;
       my $fitted = undef;

       #print STDERR "Looking at old transcript ",$oldt->id,"\n";

       foreach my $newt ( @newt ) {
	   if( exists $fitted{$newt->id} ) {
	       next;
	   }
	   my ($tscore,$perfect) = Bio::EnsEMBL::Pipeline::GeneComp::overlap_Transcript($oldt,$newt);
	   #print STDERR "Looking at newt ",$newt->id," score $tscore\n";

	   if( $perfect ) {
	       $fitted = $newt;
	       last;
	   }
	   if( $tscore > $score ) {
	       $fitted = $newt;
	   }
       }

       if ( defined $fitted ) {

	   print $logfile "TRANSCRIPT MAP : ",$fitted->id," ",$oldt->id,"\n";

	   $fitted->id($oldt->id);
	   $fitted{$fitted->id} = 1;
	   if( Bio::EnsEMBL::Pipeline::GeneComp::increment_Transcript($oldt,$fitted) == 1 ) {
	       $fitted->version($oldt->version()+1);
	       $fitted->modified($now);
	       $should_change = 1;
	   } else {
	       $fitted->version($oldt->version());
	   }
       } else {
	   # it is dead ;)
	   $should_change = 1;
	   push(@dead,$oldt->id);
	   print $logfile "TRANSCRIPT MAP: KILLED ",$oldt->id,"\n";
       }
   }

   foreach my $newt ( @newt ) {
       if( exists $fitted{$newt->id} ) {
	   next;
       }
       my $tempid = $newt->id;
       $should_change = 1;
       $newt->id($vc->dbobj->gene_Obj->get_new_TranscriptID);
       $newt->version(1);
       $newt->created($now);
       $newt->modified($now);
       print $logfile "TRANSCRIPT MAP: NEW ",$tempid," ",$newt->id,"\n";
   }

   return ($should_change,\@dead);

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


1;
