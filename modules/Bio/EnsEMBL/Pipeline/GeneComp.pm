
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
             Bio::EnsEMBL::Pipeline::GeneComp->map_temp_Exons_to_real_Exons($dbobj,
									    $timdb,
									    @newexons);


    # $mapped - reference to array of tempexons mapped to their new ids, with versions
    # and modified stamps correctly placed

    # $new - reference to array of new exons, with new ids, versions set to 1 and
    # modified/created time stamps.

    # $untransfered - reference to array of old exons, with old ids which although were 
    # remapped did not have exons in the new database to map to.


=head1 DESCRIPTION

This is a methods bag, not a real object. It deals with mapping exons, transcript and genes
from old versions through to new versions. This is where calls to get_new_ExonID etc are
actually made, and where the version logic happens. To do the mapping we need to get the
old exons out in the new coordinates (remapping). This currently is hidden behind the method
call get_old_Exons() on a contig object. This call returns old exon objects in the new coordinates,
with the method ->identical_dna set to true or not. 

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::GeneComp;

use strict;


=head2 map_temp_Exons_to_real_Exons

 Title   : map_temp_Exons_to_real_Exons
 Usage   : @mappedexons = Bio::EnsEMBL::Pipeline::GeneComp->map_temp_Exons_to_real_Exons($dbobj,$tim,@tempexons);
 Function:
 Example :
 Returns : exon objects with valid ids 
 Args    : database object to call get_new_ExonID()
           database object whoes contig have get_old_Exons();
           a list of temporary exons


=cut

sub map_temp_Exons_to_real_Exons{
   my ($dbobj,$timdb,@tempexons) = @_;

   if( !ref $dbobj || !$dbobj->isa('Bio::EnsEMBL::Pipeline::DB::ObjI') ) {
       die "This is **dreadful** I don't even have a database object to throw an exception on!";
   }

   if( !ref $timdb || !$timdb->isa('Bio::EnsEMBL::DB::ObjI') ) {
       $dbobj->throw("No second DB::ObjI provided - remember you need to provide two dbobjects!");
   }
   
   if( scalar(@tempexons) == 0 ) {
       $dbobj->warn("No temporary exons passed in for mapping. Can't map! - returning an empty list");
       return ();
   }

   # ok - we are ready to rock and roll...

   # we need some internal data structures during the mapping process.

   my %contig;   # hash of contig objects
   my %oldexons; # a hash of arrays of old exon objects
   my %moved;    # shows which old exons we have moved (or not).
   my %oldexonhash; # direct hash of old exons for final untransfered call
   my %ismapped; # tells us whether temporary exons have been mapped or not
   my @mapped;   # exons we have mapped.
   my @new;      # new exons
  
   # initialise ismapped to 0 for all temp exons and
   # build a hash of contig objects we want to see and exons.
   # better to do this is a separate loop as we have many exons to one contig

   foreach my $tempexon ( @tempexons ) {
       if( exists $ismapped{$tempexon->id} ) {
	   $dbobj->throw("you have given me temp exons with identical ids. bad...");
       }

       $ismapped{$tempexon->id} = 0;
       my $tempcid = $tempexon->contig_id;

       if( !exists $contig{$tempcid} ) {
	   # this will throw an exception if it can't find it
	   $contig{$tempcid} = $timdb->get_Contig($tempcid);
	   $oldexons{$tempcid} = [];
	   push(@{$oldexons{$tempcid}},$contig{$tempcid}->get_old_Exons);
	   
	   # set moved hash to zero
	   foreach my $oldexon ( @{$oldexons{$tempcid}} ) {
	       $moved{$oldexon->id} = 0;
	       $oldexonhash{$oldexon->id} = $oldexon;
	   }

       }
   }

   # get out one date for this mapping....

   my $time = time();

   
   # go over each exon and map old->new...

   TEMPEXON :

   foreach my $tempexon ( @tempexons ) {
       
       foreach my $oldexon ( @{$oldexons{$tempexon->contig_id}} ) {
       
	   # if start/end/strand is identical, it is definitely up for moving.
	   
	   if( $tempexon->start == $oldexon->start &&
	       $tempexon->end   == $oldexon->end   &&
	       $tempexon->strand == $oldexon->strand ) {

	       # ok - if we have already moved this - ERROR

	       if( $moved{$oldexon->id} == 1 ) {
		   $dbobj->throw("attempting to move old exon twice with identical start/end/strand. Not clever!");
	       }

	       # set ismapped to 1 for this tempexon

	       $ismapped{$tempexon->id} = 1;

	       # we move the id. Do we move the version?
	       $tempexon->id($oldexon->id);
	       push(@mapped,$tempexon);

	       if( $oldexon->has_identical_sequence == 1) {
		   # version and id
		   $tempexon->version($oldexon->version);
		   # don't update modified
	       } else {
		   $tempexon->version($oldexon->version()+1);
		   $tempexon->modified($time);
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
	   $tempexon->id($biggestoverlap->id);
	   $tempexon->version($biggestoverlap->version()+1);
	   $tempexon->modified($time);
	   push(@mapped,$tempexon);
	   
	   $moved{$biggestoverlap->id} = 1;
	   next TEMPEXON;
       } else {
	   # ok - new Exon
	   $tempexon->id($dbobj->gene_obj->get_new_ExonID);
	   $tempexon->created($time);
	   $tempexon->modified($time);
	   $tempexon->version(1);
	   push(@new,$tempexon);
	   next TEMPEXON;
       }
    
       $dbobj->throw("Error - should never reach here!");
   }
	      
   # find exons which have not moved, and push on to untransfered array
   
   my @untransfered;

   foreach my $oldexonid ( keys %moved ) {
       if( $moved{$oldexonid} == 0 ) {
	   push(@untransfered,$oldexonhash{$oldexonid});
       }
   }

   return (\@mapped,\@new,\@untransfered);

}


