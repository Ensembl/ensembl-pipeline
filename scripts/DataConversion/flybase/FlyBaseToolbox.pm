#!/usr/local/ensembl/bin/perl  -w

package FlyBaseToolbox;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate verbose);
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Exporter;
our(@EXPORT,@ISA);
@ISA = ('Exporter');

@EXPORT =qw(add_TranslationObject convert_to_exon_coords setExonPhases);

$|=1;


# alternative tools / methods
#############################


sub add_TranslationObject{
  my ($cds_id,$nw_transcript,$cds_start,$cds_end,$cds_strand) = @_;

  my $seq_start = -1;
  my $seq_end = -1 ;
  my ($first_coding_exon, $last_coding_exon);

  # getting first and last coding exon
  #################################################

  my $i=-1;
  my @all_transcript_exons = @{$nw_transcript->get_all_Exons() } ;
  for my $ex (@all_transcript_exons) {
    $i++;

    #print "$cds_end\tHAVING EXON NR $i: " ; 
    #printObject($ex);

    # forward strand
    ####################################

    if ($nw_transcript->strand == 1) {
      if ($ex->start <= $cds_start  && $cds_start <= $ex->end) {
        # EXON CONTAINS TSS, coord-conversion
        $seq_start = convert_to_exon_coords($cds_start, $ex->start,$ex->end,$ex->strand );
        $first_coding_exon = $ex;
        # set start-phase of first coding exon
        $first_coding_exon->phase("0");
      }

      if ($ex->start <= $cds_end  && $cds_end <= $ex->end){

        # CDS ENDS IN THIS EXON
	# check if STOP-CODON (CDS+3) is also in this EXON or if we have to do  a jump to the next exon 

	# cds = 100 
	#
	#   50           100                       455              674
	#   |             |                         |                |
	#   +++++++++++++++                         ++++++++++++++++++
	#               CCC                         TGA
	#   cds ends at 100, but stop-codon is in next exon at 455 ! so translationobject has to be adjusted !!!!
	#

        $seq_end = convert_to_exon_coords($cds_end, $ex->start,$ex->end,$ex->strand ) + 3 ;
        $last_coding_exon = $ex;

	if($seq_end > $ex->length ){
	  #
	  # the ATG-stopcodon lies outside this exon, we have to shift the translation
	  # to the next exon

	  my $new_ex = $all_transcript_exons[$i+1];
	  if($new_ex){
	    $seq_end = $seq_end - $ex->length ;
	    ##$seq_end = convert_to_exon_coords($cds_end, $new_ex->start,$new_ex->end,$new_ex->strand ) + 3 ;
	    $last_coding_exon = $new_ex;
	    print STDERR "PLUS: translation jumped to next exon for plus-strand: seq_end-now: $seq_end\n";
	  }else{
	    print STDERR "WARN : cant extend to stop-codon, using flybase-cds-end instead\n" ;
	    $seq_end = convert_to_exon_coords($cds_end, $ex->start,$ex->end,$ex->strand ) ;
	    if($seq_end>$ex->length){
	      throw("translation-end outside of exon");
	    }
	  }
	}
      }
    }

    # reverse strand
    #####################################
    if ($nw_transcript->strand eq "-1") {

      #
      # adjust translation-START
      #
      if ( ($ex->start <= $cds_end) && ($cds_end <= $ex->end)) {
        $first_coding_exon = $ex;
        $seq_start = ($ex->end - $cds_end) + 1 ;
      }
      #
      # adjust translation-END (watch out for stop-codon !)
      #
      if ( ($ex->start <= $cds_start) &&  ($cds_start <= $ex->end)) {
        $last_coding_exon = $ex;
        $seq_end = ($last_coding_exon->end - $cds_start)+ 1 + 3 ;

	if($seq_end > $ex->length){
	  my $new_ex = $all_transcript_exons[$i+1];
	  if($new_ex){
	    $seq_end =  $seq_end -$ex->length;
	    #$seq_end = convert_to_exon_coords($cds_end, $ex->start,$ex->end,$ex->strand ) + 3 ;
	    $last_coding_exon = $new_ex;
	    print STDERR "translation jumped to next exon for MINUS-strand. SEQ-END now:$seq_end\n";
	  }else{
	    warn("WARN: exon ends without stop-codon $i and there is no next exon which contains the stop-codon\n" ) ;
	    print STDERR "WARN: SUSPICIOUS : CDS ends in exon, but exon is too short to add 3 bases for STOP-codon \n"; 
	    print STDERR "WARN:  and it's not possible to jump over into next exon because this is alreayd the last exon\n";
	    print STDERR "WARN: " .  printObject($nw_transcript) ;


	    $seq_end = ($last_coding_exon->end - $cds_start)+ 1 ;
	   if($seq_end>$ex->length){
	     throw("translation-end outside of exon");
	   }
	  }
	}
      }
    }
  }


  my $tl = Bio::EnsEMBL::Translation->new(
                                          -START_EXON => $first_coding_exon,
                                          -END_EXON   => $last_coding_exon,
                                          -SEQ_START  => $seq_start,
                                          -SEQ_END    => $seq_end,
                                          -STABLE_ID => $cds_id,
                                          -VERSION => "3",
                                         );
  # add translation
  $nw_transcript->translation($tl);

  # set end-phase for case if single exon contains 2 UTR-regions at end
  throw("There is a CDS-Start or CDS-End, but I can't find an exon which spans these points\n")   unless ($first_coding_exon && $last_coding_exon );

  if ($first_coding_exon eq $last_coding_exon) {
    my $tl_end =  $nw_transcript->translation->end;
    if ($first_coding_exon->length - $tl_end > 0) {
      # single-exon gene with UTR
      $first_coding_exon->end_phase("-1");
    } else {
      my $endphase = ($first_coding_exon->length - $seq_start) % 3 ;
      $first_coding_exon->end_phase($endphase);
    }
  }
  # consistency checks
  throw("There seems to be no cds-region-start in CDS-ID $cds_id\n")  if ($seq_start  eq "-1");
  throw("There seems to be no cds-region-end in CDS-ID $cds_id\n")   if ($seq_end  eq "-1") ;

  return $nw_transcript;
}

sub printObject{
  my ($o)= @_;
  print "START: " . $o->seq_region_start . 
         "\tEND: " . $o->seq_region_end . 
         "\tSTRD: " . $o->seq_region_strand . "\n";
}



sub setExonPhases{
  my ( $self,$transcript) = @_;
  my $fce = $transcript->translation->start_Exon;
  my $lce = $transcript->translation->end_Exon;
  my @all_exons = @{$transcript->get_all_Exons};
  my $end_phase;

  for my $exon (@all_exons) {
      ##### start processing forward strand #####
      if ($exon->strand == 1 ) {

        if ($exon->start < $fce->start && $exon->end < $fce->end) {
          # 5'-complete UTR
          $exon->phase("-1");
          $exon->end_phase("-1");
        }

        if ($exon->start >= $fce->start && $exon->end <=$lce->end) {

          # exon is fce
          if ($exon eq $fce) {
            # set startphase of first exon
            $exon->phase( 0 ) ;
            # set end-phase
            $end_phase = ($exon->length - $transcript->translation->start + 1 )%3;
            $exon->end_phase($end_phase);
          }

          # exon is middle exon
          if ( $exon ne $fce  && $exon ne $lce) {
            $exon->phase($end_phase);
            $end_phase =( $exon->phase + $exon->length ) % 3;
            $exon->end_phase($end_phase);
          }

          # exon is lce
          if ($exon eq $lce) {
            # only change phase if single-exon-gene
            if ($exon ne $fce) {
              $exon->phase($end_phase);
            }

            # check if exon contains UTR-tailing end
            if ( ($exon->length - $transcript->translation->end  ) == 0) {
              $exon->end_phase(  ($end_phase + $exon->length) % 3);
            } else {
              # 3'UTR
              $exon->end_phase("-1");
            }
          }
          if ($exon->start > $lce->start) {
            $exon->end_phase("-1");
            $exon->phase("-1");
          }
        }
      } else {

        ##### start processing reverse strand #####

        # first coding exon
        if ($exon eq $fce) {
          $exon->phase( 0 ) ;
          $end_phase = ($exon->length - ($transcript->translation->start -1) ) % 3;
          $exon->end_phase($end_phase);
        }

        # middle exon
        if ($exon ne $fce && $exon ne $lce) {
          $exon->phase($end_phase);
          $end_phase = ($exon->phase + $exon->length) % 3;
          $exon->end_phase($end_phase);
        }

        # exon is lce
        if ($exon eq $lce) {
          if ($exon ne $fce) {
            $exon->phase($end_phase);
          }

          # check if exon contains UTR-tailing end
          if ( ($exon->length - $transcript->translation->end  ) == 0) {
            $exon->end_phase(  ($end_phase + $exon->length) % 3);
          } else {
            # exon has 3' UTR-end
            $exon->end_phase("-1");
          }
        }

        # now the UTR's
        if ($exon->start > $fce->start) {
          # 5' UTR
          $exon->phase("-1");
          $exon->end_phase("-1");
        }

        # 3'UTR
        if ($exon->start < $lce->start) {
          $exon->phase("-1");
          $exon->end_phase("-1");
        }
      } #end REVSTRAND
    } #allexons

  return $transcript;
}




#=pod

#=head3 alternative methods

#=head2 add_DBEntry

#  Usage       : adds a DBEntry to the object
#  Function    :
#  Return-Value:
#=cut

#sub add_DBEntry {
#  my ($self,$obj) = @_;

#  my $dbentry = new Bio::EnsEMBL::DBEntry(
#                                          -adaptor => $adaptor,
#                                          -primary_id => $pid,
#                                          -dbname  => $dbname,
#                                          -release => $release,
#                                          -display_id => $did,
#                                          -description => $description
#                                         );
#}






=pod

=head2 convert_to_exon_coords

  Usage       : convert_to_exon_coords ($cds_coord, $exon_start, $exon_end )
  Function    : checks if the cds-coordiante (start or end) lies in an exon
  Return-Value: returns the coordiante of the TSS in exon_coordiantes if given coordinate is inside the exon, and  "-1" if cds_coord is not
=cut


sub convert_to_exon_coords{
  my ( $cds_coord, $exon_start, $exon_end, $exon_strand ) = @_;
  my $cds_coord_converted = "-1";
  warn_inconsistency("Exon-start-coord >= exon_end_coord: $exon_start\t$exon_end\n")  if ($exon_start >= $exon_end);

  print STDERR "\nCDS: cds_coord: $cds_coord :: Exon-start:  $exon_start ::: $exon_end ::: ".($exon_end-$exon_start)."\n" ;
      # cds_coord is in exon, convert coordiantes to exon-coords

  # check if cds-coordinate is inside exon
  if ($cds_coord >= $exon_start) {
    if ($cds_coord <= $exon_end) {
      #print STDERR "\nCDS: cds_coord: $cds_coord >= $exon_start AND cds_coord: $cds_coord <= $exon_end \n" ; 
      # cds_coord is in exon, convert coordiantes to exon-coords

      if ($exon_strand eq "+1") {
        $cds_coord_converted = $cds_coord - $exon_start + 1 ;

      } elsif ( $exon_strand eq "-1") {
        # on rev-strand, $exon_end is the start and $exon_start is end
        $cds_coord_converted = $cds_coord - $exon_start + 1 ;

      } else {
        # exon is on unkown strand... what now ?
        warn_inconsistency("Exon is on unknown strand, coords: $exon_start:$exon_end:$exon_strand\n");
        $cds_coord_converted =-1;
      }

    }
  }

  throw ("coord-error\n")  if ($cds_coord_converted eq "0");
  throw ("coord-error\n")  if ($cds_coord_converted eq "-1");
  return $cds_coord_converted;
}
