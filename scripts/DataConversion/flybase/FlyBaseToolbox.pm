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

# methods to export
@EXPORT =qw(addTranslationObject convert_to_exon_coords setExonPhases);

$|=1;


#
# alternative tools / methods
#
 ################################################################################


sub addTranslationObject{
  my ($cds_id,$nw_transcript,$cds_start,$cds_end,$cds_strand) = @_;

  my $seq_start = -1;
  my $seq_end = -1 ;
  my ($first_coding_exon, $last_coding_exon);

#  print "processing $cds_id\n";

  # getting first and last coding exon
  #################################################

  for my $ex (@{$nw_transcript->get_all_Exons()}) {


    # forward strand
    ####################################

    if($nw_transcript->strand == 1){
      if ($ex->start <= $cds_start  && $cds_start <= $ex->end) {
        # EXON CONTAINS TSS, coord-conversion
        $seq_start = convert_to_exon_coords($cds_start, $ex->start,$ex->end,$ex->strand );
        $first_coding_exon = $ex;

        # set start-phase of first coding exon
        $first_coding_exon->phase("0");
      }
      if ($ex->start <= $cds_end  && $cds_end <= $ex->end) {
        # EXON CONTAINS TES, coord-conversion and add the stop-codon (3bp.) to the seq-end
        $seq_end = convert_to_exon_coords($cds_end, $ex->start,$ex->end,$ex->strand ) + 3 ;
        $last_coding_exon = $ex;
      }
    }

    # reverse strand
    #####################################
    if ($nw_transcript->strand eq "-1"){

      if( ($ex->start <= $cds_end) && ($cds_end <= $ex->end)){
        print "first coding exon is ".$ex->stable_id." with start ".$ex->start ."\n";
        $first_coding_exon = $ex;
        $seq_start = ($ex->end - $cds_end) + 1 ;
      }

      if ( ($ex->start <= $cds_start) &&  ($cds_start <= $ex->end)){
        $last_coding_exon = $ex;
        $seq_end = ($last_coding_exon->end - $cds_start)+ 1 + 3 ;
      }
    }
  }

# D E B U G
#  print "\n"x2;print "----translation_data of $cds_id ------\n";
#  print "$cds_id\t".$nw_transcript->strand."\t$cds_start\t$cds_end\tFirst-coding-exon-start:".$first_coding_exon->start."\tfce-offset: ".$seq_start."\n";
#  print "$cds_id\t".$nw_transcript->strand."\t$cds_start\t$cds_end\tlast-coding-exon-start:".$last_coding_exon->start."\tLce-offset: ".$seq_end."\n";
#  print "-------XXXXXXXXXXXXXXXXXXXXX------\n";



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
  unless ($first_coding_exon && $last_coding_exon ) {
    throw("There is a CDS-Start or CDS-End, but I can't find an exon which spans these points\n");
  }

  if ($first_coding_exon eq $last_coding_exon){
    my $tl_end =  $nw_transcript->translation->end;
    if($first_coding_exon->length - $tl_end > 0){
      # single-exon gene with UTR
      $first_coding_exon->end_phase("-1");
    }else{
      my $endphase = ($first_coding_exon->length - $seq_start) % 3 ;
      $first_coding_exon->end_phase($endphase);
    }
  }



  # consistency checks
  throw("There seems to be no cds-region-start in CDS-ID $cds_id\n")  if ($seq_start  eq "-1");
  throw("There seems to be no cds-region-end in CDS-ID $cds_id\n")   if ($seq_end  eq "-1") ;

  return $nw_transcript;
}



sub setExonPhases{
  my ( $self,$transcript) = @_;

  my $fce = $transcript->translation->start_Exon;
  my $lce = $transcript->translation->end_Exon;

  my @all_exons = @{$transcript->get_all_Exons};

  my $end_phase;

  ALL_EXON: for my $exon (@all_exons){

      # plus-strand
      if($exon->strand == 1 ){

        if ($exon->start < $fce->start && $exon->end < $fce->end){
          # 5'-complete UTR
          $exon->phase("-1");
          $exon->end_phase("-1");
        }

        if ($exon->start >= $fce->start && $exon->end <=$lce->end){
          # exon is fce
          if($exon eq $fce){
            print $exon->stable_id()." : first_coding_exon\t";
            # set startphase of first exon
            if($self->processed_exon($exon)){
              print "processed\n";
             $exon->phase($self->processed_exon($exon)->phase);
            }else{
              print "not_processed\t";
              $exon->phase( 0 ) ;
            }
            # set end-phase
            print "tss:".$transcript->translation->start."\t";
            print "elength:".$exon->length."\t";
            $end_phase = ($exon->length - $transcript->translation->start + 1 )%3;
            $exon->end_phase($end_phase);
            print "sp: 0\tep: $end_phase\n";
          }


          # exon is middle exon
          if ( $exon ne $fce  && $exon ne $lce){
            print $exon->stable_id()." : middle_exon\t";
            $exon->phase($end_phase);
            $end_phase =( $exon->phase + $exon->length ) % 3;
            $exon->end_phase($end_phase);
            print "sp: ".$exon->phase . "\tep: $end_phase\n";
          }

          # exon is lce
          if($exon eq $lce){
            print $exon->stable_id()." : last_coding_exon\t";
            # only change phase if single-exon-gene
            if ($exon ne $fce){
              $exon->phase($end_phase);
            }

            # check if exon contains UTR-tailing end
            if ( ($exon->length - $transcript->translation->end  ) == 0){
              $exon->end_phase(  ($end_phase + $exon->length) % 3);
            }else{
              print "3'UTR";
              $exon->end_phase("-1");
            }
          }
          if ($exon->start > $lce->start){
            $exon->end_phase("-1");
            $exon->phase("-1");
          }
        }
        # end PLUS
      }else{
        # start REVSTRAND
        if ($exon eq $fce){
          #          print "exon equals fce\t".$exon->stable_id . "\t".$fce->stable_id."\n";

          # set start-phase of first coding exon
          if($self->processed_exon($exon)){
            # exon is also member of other transcript->use the already calculated startphase
            $exon->phase($self->processed_exon($exon)->phase);
          }else{
            $exon->phase( 0 ) ;
          }

          $end_phase = ($exon->length - ($transcript->translation->start -1) ) % 3;
          $exon->end_phase($end_phase);
          #          print "change:".$exon->stable_id."\tep:".$exon->end_phase."\n";

          # consistency-check: check endphase of "new" exon against end-phase of known exon
          if ($self->processed_exon($exon)->end_phase() ne $end_phase){
            warn("mismatch between phases of exons used by different transcripts\t".$exon->stable_id." ep:".$self->processed_exon($exon)->end_phase."vs $end_phase\n");
          }
        }

        # middle exon
        if($exon ne $fce && $exon ne $lce){
          $exon->phase($end_phase);
          $end_phase = ($exon->phase + $exon->length) % 3;
          $exon->end_phase($end_phase);
        }
        # exon is lce
        if ($exon eq $lce){
          $exon->phase($end_phase);

          # check if exon contains UTR-tailing end
          if ( ($exon->length - $transcript->translation->end  ) == 0){
            $exon->end_phase(  ($end_phase + $exon->length) % 3);
          }else{
            # exon has 3' UTR-end
            $exon->end_phase("-1");
          }
        }

        # processing the UTR's:
        ########################################

        if ($exon->start > $fce->start){
          # 5' UTR
          $exon->phase("-1");
          $exon->end_phase("-1");
        }

        # 3'UTR
        if($exon->start < $lce->start){
          $exon->phase("-1");
          $exon->end_phase("-1");
        }
      } #end REVSTRAND



      $self->processed_exon($exon);
    }




  print "\n----summary of calcualted values for transcript ".$transcript->stable_id."-------------------------\n";
  for( @{$transcript->get_all_Exons}){
    print "EXON:".$_->stable_id ."\t". $_->hashkey . "\tSphase:".$_->phase."\tEphase:".$_->end_phase."\n";
  }
  print "\n"x4;
  return $transcript;
}


 # ALL_EXON
#  }

#      if($fce eq $lce){
#        $prev_end_phase = $exon->end_phase;
#        $exon->phase($prev_end_phase);
#        last ALL_EXON;

#      }elsif($exon eq $fce){
#        # first coding exon
##        print "startexon set phase to $prev_end_phase\n";
#        $prev_end_phase = $exon->end_phase;

#      }elsif( $exon eq $lce){
#        # last coding exon
##        print "LCE\n";
#        $exon->phase($prev_end_phase);

#        # end-phase -1 if UTR
#        if( abs( $exon->length - ($exon->end - $exon->start)) ne 0 ){
#          $exon->end_phase(-1);
#        }else{
#          $exon->end_phase( ($exon->phase + $exon->length) %3 );
#        }
##        print "lastexon ".$exon->phase."\t".$exon->end_phase."\n";
#        last ALL_EXON;
#      }else{
#        # framend exon
#        $exon->phase($prev_end_phase) ;
#        $exon->end_phase( ($exon->phase  + $exon->length) % 3  );
#        $prev_end_phase=$exon->end_phase;
##        print "middle ".$exon->phase ."\t$prev_end_phase\n";
#      }
#    }








=pod

=head3 alternative methods

=head2 convert_to_exon_coords

  Usage       : convert_to_exon_coords ($cds_coord, $exon_start, $exon_end )
  Function    : checks if the cds-coordiante (start or end) lies in an exon
  Return-Value: returns the coordiante of the TSS in exon_coordiantes if given coordinate is inside the exon, and  "-1" if cds_coord is not
=cut

sub convert_to_exon_coords{
  my ( $cds_coord, $exon_start, $exon_end, $exon_strand ) = @_;
  my $cds_coord_converted = "-1";
  warn_inconsistency("Exon-start-coord >= exon_end_coord: $exon_start\t$exon_end\n")  if($exon_start >= $exon_end);

  # check if cds-coordinate is inside exon
  if ($cds_coord >= $exon_start){
    if ($cds_coord <= $exon_end){

      # cds_coord is in exon, convert coordiantes to exon-coords
      if ($exon_strand eq "+1"){
        $cds_coord_converted = $cds_coord - $exon_start +1 ;

      }elsif($exon_strand eq "-1"){
        # on rev-strand, $exon_end is the start and $exon_start is end
#        $cds_coord_converted = $exon_end - $cds_coord + 1;
        # changend
        $cds_coord_converted = $cds_coord - $exon_start +1 ;
      }else{
        # exon is on unkown strand... what now ?
        warn_inconsistency("Exon is on unknown strand, coords: $exon_start:$exon_end:$exon_strand\n");
        $cds_coord_converted =-1;
      }

    }
  }
  throw ("coord-error\n")  if ($cds_coord_converted eq "0");
  return $cds_coord_converted;
}
