#!/usr/local/ensembl/bin/perl  -w

# This module currently contains a parser does not parse all information possible out of the file.

package FlyBaseGff;

use strict;
use warnings;
use FlyBaseConf;
use FlyBaseToolbox;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate verbose);
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;

use Data::Dumper;

$|=1;

=pod

=head2 new 

  Title       : new 
  Usage       : 
  Function    : new GFF object 
  Arguments   : 
  Return-Val  :

=cut

# construct newFlyBaseGff-object (GFF Version 3) of flybase gff
# desc of gff3: http://song.sourceforge.net/gff3.shtml
sub new  {
  my ($class,@args) = @_;
  my $self={};
  bless $self,$class;

  my ($db,$gff,$slice) = rearrange(
                            [qw(
                                DB
                                GFF
                                SLICE
                               )],@args);

  $self->db($db);
  $self->gff($gff);
  $self->slice($slice);
  $self->{_gene}=[];
  $self->{exons}={};
  $self->{db_entry}={};
  $self->parse_gff();
  print "\n* * *\ngff $gff parsed\n* * * \n";
  return $self;
}

=pod

=head2 db_entry 

  Title       : db_entry 
  Usage       :
  Function    : for xrefs 
  Arguments   :
  Return-Val  :

=cut
sub db_entry{
  my ($self, $tr,$db_entry) = @_;
  if($db_entry){
    $ {$self->{db_entry}} {$tr}=$db_entry;
  }
  return $ { $self->{db_entry}}{$tr};
}


=pod


=head2 get_all_attributes_by_id

  Title       : get_all_attributes_by_id
  Usage       : $obj->get_all_attributes_by_id('ID')
  Function    :  access all stored attributes of a feature-id
  Arguments   : ID
  Return-Val  : arrayref to arrayrefs of attributes of a given ID

=cut


sub get_childrenID_from_parentID {
   my ($self,$id) = @_;
   my $children_ref = {}; 
   if (${$self->parent2child}{$id}) {
     $children_ref = ${$self->parent2child}{$id};
   }
   return $children_ref;
}

sub get_features_from_ID {
   my ($self,$id) = @_;
   return ${$self->featurehash}{$id};
}



#
# Methods used for storing attributes, genes, SimpleFeatures
 ################################################################################



=pod

=head3 storing attributes, genes, SimpleFeatures

=head2 store_as_gene_object

  Title       :  store_as_gene_object
  Usage       :  store_as_gene_object ('gene', 'mRNA','flybase_gene')
  Function    :  looks in gff for entity of type 'gene' which has childs 'mRNA' and 'exon' and stores them 
                 as gene of type 'flybase_gene' in the database
  Decsription :  Given a slice, parent type (eg. Gene), child type (eg. mRNA, tRNA, pseudogene) 
                 this goes through the gff hash and finds all parents of the given type and 
                 makes them into Genes. If they have associates Transcripts and Exons then 
                 these are added too. They are all stored in the database.
  Arguments   :
  Return-Val  :

=cut

sub store_as_gene_object {
  my ($self, $parent_type, $child_type, $ensembl_gene_type) = @_;
  my %all_processed_exons;

  print "Parent feature_type=$parent_type and Child_feature_type=$child_type\n";

  # test if $parent_type is annotated in gff
  my @gene_ids = @{$self->get_ids_by_type($parent_type)};
  if (!@gene_ids || scalar(@gene_ids) < 1) {
    warning("No parent type '$parent_type' found in gff file ".$self->gff()." (child type $child_type)\n");
    return;
  } else {   
  # Get all gene ids for this type (gene) 
    print "Have ".scalar(@gene_ids)." genes of parent type $parent_type\n";
  }

  # want to timestamp the stable_ids
  my $time = time();

 GENE:  foreach my $gene_id ( @gene_ids ) {
    print "Looking at gene $gene_id\n";
    # note that genes with no transcripts are not stored
    # eg. mt:ori gene     

    # # #
    # Transcripts
    # # #
    # start looking for transcripts of this parent gene
    my $all_new_transcripts;
    # get a list of child ids for the parent
    # ie. get a list of transcript IDs for the Gene
    my @transcript_ids = keys %{$self->get_childrenID_from_parentID($gene_id)};
    my @keep;
    if (scalar(@transcript_ids) == 0) {
      # there are no children of this type for the parent
      print "   x WARNING: No child of type $child_type for $parent_type with id $gene_id\n";
      # make a gene without transcripts?
      next GENE;
    } else {
      # the parent has children of type $child_type
      # Loop thru the alternate transcripts
      foreach my $t (@transcript_ids) {
        if (${$self->get_childrenID_from_parentID($gene_id)}{$t}[1] eq $child_type) {
          push @keep, $t;
        }
      }
      if (scalar(@keep) == 0) {
        print "   x WARNING: No child of type $child_type for $parent_type with id $gene_id\n";
        next GENE;
      } else {
        print "\n* * *\n>> Doing gene $gene_id \n* * *\n";
        $all_new_transcripts = make_Transcripts($self, \@keep, $child_type, $time);      
      }
    } # if (scalar(@transcript_ids) == 0) {  

    print "Renaming exons\n";
    # Rename exons where necessary
    my ($rt) = rename_exons($gene_id, $all_new_transcripts);
    my @renamed_transcripts = @{$rt};
      
    # # #
    # Make Genes
    # # #
    # got all alternative transcripts and make genes
    my $gene = Bio::EnsEMBL::Gene->new(
                                     -START     => ${$self->get_features_from_ID($gene_id)}{'start'}, 
                                     -END       => ${$self->get_features_from_ID($gene_id)}{'end'}, 
                                     -STRAND    => ${$self->get_features_from_ID($gene_id)}{'strand'}, 
                                     -SLICE     => $self->slice,
                                     -ANALYSIS  => $self->create_analysis($LOGIC_NAME_EXON),
                                     -STABLE_ID => $gene_id,
                                     -VERSION   => $ATTRIBUTE_VERSION,
	                             -TYPE      => $ensembl_gene_type,
                                     -CREATED_DATE => $time,
                                     -MODIFIED_DATE => $time,
                                     -SOURCE    => ${$self->get_features_from_ID($gene_id)}{'source'}, 
                                    );

    # add Transcripts to Gene
    foreach my $tr (@renamed_transcripts) {
      $gene->add_Transcript($tr);
    }

    # # #
    # Store
    # # #
    # check if we already have an exon-stable-id/gene-stable-id/translation-stable-id/transcript-stable-id for this gene
    print "\nNow trying to store gene $gene_id....." ;
    my $gene_DB_id = $self->db->get_GeneAdaptor->store($gene);
    print "STORED gene $gene_DB_id \n\n" ;
    
#    # Xrefs 
#    foreach my $dbe (@{$self->get_dbxrefs($gene_id, $parent_type)}) {
#      $self->db->get_DBEntryAdaptor->store($dbe,$gene_DB_id, 'Gene');
#      print "  *stored xref with dbname ".$dbe->dbname." for ".ref($gene)." with id $gene_DB_id\n";
#    }

    # reprocess the db_entry for transcripts
    #for my $tr (@renamed_transcripts){
#    for my $tr (@{$gene->get_all_Transcripts}) {
#      foreach my $dbe (@{$self->get_dbxrefs($tr->stable_id, $child_type)}) {
#       # my $dbe = $self->db_entry($tr);
#        if($dbe){
#          $self->db->get_DBEntryAdaptor->store($dbe, $tr->dbID,'Transcript');
#          print "  *stored xref with dbname ".$dbe->dbname." for ".ref($tr)." with id ".$tr->dbID."\n";
#        }
#      }
#      if ( $tr->translation() ) {
#        my $tln = $tr->translation() ;
#        # we have a protein-coding gene and can gets xrefs for the protein
#        foreach my $dbe (@{$self->get_dbxrefs($tln->stable_id, 'protein')}) {
#          if ($dbe) {
#            $self->db->get_DBEntryAdaptor->store($dbe, $tln->dbID,'Translation');
#            print "  *stored xref with dbname ".$dbe->dbname." for ".ref($tln)." with id ".$tln->dbID."\n";
#          }
#        }
#      }
#    }

  } # for my $gene_id ( @{$self->get_ids_by_type($parent_type)} ) {
} # end sub store_as_gene_object {

sub make_Transcripts {
  my ($self, $tref, $biotype, $time) = @_;
  my @transcript_ids = @{$tref};
  my @transcripts;

  for my $transcript_id (@transcript_ids) {
    # make an EnsEMBL Transcript
    my $new_transcript = new Bio::EnsEMBL::Transcript(
                                                       -STABLE_ID => $transcript_id,
                                                       -VERSION => $ATTRIBUTE_VERSION,
                                                       -TYPE => $biotype, 
                                                       -CREATED_DATE => $time,
                                                       -MODIFIED_DATE => $time,
                                                       -ANALYSIS  => $self->create_analysis($LOGIC_NAME_EXON),
                                                      );
#    print "   Have new Transcript $transcript_id of biotype $biotype. Getting xrefs...\n";
#
#    # It takes the first Dbxref entry and stores it
#    # my @dbx = @{${$self->featurehash}{$transcript_id}{'Dbxref'}};
#    my @dbxrefs = @{$self->get_features_from_ID($transcript_id)->{'Dbxref'}};
#    if (scalar @dbxrefs) {
#      print "Found xrefs @dbxrefs";
#      # not sure if we should use all of them or just the first
#      $self->db_entry($new_transcript, ${$self->get_features_from_ID($transcript_id)}{'Dbxref'}[0]);
#    } else {
#      print "WARNING: No db entry transcript $transcript_id\n";
#    }

    # # #
    # Exons
    # # #
    my @exon_ids = keys %{$self->get_childrenID_from_parentID($transcript_id)};
    if (scalar(@exon_ids) == 0) {
      print "\n   x No children for transcript $transcript_id, must be non-coding\n";
    } else {
      print "\nTranscript $transcript_id has ".@exon_ids." possible  exons\n";
      my $transcript = make_Exons($self, $new_transcript, \@exon_ids, $time);
      my $transcript_with_tln = add_Translation($self, $transcript);  
      push @transcripts, $transcript_with_tln;
    }
  } # for my $transcript_id (@transcript_ids) {
  return \@transcripts;
}

sub make_Exons {
  my ($self, $transcript, $exon_ids, $time ) = @_;
  my @exons;

  # Loop thru exons and Make EnsEMBL objects
  EXON: foreach my $exon_id (@$exon_ids) {
    if (${$self->get_features_from_ID($exon_id)}{'type'} ne 'exon') {
      next EXON;
    }
    my $ex = Bio::EnsEMBL::Exon->new(
                                 -START     => ${$self->get_features_from_ID($exon_id)}{'start'}, 
                                 -END       => ${$self->get_features_from_ID($exon_id)}{'end'}, 
                                 -STRAND    => ${$self->get_features_from_ID($exon_id)}{'strand'}, 
                                 -PHASE     => 0,
                                 -END_PHASE => 0,
                                 -SLICE     => $self->slice,
                                 -ANALYSIS  => $self->create_analysis($LOGIC_NAME_EXON),
                                 -STABLE_ID => $exon_id, 
                                 -CREATED_DATE => $time,
                                 -MODIFIED_DATE => $time,
                                 -VERSION   => $ATTRIBUTE_VERSION
                                );
    print "  For Transcript ".$transcript->stable_id.", have Exon ".$ex->stable_id." with ".
          "start ".$ex->start." end ".$ex->end." strand ".$ex->strand."\n" ;
    # add Exon to Transcript
    $transcript->add_Exon($ex);
  }
  return $transcript; 
}

sub add_Translation {
  my ($self, $transcript) = @_;
  my @protein_ids;

  foreach my $child_id (keys %{$self->get_childrenID_from_parentID($transcript->stable_id)}) {
    if (${$self->get_features_from_ID($child_id)}{'type'} eq 'protein') {
      push @protein_ids, $child_id;
      print "   Found protein for transcript ".$transcript->stable_id." (strand ".$transcript->strand.") : ";
    }
  }
  if (scalar(@protein_ids) > 1) {
    print "WARNING: More than one protein for transcript ".$transcript->stable_id." @protein_ids";
  }
  if (defined $protein_ids[0]) {
    print $protein_ids[0]."\n";

    # get 3,4 and 6th attribute of gff-line
    my $cds_start = ${$self->get_features_from_ID($protein_ids[0])}{'start'};
    my $cds_end = ${$self->get_features_from_ID($protein_ids[0])}{'end'};
    my $cds_strand = ${$self->get_features_from_ID($protein_ids[0])}{'strand'};
    # add Translatio to transcript
    print "Adding tln with protein ".$protein_ids[0]." having start $cds_start end $cds_end strand $cds_strand\n";
    $transcript = add_TranslationObject($protein_ids[0], $transcript, $cds_start,$cds_end,$cds_strand);
    $transcript =  $self->setExonPhases($transcript);
    # set phases of exons
  } else {
    print "No protein for transcript ".$transcript->stable_id." and phases not set\n";
  }
  return $transcript;
}

    
sub rename_exons {
  my ($gene_id, $all_new_transcripts) = @_;
  my @transcripts = @$all_new_transcripts;
  my %all_processed_exons;

  # renaming of exons if the have the same coordiantes but different phases
  my (%hk_name,%exoncheck);
  # check for exons which have same name but different phases 
  # by getting all exon-hashkeys of one flybase-exon-identifer CG100: hk1, hk2,..
  foreach my $trans ( @transcripts ) {
    print "Finding exons for transcript $trans with stable_id ".$trans->stable_id."\n";
    foreach my $e ( @{$trans->get_all_Exons()} ) {
      if (exists $exoncheck{$e->stable_id}) {
        ${$exoncheck{$e->stable_id}}{$e->hashkey} = $e->stable_id;
        # $e->hashkey returns a unique hashkey that can be used to uniquely identify
        # this exon.  Exons are considered to be identical if they share
        # the same seq_region, start, end, strand, phase, end_phase.
        # Note that this will consider two exons on different slices
        # to be different, even if they actually are not. 
        # Returntype : string formatted as slice_name-start-end-strand-phase-end_phase
      }else{
        #warning("No stable_id for exon ".$e->hashkey."\n");
        ${$exoncheck{$e->stable_id}}{$e->hashkey} = $e->stable_id;
      }
      $hk_name{$e->hashkey}= $e->stable_id;
      # Get slice_name, start, end, strand, phase, end_phase
      my @line = split /\-/,$e->hashkey;
      ${$all_processed_exons{"$line[0]-$line[1]-$line[2]"}}{$gene_id}=();
    }
  }
  
  # find exons which have the same id (CG923:1) but different hashkeys
  my @multiple_hashkey;
  for (keys %exoncheck) {
    if (scalar (keys %{$exoncheck{$_}} )  > 1) {
      # more than one hashkey for this exons_stable_id
      for my $hashkey (keys %{$exoncheck{$_}}) {
        push @multiple_hashkey, $hashkey;
      }
    }
  }
  
  # make new names for exon with same bounds but different phases
  my %same_bounds;
  for my $mh(@multiple_hashkey) {
    my @line = split /\-/,$mh;
    ${$same_bounds{"$line[0]-$line[1]-$line[2]"}}{$mh}=();
  }
  
  # construct new name for each of the diffrent hahskeys
  for(keys %same_bounds) {
    my %exon_hk = %{$same_bounds{$_}};
    my $cnt = "A";
    for my $ehk(keys %exon_hk) {
      $hk_name{$ehk} = $hk_name{$ehk}."-$cnt";
      $cnt++;
    }
  }
  
  # rename the exons with diff. hashkeys
  for my $mh(@multiple_hashkey) {
    foreach my $trans ( @transcripts ) {
      foreach my $e ( @{$trans->get_all_Exons} ) {
        # change stable_id
        if ($e->hashkey eq $mh) {
          $e->stable_id($hk_name{$mh});
        }
      }
    }
  }
  
  # renaming exons which have different phases and which are shared by different genes
  foreach my $trans ( @transcripts ) {
    foreach my $e ( @{$trans->get_all_Exons} ) {
      # rename exons which are shared by different genes
      my @line = split /\-/,$e->hashkey; # get "root" of exon-hashkey without phase and strand-information
      # get all unique genes which share this exon
      my %genes_of_exon = %{$all_processed_exons{"$line[0]-$line[1]-$line[2]"} };
      for my $gid (keys %genes_of_exon) {
        if ($gid ne $gene_id) {
          # if the exon is shared by another gene than the actual one rename the exon
          # processed_gene : CG233
          # name of exon CG100:4
          # new name of exon: CG100:4
          my $new_name = $e->stable_id . "--" . $gene_id ;
          $e->stable_id ($new_name);
        }
      }
    }
  }
  return \@transcripts;
}  

=head2 prune_Exons

  Arg [1]   : Bio::EnsEMBL::Gene
  Function  : remove duplicate exons between two transcripts
  Returntype: Bio::EnsEMBL::Gene


=cut


sub prune_Exons {
  my ($gene) = @_;
  my @unique_Exons;

  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_Exon and translation->end_Exon

  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_Exons) {
        if ($uni->start  == $exon->start  &&
            $uni->end    == $exon->end    &&
            $uni->strand == $exon->strand &&
            $uni->phase  == $exon->phase  &&
            $uni->end_phase == $exon->end_phase
           ) {
          $found = $uni;
          last UNI;
        }
      }
      if (defined($found)) {
        push(@newexons,$found);
        if ($exon == $tran->translation->start_Exon) {
          $tran->translation->start_Exon($found);
        }
        if ($exon == $tran->translation->end_Exon) {
          $tran->translation->end_Exon($found);
        }
      } else {
        push(@newexons,$exon);
        push(@unique_Exons, $exon);
      }
    }
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  return $gene;
}


=pod


=head2 store_as_simple_feature

  Title       : store_as_simple_feature
  Usage       : $obj->store_as_simple_feature($analysis_adaptor, $type, $logic_name)
  Function    : 1st creates an analysis-object using the submitted analysis-adaptor and 
                stores all features of gff of class $type (indicated by 3rd column in gff) as SimpleFeature of a given Analysis
  Arguments   : AnalysisAdaptor-Object, $type, $logic_name
  Return-Val  : none


=cut

sub store_as_simple_feature{
  my ($self, $sf_adaptor, $type) = @_;

  # check if there is such a feature in gff
  if ($self->get_ids_by_type($type)) {
    # get all ids of simplefeaturs of a given type
    my @ids_of_given_type = @{$self->get_ids_by_type($type)}; 
    print "\n* * *\n>> Doing type $type\n* * *\n";
    print "Found ".scalar(@ids_of_given_type)." ids of type $type\n";
    # check that these ids are unique... they should be!
    my %unique_ids;
    for (@ids_of_given_type){
      if (!exists $unique_ids{$_}) {
        $unique_ids{$_} = 1;
      } else {
        warning("Non-unique entry of type $type : ID ".$_);
      }
    }
    print "Found ".scalar(keys %unique_ids)." UNIQUE ids of type $type\n";

    # ok , now we get the simple features
    my @simple_features;
    my %analyses;

    SIMPLE: foreach my $unique_id (keys %unique_ids) {                                 # process each simple_feature
      if (${$self->get_features_from_ID($unique_id)}{'type'} eq $type) {
        # make an analysis, if it does not yet exists, from the type 
        # it would be nice to use the source too, but i'm not quite sure how to bring it in to our current db schema
        my $analysis_logic = ${$self->get_features_from_ID($unique_id)}{'type'};
        if (length($analysis_logic) > 40) {
          my $old = $analysis_logic;
          $analysis_logic = substr($analysis_logic,0,40);
          warning("Shortened analysis logic_name from $old to $analysis_logic");
        }

        if (!exists $analyses{$analysis_logic}) {
          $analyses{$analysis_logic} = $self->create_analysis($analysis_logic);
        }

      # my $cds_start = ${$self->get_features_from_ID($protein_ids[0])}{'start'};
        my $score = ${$self->get_features_from_ID($unique_id)}{'score'};
        my $display = ${$self->get_features_from_ID($unique_id)}{'ID'}[0];
        if (${$self->get_features_from_ID($unique_id)}{'Name'}[0]) {
          $display = ${$self->get_features_from_ID($unique_id)}{'Name'}[0];
        }
        # we could also use the name field as a display name but for now the ID is unique and more descriptive
        my $sf = Bio::EnsEMBL::SimpleFeature->new(
                                                    -start    => ${$self->get_features_from_ID($unique_id)}{'start'}, 
                                                    -end      => ${$self->get_features_from_ID($unique_id)}{'end'}, 
                                                    -strand   => ${$self->get_features_from_ID($unique_id)}{'strand'}, 
                                                    -slice    => $self->slice,
                                                    -analysis => $analyses{$analysis_logic},
                                                    -score    => ($score eq '.' ? '0.0' : $score),
                                                    -display_label => $display,
                                                    -dbID     => undef,
                                                    -adaptor  => undef
                                                   );
        push @simple_features, $sf;
      } else {
        next SIMPLE;
      }
    } # end for
    if (scalar(@simple_features > 0)) {    
      $sf_adaptor->store(@simple_features);
      return scalar @simple_features;
    } else {
      return 0; 
    }
  } else {
    warning("Can't find any simple features of type $type");
    return 0;
  }
} # eos





=pod

=head2 create_analysis

  Title       : create_analysis
  Usage       : $obj->create_analysis ( $analysis_adaptor, $logic_name)
  Function    : 1st creates an analysis-object of the values in the FlyBaseConf-file and
                stores the analysis by using the submitted AnalysisAdaptor (normally created out of the registry-file)
  Arguments   : AnalysisAdaptor-Object, $type, $logic_name
  Return-Val  : none

=cut



sub create_analysis{
  my ($self,$logic_name, $db_file) = @_; 

  my $ana_obj =Bio::EnsEMBL::Analysis->new(
                                           -logic_name =>$logic_name,
                                           -db => $ANA_DB,
                                           -db_version => $ANA_DB_VERSION,
                                           -db_file => $ANA_DB_FILE,
                                           -program => $ANA_PROGRAM,
                                           -program_version => $ANA_PROGRAM_VERSION,
                                           -program_file => $ANA_PROGRAM_FILE,
                                           -gff_source => $ANA_GFF_SOURCE,
                                           -gff_feature => $ANA_GFF_FEATURE,
                                           -module => $ANA_MODULE,
                                           -module_version => $ANA_MODULE_VERSION,
                                           -parameters => $ANA_PARAMETERS,
                                           -created => $ANA_CREATED
                                          );
  if ($db_file) {
    $ana_obj->db_file($db_file);
  }
  return $ana_obj;
}


sub parse_gff {
  my ($self) = @_;
  my %featurehash;
  my %parent2child;
  my %m_type2id;

  open(GFF,$self->gff)|| die "Could not open ".$self->gff."\n";
  LINE: while(<GFF>) {
    # read in the gff file, line by line
    chomp;
    last if /^\#\#FASTA/;
    next if /^\#/;
    my @line = split/\t/;    # line consists of 9 tab-delimited fields in GFF Version v3.0

    # split the line up into the tab-delimited fields and store them
    my ($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attrib) = @line;
         
    # # #
    # check the strand
    # # #
    if ($strand eq '-') {
      $strand = $line[6] = "-1";
    }
    if ($strand eq '+') {
      $strand = $line[6] = "+1";
    }
    if ($strand eq '.') {
      $strand = $line[6] = "0";
    }
    if ($strand eq '?') {
      $strand = $line[6] = "0";
    }

    # # #
    # attributes
    # # #
    my $unique_id;
    my @attributes = split/;/,$attrib if $attrib;
    my $seenID;
    # possible tags are: ID, Parent, Name, Dbxref, gbunit, cyto_range etc etc
    for (@attributes) {     #process all values of an attribute Parent=ID1,ID2,ID3 and write them in an attribute_hash
      my ($attrib_name, $attrib_values ) = split /=/;
      my @values = split /,/,$attrib_values;
      # we expect the ID to come first
      if ($attrib_name eq 'ID') {
        $unique_id = $values[0];
        $seenID = 1;
        if (exists($featurehash{$unique_id})) {
          warning("Already seen ID $unique_id in line @line");
        }
      } elsif (!$seenID) {
        throw("parsing problem, need ID @line");
      }
      
      # add to our featurehash as an array
      foreach my $value (@values) {
        push @{$featurehash{$unique_id}{$attrib_name}}, $value;
        if ($attrib_name eq 'Parent' || $attrib_name eq 'Derives_from') {
          # parent --> child
          my @parents = @values;
          foreach my $parent (@parents) {
            push @{$parent2child{$parent}{$unique_id}}, ($attrib_name, $type);
            # as in, a protein derives from a transcript
          }
        }
      }
    }# for (@attributes) { 
    push @{$m_type2id{$type}}, $unique_id;
    # # # 
    # hash for this feature
    # # #
    # add to our featurehash as a scalar
    $featurehash{$unique_id}{'seqid'} = $seqid;
    $featurehash{$unique_id}{'source'} = $source;
    $featurehash{$unique_id}{'type'} = $type;
    $featurehash{$unique_id}{'start'} = $start;
    $featurehash{$unique_id}{'end'} = $end;
    $featurehash{$unique_id}{'score'} = $score;
    $featurehash{$unique_id}{'strand'} = $strand;
    $featurehash{$unique_id}{'phase'} = $phase;

    if ($featurehash{$unique_id}{'type'} eq 'orthologous_to') {
      my $new_unique_id = $unique_id."_".$featurehash{$unique_id}{'to_id'}[0];
      if (exists $featurehash{$unique_id}{'to_id'}[0]) {
        # change the featurehash
        foreach my $old_attrib_name (keys %{$featurehash{$unique_id}}) {
            $featurehash{$new_unique_id}{$old_attrib_name} = $featurehash{$unique_id}{$old_attrib_name};
        } # foreach my $old_attrib_name (keys %{$featurehash{$unique_id}}) {
        delete $featurehash{$unique_id};
        # change the typeid hash
        for (my $x=0; $x < scalar(@{$m_type2id{$type}}); $x++) {
          if (${$m_type2id{$type}}[$x] eq $unique_id) {
            ${$m_type2id{$type}}[$x] = $new_unique_id;
            next;
          }
        } 
        # change the parent2childhash
        foreach my $p (keys %parent2child) {
          if (exists $parent2child{$p}{$unique_id}) {
            @{$parent2child{$p}{$new_unique_id}} = @{$parent2child{$p}{$unique_id}};
            #push @{$parent2child{$p}{$new_unique_id}}, @{$parent2child{$p}{$unique_id}};
            delete $parent2child{$p}{$unique_id};
          }
        }
      } else {
        throw("Unable to make ID ".$featurehash{$unique_id}{'ID'}[0]." unique");
      } # if (exists $featurehash{$unique_id}{'to_id'}[0]) {
    } # if ($featurehash{$unique_id}{'type'} eq 'orthologous_to') {
  } # LINE: while(<GFF>) {

  $self->featurehash(\%featurehash);
  $self->parent2child(\%parent2child);
  $self->type2id(\%m_type2id);
  close GFF;     
}


sub parent2child {
  my ($self,$hashref) = @_;
  $self->{parent2child} = $hashref if $hashref;
  return $self->{parent2child};
}

sub featurehash {
  my ($self,$hashref) = @_;
  $self->{featurehash} = $hashref if $hashref;
  return $self->{featurehash};
}


#
#  alternative getter/setters
#
 ################################################################################



=pod

=head3 alternative getter/setters

=head2 slice()

  Usage       : $obj->slice || $obj->slice($slice)
  Function    : sets/gets the actual slice on which the features/genes are stored

=cut

sub slice{
  my ($self,$slice) = @_;
  $self->{'slice'} = $slice if $slice;
  return $self->{'slice'};
}


=pod

=head2 db()

  Usage       : $obj->db() || $obj->db($database_adaptor)
  Function    : sets/gets the DatabaseAdaptor

=cut

sub db{
  my ($self,$db) = @_;
  $self->{'db'} = $db if $db;
  return $self->{'db'};
}



=pod

=head2 gff()

  Usage       : $obj->gff() || $obj->gff( $gff_filename)
  Function    : sets/gets the filename of the current gff-file

=cut


sub gff{
  my ($self,$gff)=@_;
  if ($gff) {
    $self->{gff}=$gff;
    $self->throw("gff not found $gff: $!\n") unless (-e $gff);
  }
  return $self->{gff};
}



=pod

=head2 warn_inconsistency( $warning )

  Usage       : warn_inconsistency($warning)
  Function    : prints out warning-statement

=cut


sub warn_inconsistency{
  my $error = shift;
  print STDERR "\tWARNING: UNEXPECTED INCONSISTENCY IN GFF-FILE !!\t\t$error";
  return;
}


=pod

=head2 get_dbxrefs_($id, $xref_type )

  Usage       : $obj->get_field_attributes($id, $descriptor )
  Function    : returns the Dbxref-Identifier of a given ID like 'FBgn0033692' for type 'gene'
  Returnvalue : Arrayref. to all DBxrefs of a given type
=cut


sub get_dbxrefs {
  my ($self,$id, $xref_type) = @_;
  # eg. where $id is a Transcript_id and $xref_type is the child_type eg. mRNA, miRNA, pseudogene
  my @dbes; 

  if (exists ${$self->get_features_from_ID($id)}{'Dbxref'}) {
    DBXREF: foreach my $dbx (@{$self->get_features_from_ID($id)->{'Dbxref'}}) {
      my ($db_name, @primary_id_string) = split /\:/,$dbx;
      # split up into db_name and primary_id
      my $primary_id = join(':',@primary_id_string);
      if ($db_name eq 'if' && $primary_id !~ m/www/) {
        # Interactive Fly = IF
        $primary_id = "http://www.sdbonline.org/fly".$primary_id;
        # eg. http://www.sdbonline.org/fly/cytoskel/lethl2g1.htm
      }
      if ($db_name && $primary_id) {
        print "got xref $db_name with pid $primary_id\n";
        # make a db_entry
        my $db_entry = Bio::EnsEMBL::DBEntry->new(
                                              -DBNAME  => $db_name,
                                              -RELEASE => 1,
                                              -PRIMARY_ID => $primary_id,
                                              -DISPLAY_ID => $primary_id
                                             );

        push @dbes, $db_entry;
      } else {
        print "No xref db_name or primary_id for id $id ";
      } # if ($db_name && $primary_id) {
    } # DBXREF
  } else {
    warning("No xrefs for id $id of type $xref_type\n");
    @dbes = ();
  }
  if (exists ${$self->get_features_from_ID($id)}{'Name'}) { 
    foreach my $syn (@{${$self->get_features_from_ID($id)}{'Name'}}) {
      $dbes[0]->add_synonym($syn);
      print "Added synonym name $syn to ".$dbes[0]->primary_id."\n";
    }
  }
  if (exists ${$self->get_features_from_ID($id)}{'Alias'}) {
    foreach my $syn (@{${$self->get_features_from_ID($id)}{'Alias'}}) {
      $dbes[0]->add_synonym($syn);
      print "Added synonym alias $syn to ".$dbes[0]->primary_id."\n";
    }
  }

  # at the moment, we could parse out SO and GO too, but we don't 
  # eg. if (exists ${$self->get_features_from_ID($id)}{'Ontolgy_term'}) {
  # get the SO and Go things...
  # }
  return \@dbes;
}
sub type2id{
  my ($self,$type2id_hash) = @_;

  $self->{_type2id} = $type2id_hash if $type2id_hash;
  return $self->{_type2id};
}


=pod

=head3 getting attributes by their type ('gene')

=head2 get_ids_by_type

  Title       : get_ids_by_type
  Usage       : $obj->get_ids_by_type('gene')
  Function    : returns an Arrayref to all known IDs of a given type (in the scope of the gff)
  Arguments   : ID
  Return-Val  : Arrayref

=cut


sub get_ids_by_type{
  my ($self,$type) = @_;
  # get hash in which TYPES of the features are stored as keys return array-ref to all feature-identifiers of this TYPE
  my $idref = [];
  if (${$self->type2id}{$type}) {
    $idref = ${$self->type2id}{$type}
  }
  return $idref;
}
1;
