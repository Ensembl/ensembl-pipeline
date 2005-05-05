#!/usr/local/ensembl/bin/perl  -w

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

$|=1;

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
  $self->{non_unique_id}={};
  $self->{exons}={};
  $self->{db_entry}={};
  $self->parse_gff();
  print "gff $gff parsed\n";
  return $self;
}

sub db_entry{
  my ($self, $tr,$db_entry) = @_;
  if($db_entry){
    $ {$self->{db_entry}} {$tr}=$db_entry;
  }
  return $ { $self->{db_entry}}{$tr};
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
  return ${$self->type2id}{$type};
}


=pod

=head2 get_all_attributes_by_id

  Title       : get_all_attributes_by_id
  Usage       : $obj->get_all_attributes_by_id('ID')
  Function    :  access all stored attributes of a feature-id
  Arguments   : ID
  Return-Val  : arrayref to arrayrefs of attributes of a given ID

=cut


sub get_all_attributes_by_id{
   my ($self,$id) = @_;
   return ${$self->id2attributes}{$id};
}




=pod

=head2 get_source_by_id

  Title       : get_source_by_id
  Usage       : $obj->get_source_by_id('ID')
  Function    : access the second filed of the gff (source) by feature-id
  Arguments   : ID
  Return-Val  : String describing source i.e. "fgenesh" "part_of" "piecegenie" "predicted" ...

=cut


sub get_source_by_id{
   my ($self,$id) = @_;
   my @atr = @ { $ {$self->id2attributes } {$id} } ;
   if(scalar(@atr)>1){
     warn("Id $id appears in more than one gff-line and seems to be not unique\n");
   }
   return $ { $atr[0]} [1];
}



=pod

=head2 get_all_childs_of_id

  Title       : get_all_childs_of_id
  Usage       : $obj->get_all_childs_of_id'ID')
  Function    : returns child(s)-ID of a given parent-ID
  Arguments   : String describing id
  Return-Val  : arrayref

=cut


sub get_all_childs_of_id{
  my ($self,$id) = @_;
  return ${$self->parent2childs}{$id};
}


=pod

=head2 get_spec_childs_of_id

  Title       : get_spec_childs_of_id
  Usage       : $obj->get_spec_childs_of_id('ID','type')
  Function    : returns child(s)-IDs of a given parent-ID which have the specified type
  Arguments   : String describing ID and type
  Return-Val  : arrayref

=cut

 # work
sub get_spec_child_of_id{
  my ($self,$parent_id,$spec_type) = @_;
  my @matching_ids;

#    print "get_spec_child_of_id: Trying to get child of parent -$parent_id- (-$spec_type-)\n" ; 

  if ( ${ $self->parent2childs} {$parent_id} ) {
    my @poss_child_ids = @{${$self->parent2childs}{$parent_id}};
    my $parent_has_no_child_of_type = 1 ;
    for my $child (@poss_child_ids) {
      for my $typ(@{ $self->get_type_of_id($child) } ) {
        if ($typ =~ m/$spec_type/) {
          push @matching_ids,$child;
	  #print "CHILD $child\n" ;
	  $parent_has_no_child_of_type = 0 ;
        }
      }
    }
    #warn("no child for $parent_id of type $spec_type") if($parent_has_no_child_of_type);
    return \@matching_ids;
  }else{
      warn("no child for parent $parent_id of type: $spec_type  found\n");
      return \@matching_ids;
  }
}

=pod

=head2 get_type_of_id

  Title       : get_type_of_id
  Usage       : $obj->get_type_of_id('ID')
  Function    : returns type(s) of a given Id
  Arguments   : String describing ID
  Return-Val  : arrayref

=cut


sub get_type_of_id{
  my ($self,$id) = @_;
  return ${$self->id2types}{$id};
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
  Arguments   :
  Return-Val  :

=cut


sub store_as_gene_object{
  my ($self,$target, $child_type, $ensembl_gene_type) = @_;
  print "Target: =$target= -> =$child_type=\n";


  print "Having gene_ids:\n" ; 
  print join(" will be processed.\n", @{$self->get_ids_by_type($target)} ) . "\n" ;

  my %all_processed_exons;

  # test if $target is annotated in gff
  unless ($self->get_ids_by_type($target)){
    print "no $target in ".$self->gff()." found\n";
    return;
  }




  for my $gene_id ( @{$self->get_ids_by_type($target)} ) {

    if($self->get_source_by_id($gene_id) ne "."){
      print  "WARN: $gene_id will be skipped because it has unrecognized source:" . $self->get_source_by_id($gene_id) . "\n";
      next;
    }


   my @all_nw_transcripts;

    #---- start processing the mRNA-childs (transcripts) ----


#    print $self->get_spec_child_of_id($gene_id,"$child_type") . "\n" ;

    if (scalar(@{$self->get_spec_child_of_id($gene_id,"$child_type")}) ==0){
      #warn ("no child of category $child_type of gene $gene_id found");
      next;
    }else{

      for my $trans_id (@{$self->get_spec_child_of_id($gene_id,"$child_type")}) {

        my $nw_transcript = new Bio::EnsEMBL::Transcript(
                                                       -STABLE_ID => $trans_id,
                                                       -VERSION => 3,
                                                      );

       my $db_entry = $self->get_dbxref($trans_id, $child_type );

        if($db_entry){
          $self->db_entry($nw_transcript,$db_entry);
        }


        #---- start store each exon ----
        my @exons =  @{$self->get_spec_child_of_id($trans_id,"exon")};

        for my $exon (@exons) {
          my @gff_line = @{$self->get_all_attributes_by_id($exon)};
          my ($seqid,$source,$type,$ex_start,$ex_end,$score,$ex_strand,$phase,$attrib) = @{$gff_line[0]};
          warn_inconsistency("ID of exon $exon is not unique \n") if (exists $gff_line[1]);

          my $ex = Bio::EnsEMBL::Exon->new(
                                           -START     => $ex_start,
                                           -END       => $ex_end,
                                           -STRAND    => $ex_strand,
                                           -PHASE     => 0,
                                           -END_PHASE => 0,
                                           -SLICE     => $self->slice,
                                           -ANALYSIS  => $self->create_analysis($LOGIC_NAME_EXON),
                                           -STABLE_ID => ${  ${ $attrib }{'ID'} }[0] ,
                                           -VERSION   => 3
                                          );
	  #print "having exon $ex_start\t$ex_end\n" ; 

          my $esid = ${  ${ $attrib }{'ID'} }[0] ;
          $nw_transcript->add_Exon($ex);
        }


        #---- add Translation to transcript if CDS

       #  check if there is a CDS annotated for the mRNA and get coordiantes
        my $cds_id = ${$self->get_spec_child_of_id($trans_id,"CDS")}[0];

        if ($cds_id) {
          warn_inconsistency("There is a mRNA with more than one CDS\n")
           if exists ${$self->get_spec_child_of_id($trans_id,"CDS")}[1];

          # get 3,4 and 6th attribute of gff-line
          my ($cds_start,$cds_end,$cds_strand) =@{${$self->get_all_attributes_by_id($cds_id)}[0]}[3,4,6];

          # add Translatio to transcript
          $nw_transcript = add_TranslationObject($cds_id,$nw_transcript,$cds_start,$cds_end,$cds_strand);

          # set phases of other exons
          $nw_transcript =  $self->setExonPhases($nw_transcript);

      }else{
         warn("mRNA $trans_id has no CDS\n") if ($target eq "gene" && $child_type eq 'mRNA' );
      }

#       add_DBEntry($nw_transcript);

       push @all_nw_transcripts, $nw_transcript;



    } #--- end each alternative transcript ---



    # renaming of exons if the have the same coordiantes but different phases
    ################################################################################
    my (%hk_name,%exoncheck);

    # check for exons which have same name but different phases 
    # by getting all exon-hashkeys of one flybase-exon-identifer CG100: hk1, hk2,..

    foreach my $trans ( @all_nw_transcripts ) {
      foreach my $e ( @{$trans->get_all_Exons} ) {
        if (exists $exoncheck{$e->stable_id}) {
          ${$exoncheck{$e->stable_id}}{$e->hashkey} = $e->stable_id;;
        }else{
          ${ $exoncheck{$e->stable_id} }  {$e->hashkey} = $e->stable_id;
        }
        $hk_name{$e->hashkey}= $e->stable_id;
        my @line = split /\-/,$e->hashkey;
       ${$all_processed_exons{"$line[0]-$line[1]-$line[2]"}}{$gene_id}=();
     }
   }
   # find exons which have the same id (CG923:1) but different hashkeys
   my @multiple_hashkey;
   for (keys %exoncheck) {
     if (scalar (keys %{$exoncheck{$_}} )  > 1) {
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
     foreach my $trans ( @all_nw_transcripts ) {
       foreach my $e ( @{$trans->get_all_Exons} ) {
         # change stable_id
         if ($e->hashkey eq $mh) {
           $e->stable_id($hk_name{$mh});

         }
       }
     }
   }


# renaming exons which have different phases and which are shared by different genes

foreach my $trans ( @all_nw_transcripts ) {
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





################################################################################









  # got all alternative transcripts and make genes
  my ($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attrib) = @{$self->get_all_attributes_by_id($gene_id)};

  #  print "C-processing gene with id : $gene_id    " ;
 	
  my $gene = Bio::EnsEMBL::Gene->new(
                                     -START     => $start,
                                     -END       => $end,
                                     -STRAND    => $strand,
                                     -SLICE     => $self->slice,
                                     -ANALYSIS  => $self->create_analysis($LOGIC_NAME_EXON),
                                     -STABLE_ID => $gene_id,
                                     -VERSION => 3,
	                             -TYPE    =>$ensembl_gene_type,
                                    );

   map ($gene->add_Transcript($_) , @all_nw_transcripts) ;

   # check if we already have an exon-stable-id/gene-stable-id/translation-stable-id/transcript-stable-id for this gene


   print "\nNow trying to store gene $gene_id....." ;
   my $gene_DB_id = $self->db->get_GeneAdaptor->store($gene);
   print "STORED\n\n\n" ;
   my $db_entry = $self->get_dbxref($gene_id, $target);
   $self->db->get_DBEntryAdaptor->store($db_entry,$gene_DB_id,'Gene');

   # reprocess the db_entry for trs

   for my $tr (@all_nw_transcripts){
     my $dbe = $self->db_entry($tr);
     if($dbe){
       $self->db->get_DBEntryAdaptor->store($dbe , $tr->dbID() ,'Transcript');
     }
   }


  } #end foreach gene
 } #end unless...else
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
  my ($self, $sf_adaptor, $type, $logic_name, $label) = @_;

  # add analysis to db (the storage of the analysis is handeled by SimpleFeatureAdaptor
  my $analysis = $self->create_analysis($logic_name);

  # check if there is such a feature in gff
  if ($self->get_ids_by_type($type)) {

    my @ids_of_given_type = @{$self->get_ids_by_type($type)};   # get all ids of simplefeaturs of a given type

    my @simple_features;

    for(@ids_of_given_type) {                                 # process each simple_feature
      my @attributes = @{${$self->id2attributes}{$_}};    # get attributes of SimpleFeature
      warn_inconsistency("Featureid $_ of TYPE <$type> is not unique\n") if (scalar(@attributes) > 1);

      # creation and storing of simpe-feature
      my ($seqid,$source,$type,$start,$end,$score,$strand,$phase,$href) = @{$attributes[0]};

      # process attributes of the SimpleFeature to get the DisplayId if no $label is given
      my $display_label;
      unless ($label) {   # use id as label if no label has been supplied
        my %vals = %$href;
        $display_label = ${$vals{ID}}[0];
      }else{
       $display_label = $label;
      }

      my $sf = Bio::EnsEMBL::SimpleFeature->new(
                                                -start    => $start,
                                                -end      => $end,
                                                -strand   => $strand,
                                                -slice    => $self->slice,
                                                -analysis => $analysis,
                                                -score    => $score,
                                                -display_label => $display_label,
                                                -dbID     => undef,
                                                -adaptor  => undef
                                               );
     push @simple_features, $sf;

   } # end for
   $sf_adaptor->store(@simple_features);
   return scalar @simple_features;
  }else{
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
  my ($self,$logic_name) = @_; 

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
  return $ana_obj;
}





# 
# 2R	.	tRNA	1281033	1281105	.	+	.	ID=CR30304;Name=tRNA:R2:42Ad;Dbxref=FlyBase:FBan0030304,FlyBase:FBgn0003757,FlyBase:FBgn0011958,FlyBase:FBgn0050304;cyto_range=42A12-42A12;gbunit=AE003784
# 2R	.	mRNA	1281033	1281105	.	+	.	ID=CR30304-RA;Name=tRNA:R2:42Ad-RA;Dbxref=FlyBase:FBtr0085979,FlyBase:CR30304-RA,FlyBase:FBgn0003757;Parent=CR30304
# 2R	trnascan	tRNA	1281033	1281105	.	.	.	ID=NULL:1954652
# 2R	.	exon	1281240	1281312	.	-	.	ID=NULL:1958261;Parent=NULL:1958260,CR32837-RA
# 2R	.	tRNA	1281240	1281312	.	-	.	ID=CR32837;Name=tRNA:K2:42Ae;Dbxref=FlyBase:FBan0032837,FlyBase:FBgn0003795,FlyBase:FBgn0011891,FlyBase:FBgn0052837;cyto_range=42A12-42A12;gbunit=AE003784
# 2R	.	mRNA	1281240	1281312	.	-	.	ID=CR32837-RA;Name=tRNA:K2:42Ae-RA;Dbxref=FlyBase:FBtr0086017,FlyBase:CR32837-RA,FlyBase:FBgn0003795;Parent=CR32837
# 2R	trnascan	tRNA	1281240	1281312	.	.	.	ID=NULL:1958260
# 2R	.	exon	1284763	1284885	.	+	.	ID=NULL:1954611;Parent=NULL:1954606
# 2R	genscan	mRNA	20276539	20291074	.	.	.	ID=NULL:310788
# 2R	.	exon	20276539	20276718	.	+	.	ID=NULL:310789;Parent=NULL:310788
# 2R	.	gene	20290099	20291074	.	+	.	ID=CG30429;Dbxref=FlyBase:FBan0030429,FlyBase:FBgn0050429;cyto_range=60F5-60F5;gbunit=AE003466
# 2R	.	mRNA	20290099	20291074	.	+	.	ID=CG30429-RA;Dbxref=FlyBase:FBtr0072450,FlyBase:CG30429-RA,FlyBase:FBgn0050429;Parent=CG30429
# 2R	.	exon	20290099	20290576	.	+	.	ID=CG30429:1;Parent=CG30429-RA
# 2R	.	exon	20290225	20290423	.	+	.	ID=NULL:310790;Parent=NULL:310788
# 2R	.	exon	20290635	20291074	.	+	.	ID=CG30429:2;Parent=CG30429-RA,NULL:310788
# 2R	genscan	mRNA	20295508	20298047	.	.	.	ID=NULL:310736
# 2R	.	exon	20295508	20295539	.	+	.	ID=NULL:310737;Parent=NULL:310736
# 2R	.	gene	20296884	20298047	.	+	.	ID=CG30428;Dbxref=FlyBase:FBan0030428,FlyBase:FBgn0050428,FlyBase:FBgn0064772;cyto_range=60F5-60F5;gbunit=AE003466
# 2R	.	mRNA	20296884	20298047	.	+	.	ID=CG30428-RA;Dbxref=FlyBase:FBtr0072451,FlyBase:CG30428-RA,FlyBase:FBgn0050428;Parent=CG30428
# 2R	.	exon	20296884	20297010	.	+	.	ID=CG30428:1;Parent=CG30428-RA
# 2R	.	exon	20297151	20297331	.	+	.	ID=CG30428:2;Parent=CG30428-RA,NULL:310736
# 2R	.	exon	20297390	20297648	.	+	.	ID=CG30428:3;Parent=CG30428-RA,NULL:310736
# 2R	.	exon	20297698	20298047	.	+	.	ID=CG30428:4;Parent=CG30428-RA,NULL:310736










#
# GFF -  P A R S E R
#
 ################################################################################


sub parse_gff{
  my ($self) = @_;

    #
    # inconsinstency !! ID=Q9W0V7 appears twotimes in gff dmel_3L_r3.2.1.gff, also id for repeat_region ID=yoyo
    #


  # local (line)
  my %id2types;             # stores the types associated with an ID-entry   key:ID ==> value = types (gene,mRNA,EST,..)

  # global (gff)
  my %m_parent2childs;                          # hash containing the parent's id and the id of their childs
  my %m_type2id;                                   # stores all feature_identifiers of the same type (mRNA,gene) in Hash
  my %non_unique_feature_identifiers;     # collection of non-unique identifiers

  my %id2attributes;                             # collection of all attributes of a given feature_id;

################################################################################
my $cnt=0;
  open(GFF,$self->gff)|| die "Could not open ".$self->gff."\n";
  while(<GFF>) {
    chomp;
    next if /^\#/;
    my @line = split/\t/;    # line consists of 9 tab-delimited fields in GFF Version v3.0
    my ($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attrib) = @line;

    # ENSEMBL-STRAND-CONVERSION
    if ($strand eq '-') {
      $strand = $line[6] = "-1";
    }


    if ($strand eq '+') {
      $strand = $line[6] = "+1";
    }


    if ($strand eq '.') {
      $strand = $line[6] = "0";
    }


    # split semicolon-separated attributes-field
    my @attributes = split/;/,$attrib if $attrib;

    # convert list of attributes to hash, whereas the key is the tag-pair i.e. tag=value ID=200
    # possible tags are: ID, Parent, Name, Dbxref, gbunit, cyto_range

    my (%tags) ;

    for(@attributes) {     #process all values of an attrribute Parent=ID1,ID2,ID3 and write them in an attribute_hash
      my ($id_name, $id_values ) = split /=/;
      # $id_name = Parent or ID....

      my @values = split /,/,$id_values;
      push @{$tags{$id_name}}, @values;
    }


    # GET (SEMI-UNIQUE) ID OF LINE/OBJECT
    my $semi_unique_id;
    if (exists $tags{'ID'}) {
      $semi_unique_id =${$tags{'ID'}}[0];
      # STORE ALL FEATURE_IDENTIFIERS OF THE SAME TYPE IN A HASH
      push @{$m_type2id{$type}} , $semi_unique_id ;

      $non_unique_feature_identifiers{$semi_unique_id}++;
  }


    # IF THE ITEM PROCESSED HAS PARENTS it has also an ID
    if (exists $tags{'Parent'}) {
      my $child_id = $semi_unique_id;
      my @parents = @{$tags{'Parent'}};

      foreach my $parent_id (@parents) {
        warn_inconsistency("Child has no ID but parent: $parent_id\n")  unless ($child_id);
        push @{$m_parent2childs{$parent_id}} , $child_id;
      }
    }


    #
    # line-depentend tasks
    #
     ######################################################################


    # all types of an ID are stored in a hash (key: ID, values = types)
    # to push all in a temp hash and store the temp-hash is faster than doing it in an object-method
     push @{$id2types{$semi_unique_id}}, $type;


    # store all of line in array for whole gff
    $self->gff_data( [@line[0..7], \%tags ]);

    # store all attributes of a feature by the feature-id
    push @{$id2attributes{$semi_unique_id}}, [@line[0..7], \%tags ];

#    print "-->SEMI-unique_id $semi_unique_id: " ;
#    print join("\t",@line[0..7]) ."\n" ; 
    #    # store all IDs of gff-entryis of type "gene" in an array
    #    $self->genes($semi_unique_id) if ($type =~ m/^gene$/);

  }
  close(GFF);


  #
  #  store items with non-unique identifers
  #
   ################################################################################


  for(keys %non_unique_feature_identifiers) {
    if ($non_unique_feature_identifiers{$_}>1) {
      $self->non_unique_ids($_);
    }
  }



  # STORE ALL ITEMS
  ###############################################################################

  # store features which are acting as parents and their childs
  $self->parent2childs(\%m_parent2childs);


  # store a refrence to all features of the same TYPE (features are identified by their unique_id_nr
  $self->type2id(\%m_type2id);


  # store all types of an ID (key=ID, value = types (mRNA,gene)
  $self->id2types(\%id2types);

  # make all attributes of a featureid accessible
  $self->id2attributes(\%id2attributes);

} # end sub




#
# Methods used for storing parsed data in Hashes or Arrays
#
 ################################################################################



=pod

=head3 Methods used for storing parsed data in Hashes or Arrays

=head2 non_unique_ids

  Title       : non_unique_ids
  Usage       : $obj->non_unique_ids() or $obj->non_unique_ids('ID89234')
  Function    : stores all non-unique ids in a hash and returns this hash if called without argument
  Arguments   : FeatureIdentifier, none
  Return-Val  : Ref. to hash with all non-unique feature-ids as keys (and null as  value)

=cut

sub non_unique_ids{
  my ($self,$id) = @_;
  ${$self->{non_unique_id}}{$id}=() if $id;
  return $self->{non_unique_id};
}






=pod

=head2 id2types

  Title       : id2types
  Usage       : $obj->id2types( \%HASH) || $obj->id2types()
  Function    : returns all ids and their associated types (because of inconsistency!)
  Arguments   : %HASHREF, none
  Return-Val  : Hashref

=cut


sub id2types{
  my $self = shift;
  $self->{_diff_types} = shift if @_ ;
  return $self->{_diff_types};
}




=pod

=head2 id2attributes

  Title       : id2attributes
  Usage       : $obj->id2attributes() || $obj->id2attributes(\%HASH)
  Function    : returns ref to hash containing feature-ids as keys and associated attributes as values
  Arguments   : none | hashref
  Return-Val  : Hashref

=cut


sub id2attributes{
  my ($self,$hashref) = @_;

  $self->{id2attributes} = $hashref if $hashref;
  return $self->{id2attributes};
}



=pod

=head2 type2id

  Title       : type2id (used by get_ids_by_type)
  Usage       : $obj->type2id("TYPE") i.e. $obj->type2id("gene")
  Function    : gets/sets a ref to a hash which contains all known feature-types as keys and their representatives (Ids) as values
  Arguments   : HASH, none
  Return-Val  : Hashref.

=cut

sub type2id{
  my ($self,$type2id_hash) = @_;

  $self->{_type2id} = $type2id_hash if $type2id_hash;
  return $self->{_type2id};
}



=pod

=head2 parent2childs (used by get_all_childs_of_id)

  Title       : parent2childs
  Usage       : $obj->parent2childs( \%HASH ) || $obj->parent2childs()
  Function    :  parent/child-relation  contains id of all items which are refenced as "PARENT" in the gff and their associated childs
  Arguments   : %HASH, none
  Return-Val  : Hashref

=cut

sub parent2childs{
  my $self = shift;

  $self->{_parent2childs} = shift if @_ ;
  return $self->{_parent2childs};
}



=pod

=head2 gff_data

  Title       : gff_data
  Usage       : $obj->gff_data(\@arrary), $obj->gff_data();
  Function    : stores all lines of gff-data in Array and returns them if no argument is given
  Arguments   : \@arrayref, none
  Return-Val  : Reference to an Array containig the data of the gff-lines  [8 SCALAR-VALUES and one HASHREF] in each line
                [3L,.,gene, 561966, 563439, ., + ,. ref(%HASH) ]
                the returnend HASH contains the field-descriptors of the attribute-field as keys, e.g "ID", "PARENT", ...


=cut


sub gff_data{
  my ($self,$line) = @_;
  push @{$self->{_gff_data}}, $line if $line;
  return \@{$self->{_gff_data}};
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

=head2 get_dbxref_($id, $xref_type )

  Usage       : $obj->get_field_attributes($id, $descriptor )
  Function    : returns the Dbxref-Identifier of a given ID like 'FBgn0033692' for type 'gene'
  Returnvalue : Arrayref. to all DBxrefs of a given type
=cut



sub get_dbxref{
  my ($self,$id, $xref_type) = @_;
  my ($primary_id,$db_name);


  my %fb_id = (
               'CDS' =>  'FBpp',
               'exon' => 'FBex',
               'gene' => 'FBgn',
               'mRNA' => 'FBtr',
               'ncRNA' => 'FBgn',
               'pseudogene' => 'FBgn',
               'rRNA' => 'FBgn',
               'tRNA' => 'FBgn',
               'transposable_element' => 'FBti',
               'snRNA' => 'FBgn',
               'snoRNA' => 'FBgn',
              );

  my @found = grep {/$fb_id{$xref_type}/} @ { $ { pop @{ $ { $self->get_all_attributes_by_id($id) } [0] } } {'Dbxref'} } ;

  map {s/FlyBase://g} @found;

  $primary_id = shift @found;

  if($xref_type eq "mRNA"){
    $db_name = "flybase_transcript_id";
  }else{
    $db_name = "flybase_gene_id";
  }

  if($primary_id){
    my $db_entry = Bio::EnsEMBL::DBEntry->new(
                                              -DBNAME  => $db_name,
                                              -RELEASE => 1,
                                              -PRIMARY_ID => $primary_id,
                                              -DISPLAY_ID => $primary_id
                                           );
    return $db_entry;
  }
  return 0;
}





1;
