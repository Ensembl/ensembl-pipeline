#!/usr/local/ensembl/bin/perl  -w


package FlyBaseGff;

use strict;
use warnings;
use FlyBaseConf;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate verbose);
use Bio::EnsEMBL::SimpleFeature;


# construct newFlyBaseGff-object (GFF Version 3)

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


  $self->parse_gff();

  return $self;
}



# access all stored attributes of a featureid


sub attributes_by_id{
  my ($self,$hashref) = @_;

  $self->{attrib_by_id} = $hashref if $hashref;
  return $self->{attrib_by_id};
}





=pod

=head2 store_genes

  Title       : store_genes
  Usage       :
  Function    :
  Arguments   :
  Return-Val  :

=cut


sub store_genes{
  my ($self) = @_;

  my @genes = 

  #
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
  my $analysis = $self->get_analysis($logic_name);




  # check if there is such a feature in gff
  if ($self->get_ids_by_type($type)){

    my @ids_of_given_type = @{$self->get_ids_by_type($type)};   # get all ids of simplefeaturs of a given type

    my @simple_features;

    for(@ids_of_given_type){                                 # process each simple_feature
      my @attributes = @{${$self->attributes_by_id}{$_}};    # get attributes of SimpleFeature
      warn_inconsistency("Featureid $_ of TYPE <$type> is not unique\n") if (scalar(@attributes) > 1);

      # creation and storing of simpe-feature
      my ($seqid,$source,$type,$start,$end,$score,$strand,$phase,$href) = @{$attributes[0]};

      # process attributes of the SimpleFeature to get the DisplayId if no $label is given
      my $display_label;
      unless ($label){   # use id as label if no label has been supplied
        my %vals = %$href;
        $display_label = ${$vals{ID}}[0];
        print $display_label."\n";
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

=head2 non_unique_ids

  Title       : non_unique_ids
  Usage       : $obj->non_unique_ids() or $obj->non_unique_ids('ID89234')
  Function    : stores all non-unique ids in a hash and returns this hahs if called without argument
  Arguments   : FeatureIdentifier, none
  Return-Val  : Ref. to hash with all non-unique feature-ids as keys (and null as  value)

=cut



sub non_unique_ids{
  my ($self,$id) = @_;
  ${$self->{non_unique_id}}{$id}=() if $id;
  return $self->{non_unique_id};
}




=pod

=head2 get_analysis

  Title       : get_analysis
  Usage       : $obj->get_analysis ( $analysis_adaptor, $logic_name)
  Function    : 1st creates an analysis-object of the values in the FlyBaseConf-file and
                stores the analysis by using the submitted AnalysisAdaptor (normally created out of the registry-file)
  Arguments   : AnalysisAdaptor-Object, $type, $logic_name
  Return-Val  : none

=cut


sub get_analysis{
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




#=pod

#=head2 genes

#  Title       : genes
#  Usage       : $obj->genes(), $obj->genes('ID')
#  Function    : gets / stores all IDs which have type '^gene%'
#  Arguments   : optional ID
#  Return-Val  : Reference to an array containing all gene-IDs

#=cut


#sub genes{
#  my ($self, $gene_id) = @_;

#  push @{$self->{_gene}}, $gene_id if $gene_id;
#  return \@{$self->{_gene}};
#}




=pod

=head2 get_ids_by_type

  Title       : get_ids_by_type
  Usage       : $obj->get_ids_by_type('ID')
  Function    : returns an Arrayref to all known types of a given ID (in the scope of the gff)
  Arguments   : ID
  Return-Val  : Arrayref

=cut


sub get_ids_by_type{
  my ($self,$id) = @_;
  # get hash in which the types of the features are stored as keys
  # and return array-ref to all feature-identifiers of this TYPE
  return ${$self->type2id}{$id};
}




=pod

=head2 type2id

  Title       : type2id
  Usage       : $obj->type2id("TYPE") i.e. $obj->type2id("gene")
  Function    : gets/sets a ref to a hash which contains all known feature-types as keys and their representatives as values
  Arguments   : HASH, none
  Return-Val  : Hashref.

=cut


sub type2id{
  my ($self,$type2id_hash) = @_;

  $self->{_type2id} = $type2id_hash if $type2id_hash;
  return $self->{_type2id};
}




=pod

=head2 diff_types_of_given_id

  Title       : diff_types_of_given_id
  Usage       : $obj->diff_types_of_given_id('LP02895:contig1');
  Function    : returns arrayref containing all types of a given id (because of inconsistency!)
  Arguments   : unique_identifier like 'LP02895:contig1'
  Return-Val  : Hashref of all identifiers and their associated types ('LP02895:contig1' => qw(mRNA gene EST))

=cut


sub diff_types_of_given_id{
  my $self = shift;
  $self->{_diff_types} = shift if @_ ;
  return $self->{_diff_types};
}



=pod

=head2 get_childs_of_id

  Title       : get_childs_of_id
  Usage       : $obj->child('LP02895:contig1')
  Function    : returns child(s)-ID of a given parent-ID
  Arguments   : String describing id
  Return-Val  : arrayref

=cut


sub get_childs_of_id{
  my ($self,$id) = @_;

  my %tmp = %{$self->parent2childs};
  return $tmp{$id};
}




=pod

=head2 parent2childs

  Title       : parent2childs
  Usage       : $obj->parent2childs(%HASH || none)
  Function    :  parent/child-relation : contains id of all items which are refenced as "PARENT" in the gff and their associated childs
  Arguments   : %HASH, none
  Return-Val  : hashref

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



sub warn_inconsistency{
  my $error = shift;
  print STDERR "\tWARNING: UNEXPECTED INCONSISTENCY !!\t\t$error";
  return;
}





# getter / setter of slice
sub slice{
  my ($self,$slice) = @_;
  $self->{'slice'} = $slice if $slice;
  return $self->{'slice'};
}



# getter/setter for db
sub db{
  my ($self,$db) = @_;
  $self->{'db'} = $db if $db;
  return $self->{'db'};
}

# getter/setter of gff-filename
sub gff{
  my ($self,$gff)=@_;
  if($gff){
    $self->{gff}=$gff;
    $self->throw("gff not found $gff: $!\n") unless (-e $gff);
  }
  return $self->{gff};
}






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
  my %diff_types_of_given_id;             # stores the types associated with an ID-entry   key:ID ==> value = types (gene,mRNA,EST,..)

  # global (gff)
  my %m_parent2childs;                          # hash containing the parent's id and the id of their childs
  my %m_type2id;                                   # stores all feature_identifiers of the same type (mRNA,gene) in Hash
  my %non_unique_feature_identifiers;     # collection of non-unique identifiers

  my %id2attributes;                             # collection of all attributes of a given feature_id;

################################################################################
my $cnt=0;
  open(GFF,$self->gff)|| die "Could not open ".$self->gff."\n";
  while(<GFF>){
    chomp;
    next if /\#/;
    my @line = split/\t/;    # line consists of 9 tab-delimited fields in GFF Version v3.0
    my ($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attrib) = @line;


    # ENSEMBL-STRAND-CONVERSION
    if($strand eq '-'){
      $strand = $line[6] = "-1";
    }


    if($strand eq '+'){
      $strand = $line[6] = "+1";
    }


    if($strand eq '.'){
      $strand = $line[6] = "0";
    }


    # split semicolon-separated attributes-field
    my @attributes = split/;/,$attrib if $attrib;

    # convert list of attributes to hash, whereas the key is the tag-pair i.e. tag=value ID=200
    # possible tags are: ID, Parent, Name, Dbxref, gbunit, cyto_range

    my (%tags) ;

    for(@attributes){     #process all values of an attrribute Parent=ID1,ID2,ID3 and write them in an attribute_hash
      my ($id_name, $id_values ) = split /=/;
      # $id_name = Parent or ID....

      my @values = split /,/,$id_values;
      push @{$tags{$id_name}}, @values;
    }


    # GET (SEMI-UNIQUE) ID OF LINE/OBJECT
    my $semi_unique_id;
    if (exists $tags{'ID'}){
      $semi_unique_id =${$tags{'ID'}}[0];
      # STORE ALL FEATURE_IDENTIFIERS OF THE SAME TYPE IN A HASH
      push @{$m_type2id{$type}} , $semi_unique_id ;

      $non_unique_feature_identifiers{$semi_unique_id}++;
  }


    # IF THE ITEM PROCESSED HAS PARENTS it has also an ID
    if (exists $tags{'Parent'}){
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
     push @{$diff_types_of_given_id{$semi_unique_id}}, $type;


    # store all of line in array for whole gff
    $self->gff_data( [@line[0..7], \%tags ]);

    # store all attributes of a feature by the feature-id
    push @{$id2attributes{$semi_unique_id}}, [@line[0..7], \%tags ];

    #    # store all IDs of gff-entryis of type "gene" in an array
    #    $self->genes($semi_unique_id) if ($type =~ m/^gene$/);

  }
  close(GFF);


  #
  #  store items with non-unique identifers
  #
   ################################################################################


  for(keys %non_unique_feature_identifiers){
    if ($non_unique_feature_identifiers{$_}>1){
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
  $self->diff_types_of_given_id(\%diff_types_of_given_id);

  # make all attributes of a featureid accessible
  $self->attributes_by_id(\%id2attributes);


} # end sub






1;
