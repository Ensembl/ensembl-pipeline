=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Finished_HalfwiseHMM

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Finished_HalfwiseHMM->new(-db       => $db,
					                                    -input_id => $id
                                                                           );
    $obj->fetch_input;
    $obj->run;

    my @newfeatures = $obj->output;

=head1 DESCRIPTION

runs HalfwiseHMM runnable and converts it output into genes which can be stored in an ensembl database

=head1 CONTACT

lec@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Finished_HalfwiseHMM;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Root;

use strict;

use Bio::EnsEMBL::Pipeline::Runnable::Finished_HalfwiseHMM;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Similarity qw(GB_SIMILARITY_DATABASES);
use BlastableVersion;

use Data::Dumper;

our @ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head1  new

    Arg      : all those inherited from RunnableDB
    Function : Make a new HalfwiseHMM object defining the above variables
    Exception: thrown if no input id is provided, this is found in RunnableDB
    Caller   : 
    Example  : $runnable = Bio::EnsEMBL::Pipeline::RunnableDB::HalfwiseHMM new->(-db => $db
										 -INPUT_ID => $id
										 -ANALYSIS => $analysis);

=cut


sub new {
    my ($new,@args) = @_;
    my $self = $new->SUPER::new(@args);    
           
    # db, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)
    my ($type, $threshold) = $self->_rearrange([qw(TYPE THRESHOLD)], @args);
    $self->{'_fplist'} = []; #create key to an array of feature pairs

    # SET location for Pfam_ls unless set from db analysis.db_file
    $self->analysis->db_file('/data/blastdb/Ensembl/Pfam_ls') 
        unless $self->analysis->db_file();

    return $self;
}


sub type {
    my ($self,$type) = @_;
    
    if (defined($type)) {
        $self->{_type} = $type;
    }
    return $self->{_type};
}

sub threshold {
    my ($self,$threshold) = @_;

    if (defined($threshold)) {
        $self->{_threshold} = $threshold;
    }
    return $self->{_threshold};
}


=head1  fetch_input

    Arg      : none
    Function : fetches the repeatmasked sequence and the swall features for the contig being run on and creates the HalfwiseHMM Runnable
    Exception: throws if no input_id has been provided
    Caller   : 
    Example  : 

=cut

sub fetch_input {
    my( $self) = @_;
    print STDERR "running fetch input\n"; 
    my $i_care = 1;
    $self->check_db_versions() if $i_care;

    my %ests;
    my @estseqs;
    $self->throw("No input id") unless defined($self->input_id);
      
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    print "got contig\n";
    my $genseq   = $contig->get_repeatmasked_seq;
    print "got dnaseq\n";
  
    my $alignAdaptor = $self->db->get_ProteinAlignFeatureAdaptor();
  
    foreach my $database(@{$GB_SIMILARITY_DATABASES}){
        my $fps = [];
        my $features = $alignAdaptor->fetch_all_by_RawContig_and_pid($contig, 
                                                                     $database->{'threshold'}, 
                                                                     $database->{'type'});
        
        print STDERR "Number of features matching threshold $database->{'threshold'} = " . scalar(@$features) . "\n";

        foreach my $f (@$features) {
            if ($f && $f->isa("Bio::EnsEMBL::FeaturePair") && defined($f->hseqname)) {
                push(@$fps, $f);
            }
        }
        print STDERR "got '".scalar(@$fps)."' feature pairs\n";

        my $swiss_ids = $self->swissprot_ids_from_features($features);
        my $pfam_ids  = $self->pfam_ids_from_swissprot_ids($swiss_ids);

        my $runnable  =
            Bio::EnsEMBL::Pipeline::Runnable::Finished_HalfwiseHMM->new('-query'    => $genseq,
                                                                        '-features' => $fps,
                                                                        '-hmmdb'    => $self->analysis->db_file,
                                                                        '-pfamids'  => $pfam_ids,
                                                                        #'-pfamdb'   => $self->getPfamDB(),
                                                                        '-options'  => $self->parameters(),
                                                                        '-program'  => $self->analysis->program()
                                                                        );
        #print "created HalfwiseHMM Runnable\n";  
        $self->runnable($runnable);
        #print "finshed fetching input\n";   
    }   
}

=head1 getPfamDB

hack to get a handle to a pfam database, should probably use pfam code

=cut

sub getPfamDB{
    my ($self) = @_;
    unless($self->{'_pfam_db'}){
        my $pfam_meta = $self->db->get_MetaContainer();
	my $value = $pfam_meta->list_value_by_key('pfam_db') || $self->throw("please enter pfam_db key - value into meta table\n");
	my $pfam_db_conn = $self->db->make_hash_from_meta_value($value->[0]);
#        $self->{'_pfam_db'} = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(%{$pfam_meta->get_hash_by_key('pfam_db')});
        $self->{'_pfam_db'} = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(%$pfam_db_conn);
    }
    return $self->{'_pfam_db'};
}

sub check_db_versions{
    my ($self) = @_;
    my $my_pfam = $self->db_version_searched();
    my $pfam_ls = $self->pfam_ls_version();
    unless($pfam_ls eq $my_pfam){
#        $self->throw("VERSION MISMATCH : $pfam_ls not equal to $my_pfam\n");
    }
}

sub db_version_searched{
    my ($self) = @_;    
    unless($self->{'_pfam_db_version'}){
        my $sql = "SELECT pfam_release FROM VERSION LIMIT 1";
        my $pfamDBA = $self->getPfamDB();
        my $sth = $pfamDBA->prepare($sql);
        $sth->execute();
        my ($version) = @{$sth->fetchrow_arrayref};
        $self->{'_pfam_db_version'} = $version;
    }
    return $self->{'_pfam_db_version'};
}

sub pfam_ls_version{
    my ($self) = @_;
    my $debug_this = 1; # this just shows debug info.
    my $force_dbi  = 0; # this will force a dbi call SLOW!!!!!!
    unless($self->{'_pfam_ls_version'}){
        my $db = $self->analysis->db_file;
        $BlastableVersion::debug = $debug_this;            
        warn "BlastableVersion is cvs revision $BlastableVersion::revision \n" if $debug_this;
            
        my $ver = eval { 
            my $blast_ver = BlastableVersion->new();
            $blast_ver->force_dbi($force_dbi); # if set will be SLOW.
            $blast_ver->get_version($db);
            $blast_ver;
        };
        $self->throw("I failed to get a BlastableVersion for $db") if $@;
        
        my $dbv = $ver->version();
        my $sgv = $ver->sanger_version();
        my $name = $ver->name();
        my $date = $ver->date();

        my $message = " - name <" . $name . ">\n" .
                      " - date <" . $date . ">\n" .
                      " - version <" . $dbv . ">\n" .
                      " - sanger_version <" . $sgv . ">\n";
        print STDERR $message if $debug_this;
        unless ($dbv){
            $self->throw("I know nothing about $db I tried to find out:\n$message");
        }
        
        $self->{'_pfam_ls_version'} = $dbv;
    }
    return $self->{'_pfam_ls_version'};
}

=head1  runnable

    Arg      : a Bio::EnsEMBL::Pipeline::RunnableI
    Function : Gets/sets the runnable 
    Exception: throws if argument passed isn't a runnable
    Caller   : 
    Example  :'

=cut  

sub runnable {
    my ($self,$arg) = @_;
 
    if(!defined($self->{'_seqfetchers'})) {
        $self->{'_seqfetchers'} = [];
    }
    
    if (defined($arg)) {
        $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI") unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");
	
        push(@{$self->{_runnable}}, $arg);
    }

    return @{$self->{_runnable}};
}


### 1, -1 and 0 are used to repesent strand, 1 and -1 have the same meaning as in teh ensembl databases 1 is the forward strand and -1 is the reverse 
### 0 is used if a particular id has hit on both strands

sub swissprot_ids_from_features{
    my ($self, $features) = @_;
    $features ||= [];
    #my @features = $self->all_input_features();
    my $swissprot_ids; #{}
    
    foreach my $f(@$features){
        my $hname = $f->hseqname();
        my $strand = $f->strand();
        #print "swissprot id = ".$hname."\n";
        if(!defined($swissprot_ids->{$hname})){
            $swissprot_ids->{$hname} = $strand;
        }else{
            $swissprot_ids->{$hname} = ($strand eq $swissprot_ids->{$hname} ? $strand : 0);
        }
    }
    return $swissprot_ids;
}

=head1 get_pfam_ids

    Arg      : reference to a hash of swissprotids and stands
    Function : gets all pfam ids for a particular swissprot id and puts them in a hash along with the strand the swissprot id is found on and returns its
    Exception: warns if swissprot id has no pfam domains
    Caller   : 

=cut

sub pfam_ids_from_swissprot_ids{
    my ($self, $swissprot_ids) = @_;
    my $pfam_accs;
    my $pfam_lookup;
    # CREATE TEMPORARY TABLE
    my $tbl_name = "pipeline_tmp_$$";
    my $create_table = qq{CREATE TEMPORARY TABLE $tbl_name(
                            pfamseq_acc varchar(12) NOT NULL PRIMARY KEY,
                            strand enum('1','0','-1') DEFAULT '0'
                            )TYPE = HEAP
			}; # should this be HEAP?
                            # There's never gonna be that many matches
                            # to exceed tmp_table_size = 10048576 ???
    my $db = $self->getPfamDB();
    my $sth = $db->prepare($create_table);
    $sth->execute();
    $sth->finish();
    # INSERT
    my (@binds, @values);
    my $sql = qq{INSERT IGNORE INTO $tbl_name (pfamseq_acc, strand) VALUES };
    while (my ($swiss_id, $strand) = each(%$swissprot_ids)){
        #print STDERR "$swiss_id : $strand\n";
	push(@binds, $swiss_id, $strand);
	push(@values, qq{ (?, ?)});
    }
    if(scalar(@values)){
	$sql .= join(", ", @values);
	warn $sql;
	$sth = $db->prepare($sql);
	$sth->execute(@binds);
	$sth->finish();
    }
    # SELECT
    my $select = qq{SELECT CONCAT(a.pfamA_acc, '.',a.version), t.strand, t.pfamseq_acc, a.pfamA_id, a.description
                         FROM pfamA a, pfamA_reg_full f, pfamseq p, $tbl_name t
                         WHERE p.auto_pfamseq = f.auto_pfamseq
                         && a.auto_pfamA = f.auto_pfamA
                         && p.pfamseq_acc = t.pfamseq_acc
                         && f.in_full = 1
                     };
#     my $select = qq{SELECT a.pfamA_acc, t.strand, t.pfamseq_id, a.pfamA_id, p.description
# 			FROM pfamA_reg_full f, pfamseq p, $tbl_name t,  pfamA a
# 			WHERE f.auto_pfamseq = p.auto_pfamseq
# 			&& p.pfamseq_acc     = t.pfamseq_id
# 			&& f.in_full         = 1
# 			&& a.auto_pfamA      = f.auto_pfamA;};
    $sth = $db->prepare($select);
    $sth->execute();
    my ($pfam_acc, $strand, $swall, $pfam_id, $description);
    $sth->bind_columns(\($pfam_acc, $strand, $swall, $pfam_id, $description));
    while(my $row = $sth->fetchrow_arrayref()){
        print STDERR "$swall : $strand : $pfam_acc : $pfam_id\n";
        # making pfam_lookup for later ...
        $pfam_lookup->{$pfam_id} = [$pfam_acc, $description];
        if(defined($pfam_accs->{$pfam_acc})){
            $pfam_accs->{$pfam_acc} = ($strand eq $pfam_accs->{$pfam_acc} ? $strand : 0);
        }else{
            $pfam_accs->{$pfam_acc} = $strand;
        }
    }
    $self->pfam_lookup($pfam_lookup); # store the lookup somewhere useful.
    return $pfam_accs;
}

# get/set for lookup multi-valued hash { pfam_id => [pfam_acc, pfam_desc], ... }
# can append multiple to the lookup (see { %{$self->{'_pfam_lookup'}}, %{$hash_ref} })
sub pfam_lookup{
    my ($self, $hash_ref) = @_;
    if(ref($hash_ref) eq 'HASH'){
        $self->{'_pfam_lookup'} ||= {};
        $self->{'_pfam_lookup'}   = { %{$self->{'_pfam_lookup'}}, %{$hash_ref} };
    }
    return $self->{'_pfam_lookup'};
}

sub run {
    my ($self) = @_;
 
    foreach my $runnable ($self->runnable) {
        $runnable || $self->throw("Can't run - no runnable object");
        print STDERR "using ".$runnable."\n";
        $runnable->run;
    }
    $self->_convert_output();
}

=head1  write_output

    Arg      : none
    Function : writes the converted output to the database as genes
    Exception: none
    Caller   : 
    Example  :

=cut

sub write_output{
    my($self) = @_;
    my @times = times;
    print STDERR "started writing @times \n";
    #$self->_convert_output();
    my @genes    = $self->output();
    my $db       = $self->db();
    
    my $gene_adaptor= $db->get_GeneAdaptor;
    
  GENE: foreach my $gene (@genes) {	
      # do a per gene eval...
      eval {
          #print "gene = ".$gene->type()."\n";
          $gene_adaptor->store($gene);
      };
      if( $@ ) {
          print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
      
  }
    @times = times;
    print STDERR "finished writing @times \n";
    #$self->throw("don't die just yet");
    return 1;
}


=head1  _convert_output

    Arg      : none
    Function : takes the features from the halfwise runnable and runs _make_genes to convert them into Bio::EnsEMBL::Genes with appropriately attached exons and supporting evidence
    Exception: thows if there are no analysis types
    Caller   : 
    Example  :

=cut

sub _convert_output {
    my ($self) = @_;
    #print "converting genes to features\n";
    my @genes;
    my $genetype = 'Halfwise';
    my $anaAdaptor = $self->db->get_AnalysisAdaptor;
    my @analyses = $anaAdaptor->fetch_by_logic_name($genetype);
    my $analysis;
    if(scalar(@analyses) > 1){
        $self->throw("panic! > 1 analysis for $genetype\n");
    }
    elsif(scalar(@analyses) == 1){
        $analysis = $analyses[0];
    }else{
        # make a new analysis object
        $analysis = new Bio::EnsEMBL::Analysis(
                                               -program         => 'genewise',
                                               -program_version => 1,
                                               -gff_source      => 'HalfwiseHMM',
                                               -gff_feature     => 'gene',
                                               -logic_name      => 'Halfwise',
                                               -module          => 'HalfwiseHMM',
                                               );
    }
    # make an array of genes for each runnable
    my @out;
    foreach my $runnable($self->runnable){
        push(@out, $runnable->output);
        $self->pfam_lookup($runnable->pfam_lookup) if $runnable->can('pfam_lookup');
    }
    #print "HalfwiseDB\n";
    #"converting ".scalar(@out)." features to genes\n";
    my @g = $self->_make_genes($genetype, $analysis, \@out);
    push(@genes, @g);
    
    #print STDOUT "genes = @genes\n";
    
    
    if (!defined($self->{_output})) {
        $self->{_output} = [];
    }
    
    push(@{$self->{_output}},@genes);
}

=head1  _make_genes

    Arg      : runnable being run and analysis object being used
    Function : converts the seqfeatures outputed by the runnable and actually converts them into Bio::EnsEMBL::Genes
    Exception: none
    Caller   : 
    Example  :

=cut



=head1 make_genes

  Title   :   make_genes
  Usage   :   $self->make_genes
  Function:   makes Bio::EnsEMBL::Genes out of the output from runnables
  Returns :   array of Bio::EnsEMBL::Gene  
  Args    :   $genetype: string
              $analysis_obj: Bio::EnsEMBL::Analysis
              $results: reference to an array of FeaturePairs

=cut

sub _make_genes {
    my ($self, $genetype, $analysis_obj, $results) = @_;
    my $db           = $self->db();
    my $contig       = $db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    my $dbEntryAdapt = $db->get_DBEntryAdaptor();
    my @genes;
    my $pfam_release = $self->db_version_searched();
#   print "genetype = ".$genetype."\n";
    # fetch lookup multi-valued hash { pfam_id => [pfam_acc, pfam_desc], ... }
    my $pfam_lookup = $self->pfam_lookup();
    my $pfam = "PFAM";
    my $status = 'XREF';
    $self->_check_that_external_db_table_populated($pfam_release, $pfam, $status);

    foreach my $tmp_gene (@$results) {
        my $pfam_id = $tmp_gene->seqname();
        my $dbentry = Bio::EnsEMBL::DBEntry->new(-primary_id  => $pfam_lookup->{$pfam_id}->[0],
                                                 -display_id  => $pfam_id,
                                                 -version     => 1,
                                                 -release     => $pfam_release,
                                                 -dbname      => $pfam,
                                                 -description => $pfam_lookup->{$pfam_id}->[1]
                                                 );
        $dbentry->status($status);
        # Firstly insert the external DB
        my $gene       = Bio::EnsEMBL::Gene->new();;
        my $transcript = $self->_make_transcript($tmp_gene, $contig, $genetype, $analysis_obj);
        $gene->type($genetype);
        $gene->analysis($analysis_obj);
        $gene->add_Transcript($transcript);
        $gene->add_DBEntry($dbentry);
	$gene->display_xref($dbentry);
        $dbEntryAdapt->store($dbentry, $status, $gene);
        push (@genes, $gene);
    }
    
    return @genes;
}

=head1 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::RawContig,
  $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript{
    my ($self, $gene, $contig, $genetype, $analysis_obj) = @_;
    $genetype = 'unspecified' unless defined ($genetype);
    
    unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI")){
        $self->throw("$gene must be Bio::EnsEMBL::SeqFeatureI\n");
    }

    my $transcript   = Bio::EnsEMBL::Transcript->new();
    my $translation  = Bio::EnsEMBL::Translation->new();    
    $transcript->translation($translation);
    
    my $excount = 1;
    my @exons;
    
    foreach my $exon_pred ($gene->sub_SeqFeature) {
        # make an exon
        #print Dumper $exon_pred;
        my $exon = Bio::EnsEMBL::Exon->new();
        
        $exon->id   ($contig->dbID);
        $exon->start($exon_pred->start);
        $exon->end  ($exon_pred->end);
        $exon->strand($exon_pred->strand);
    
        $exon->phase($exon_pred->phase || 0);
        $exon->end_phase(0);
        $exon->attach_seq($contig);
        
        # sort out supporting evidence for this exon prediction
        foreach my $subf($exon_pred->sub_SeqFeature){
            $subf->feature1->seqname($contig->dbID);
            $subf->feature1->score(100);
            $subf->feature1->analysis($analysis_obj);
            
            
            $subf->feature2->score(100);
            $subf->feature2->analysis($analysis_obj);
            
            $exon->add_Supporting_Feature($subf);
        }
        
        push(@exons,$exon);
        
        $excount++;
    }
    
    if ($#exons < 0) {
        print STDERR "Odd.  No exons found\n";
    }else{
        #print STDERR "num exons: " . scalar(@exons) . "\n";
        if ($exons[0]->strand == -1) {
            @exons = sort {$b->start <=> $a->start} @exons;
        }else{
            @exons = sort {$a->start <=> $b->start} @exons;
        }
        
        foreach my $exon (@exons) {
            $transcript->add_Exon($exon);
        }
        
        $translation->start_Exon($exons[0]);
        $translation->end_Exon  ($exons[$#exons]);
        
        if ($exons[0]->phase == 0) {
            $translation->start(1);
        } elsif ($exons[0]->phase == 1) {
            $translation->start(3);
        } elsif ($exons[0]->phase == 2) {
            $translation->start(2);
        }
        
        $translation->end($exons[$#exons]->end - $exons[$#exons]->start + 1);
    }
    
    return $transcript;
}



=head1 output

 Title   : output
 Usage   :
 Function: get/set for output array
 Example :
 Returns : array of Bio::EnsEMBL::Gene
 Args    :

=cut

sub output{
   my ($self,@genes) = @_;
  
   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(scalar(@genes)){
     push(@{$self->{'_output'}},@genes);
   }
   
   return @{$self->{'_output'}};
}

sub _check_that_external_db_table_populated{
    my ($self, $release, $name, $status) = @_;
    $status ||= 'XREF';
    my $db = $self->db();
    my $find_tuple_sql = qq(SELECT count(*) AS tuple_exists 
                              FROM external_db 
                             WHERE db_name = ?
                                && release = ?);
    my $sth = $db->prepare($find_tuple_sql);
    $sth->execute($name, $release);
    my $tuple = $sth->fetchrow_hashref() || {};
    $sth->finish();
    # if there is one return do nothing and the job can get on as normal
    return if $tuple->{'tuple_exists'};
    # else lock table external_db write
    $sth = $db->prepare("LOCK TABLES external_db WRITE");
    $sth->execute();
    $sth->finish();
    # get the next external_db_id
    my $max_db_id_sql = q`SELECT MAX(external_db_id) + 1 AS next_db_id from external_db`;
    $sth = $db->prepare($max_db_id_sql);
    $sth->execute();
    my $row = $sth->fetchrow_hashref || {};
    my $max_db_id = $row->{'next_db_id'} || warn "Error";
    # insert the row
    my $insert_sql = q`INSERT INTO external_db (external_db_id, db_name, release, status) VALUES(?, ?, ?, ?)`;
    $sth = $db->prepare($insert_sql);
    $sth->execute($max_db_id, $name, $release, $status);
    $sth->finish();
    # unlock tables;
    $sth = $db->prepare("UNLOCK TABLES");
    $sth->execute();
    return $max_db_id;
}

1;



