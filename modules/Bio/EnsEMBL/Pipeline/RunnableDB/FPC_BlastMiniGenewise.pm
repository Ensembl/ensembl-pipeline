##
#
# Cared for by Ensembl  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise

=head1 SYNOPSIS

my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::MiniGenewise->new(
					     -db        => $db,
					     -input_id  => $id,
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniGenewise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Analysis;


use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBNAME
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							     GB_GW_DBPORT
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Similarity qw (
							      GB_SIMILARITY_DATABASES
							      GB_SIMILARITY_COVERAGE
							      GB_SIMILARITY_MAX_INTRON
							      GB_SIMILARITY_MIN_SPLIT_COVERAGE
							      GB_SIMILARITY_GENETYPE
							      GB_SIMILARITY_MAX_LOW_COMPLEXITY 
							      GB_SIMILARITY_MASKING
							      GB_SIMILARITY_SOFTMASK
							      GB_SIMILARITY_GENETYPEMASKED
							      GB_SIMILARITY_BLAST_FILTER
							      
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted  qw (
							     GB_TARGETTED_GW_GENETYPE
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General    qw (
							     GB_INPUTID_REGEX
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts    qw (
							     GB_KILL_LIST
							    );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
      
    # make sure at least one protein source database has been defined
    
    $self->throw("no protein source databases defined in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES\n") 
      unless scalar(@{$GB_SIMILARITY_DATABASES});
    # make all seqfetchers
    foreach my $db(@{$GB_SIMILARITY_DATABASES}){
      my $type = $db->{"type"};
      my $seqfetcher =  $self->make_seqfetcher($db->{index}, $db->{seqfetcher});  
      $self->add_seqfetcher_by_type($type, $seqfetcher);
    }

    # IMPORTANT
    # SUPER creates db, which is a reference to GB_DBHOST@GB_DBNAME containing
    # features and dna
    # Here it is used as refdb only and we need to make a connection to GB_GW_DBNAME@GB_GW_DBHOST

    $GB_GW_DBHOST = $self->db->host     if (!defined($GB_GW_DBHOST) || $GB_GW_DBHOST eq '');
    $GB_GW_DBUSER = $self->db->username if (!defined($GB_GW_DBUSER) || $GB_GW_DBUSER eq '');
    $GB_GW_DBPASS = $self->db->password if (!defined($GB_GW_DBPASS) || $GB_GW_DBPASS eq '');
    $GB_GW_DBNAME = $self->db->dbname   if (!defined($GB_GW_DBNAME) || $GB_GW_DBNAME eq '');
    $GB_GW_DBPORT = $self->db->dbname   if (!defined($GB_GW_DBPORT) || $GB_GW_DBPORT eq '');


    my $genewise_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
							 '-host'   => $GB_GW_DBHOST,
							 '-user'   => $GB_GW_DBUSER,
							 '-pass'   => $GB_GW_DBPASS,
							 '-port'   => $GB_GW_DBPORT,
							 '-dbname' => $GB_GW_DBNAME,
							);
    
    
    $genewise_db->dnadb($self->db);
    $self->genewise_db($genewise_db);
    
    
    return $self; 
  }


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self,@features) = @_;

    my $gene_adaptor = $self->genewise_db->get_GeneAdaptor;
    my @genes = $self->output;
    print STDERR "have ".@genes." genes\n";
    GENE: foreach my $gene (@genes) {	
	# do a per gene eval...
	eval {
	    $gene_adaptor->store($gene);
	    print STDERR "wrote gene " . $gene->dbID . "\n";
	}; 
	if( $@ ) {
	    print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
	}
	
    }
   
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :   none

=cut

  sub fetch_input {
    my( $self) = @_;
    
    print STDERR "Fetching input id : " . $self->input_id. " \n";

    $self->throw("No input id") unless defined($self->input_id);
    my $input_id = $self->input_id;
   
    $input_id =~ /$GB_INPUTID_REGEX/;

    my $chrid    = $1;
    my $chrstart = $2;
    my $chrend   = $3;

    
    my $sla       = $self->genewise_db->get_SliceAdaptor();
    my $slice     = $sla->fetch_by_chr_start_end($chrid,$chrstart,$chrend);
    my $seq;
    $self->query($slice);
    $slice->chr_name($chrid);
    if(@$GB_SIMILARITY_MASKING){
      $seq = $slice->get_repeatmasked_seq($GB_SIMILARITY_MASKING, $GB_SIMILARITY_SOFTMASK);
    }else{
      $seq = $slice;
    }
    
    # No similarity genes will be build overlapping the gene type put in 
    # the list 'GB_SIMILARITY_GENETYPEMASKED'. If nothing is put, the default 
    # will be to take Targetted genewise genes

    my @genes;

    if (@{$GB_SIMILARITY_GENETYPEMASKED}) {
	foreach my $type (@{$GB_SIMILARITY_GENETYPEMASKED}) {
	    print STDERR "FETCHING GENE TYPE : $type\n";
	    foreach my $mask_genes (@{$slice->get_all_Genes_by_type($type)}) {
		push (@genes,$mask_genes);
	    }
	}
    }
    else {
	print STDERR "FETCHING GENE TYPE : $GB_TARGETTED_GW_GENETYPE\n";
	@genes     = @{$slice->get_all_Genes_by_type($GB_TARGETTED_GW_GENETYPE)};
    }
            
    DATABASE: foreach my $database(@{$GB_SIMILARITY_DATABASES}){
      
	printf(STDERR "Fetching features for %s with score above %d from %s\@%s...", 
	       $database->{'type'}, 
	       $database->{'threshold'}, 
	       $self->db->dbname, 
	       $self->db->host);
	
	my $pafa = $self->db->get_ProteinAlignFeatureAdaptor();
	my @features  = @{$pafa->fetch_all_by_Slice_and_score($slice, 
							      $database->{'threshold'}, 
							      $database->{'type'})};
	print STDERR "got " . @features." \n";
	next DATABASE if not @features;

	foreach my $f(@features) {
	    my $name = $f->hseqname;
	    if ($name =~ /(\S+)\.\d+/) { 
		$f->hseqname($1);
	    }
	}

	# check which TargettedGenewise exons overlap similarity features

	my @all_exons;
	foreach my $gene (@genes) {
	    foreach my $tran (@{$gene->get_all_Transcripts}) {
		foreach my $exon (@{$tran->get_all_Exons}) {
		    if ($exon->seqname eq $slice->id) {
			push @all_exons, $exon;
		    }
		}
	    }
	}
	@all_exons = sort {$a->start <=> $b->start} @all_exons;

	my %redids;
	my $ex_idx = 0;
	
        FEAT: foreach my $f (sort {$a->start <=> $b->start} @features) {
	    for( ; $ex_idx < @all_exons; ) {
		my $exon = $all_exons[$ex_idx];
		
		if ($exon->overlaps($f)) {
		    $redids{$f->hseqname} = 1;
		    next FEAT;
		}
		elsif ($exon->start > $f->end) {
		    # no exons will overlap this feature
		    next FEAT;
		}
		else {
		    $ex_idx++;
		}
	    }
	}

	# collect those features which haven't been used by Targetted Genewise
	my (%idhash);
	# reject them if they appear in the kill list
	my %kill_list = %{$self->fill_kill_list};
	
	foreach my $f (@features) {
	    # print "Feature : " . $f->gffstring . "\n";

	    if (defined $kill_list{$f->hseqname}){
		print STDERR "skipping over " . $f->hseqname . " : in kill list\n";
		next;
	    }	

	    if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
		defined($f->hseqname) &&
		$redids{$f->hseqname} != 1) {
		
		push(@{$idhash{$f->hseqname}},$f);
	    }
	}
            
	my @ids = keys %idhash;
	
	if ($GB_SIMILARITY_BLAST_FILTER) {
	    print STDERR "Performing Anopheles-stle filtering...\n";
	    print "  No. IDS BEFORE FILTER: ".scalar(@ids)."\n"; 
	    
	    my @sortedids = $self->sort_hids_by_coverage($database,\%idhash);
	    my @newids = $self->prune_features($seq,\@sortedids,\%idhash);
	    
	    print "  No. IDS AFTER FILTER: ".scalar(@newids)."\n"; 
	    
	    @ids = @newids;
	}

	my $seqfetcher =  $self->get_seqfetcher_by_type($database->{'type'});
	#print STDERR "Feature ids are @ids\n";
	print STDERR "have ".@ids." ids to pass to BMG\n"; 
	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $seq,
									       '-ids'      => \@ids,
									       '-seqfetcher' => $seqfetcher,
									       '-trim'     => 1);
	
	$self->runnable($runnable);	
	# at present, we'll only ever have one ...       	
    }
}    
  

=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   calls the run method on each runnable, and then calls convert_output
    Returns :   nothing, but $self->output contains results
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    # is there ever going to be more than one?
    foreach my $runnable ($self->runnable) {
	$runnable->run;
    }
    
    $self->convert_output;
    # print STDERR "HAVE CONVERTED OUTPUT\n";
}

=head2 convert_output

  Title   :   convert_output
  Usage   :   $self->convert_output
  Function:   converts output from each runnable into gene predictions
  Returns :   nothing, but $self->output contains results
  Args    :   none

=cut

sub convert_output {
  my ($self) =@_;
  my $trancount = 1;
  my $genetype = $GB_SIMILARITY_GENETYPE;
  #print STDERR "in CONVERT OUTPUT\n";
  foreach my $runnable ($self->runnable) {
    $self->throw("I don't know what to do with $runnable") unless $runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise");
										 
    if(!defined($genetype) || $genetype eq ''){
      $genetype = 'similarity_genewise';
      $self->warn("Setting genetype to $genetype\n");
    }

    my $anaAdaptor = $self->db->get_AnalysisAdaptor;
    # print STDERR $anaAdaptor . "\n";

    my $analysis_obj = $anaAdaptor->fetch_by_logic_name($genetype);

    if ( !defined $analysis_obj ){
      # make a new analysis object
      $analysis_obj = new Bio::EnsEMBL::Analysis
	(-db              => 'NULL',
	 -db_version      => 1,
	 -program         => $genetype,
	 -program_version => 1,
	 -gff_source      => $genetype,
	 -gff_feature     => 'gene',
	 -logic_name      => $genetype,
	 -module          => 'FPC_BlastMiniGenewise',
      );
    }

    my @results = $runnable->output;
    my $genes = $self->make_genes($genetype, $analysis_obj, \@results);

    my @remapped = @{$self->remap_genes($genes)};
    #print STDERR "HAVE ".@remapped." remapped genes\n";
    $self->output(@remapped);

  }
}


=head2 make_genes

  Title   :   make_genes
  Usage   :   $self->make_genes
  Function:   makes Bio::EnsEMBL::Genes out of the output from runnables
  Returns :   array of Bio::EnsEMBL::Gene  
  Args    :   $genetype: string
              $analysis_obj: Bio::EnsEMBL::Analysis
              $results: reference to an array of FeaturePairs

=cut

sub make_genes {
  my ($self, $genetype, $analysis_obj, $results) = @_;
  my @genes;

  $self->throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") 
    unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");

  my $count = 0;
  foreach my $tmpf (@{$results}) {
      $count++;

      my $unchecked_transcript = 
	  Bio::EnsEMBL::Pipeline::Tools::GeneUtils->SeqFeature_to_Transcript($tmpf, 
									     $self->query, 
									     $analysis_obj, 
									 $self->genewise_db, 
									     0);

      if (not defined $unchecked_transcript) {
	  print STDERR "   Transcript $count : REJECTED because could not make from SeqFeatures\n";
	  next;
      }

      
      # validate transcript
      my @seqfetchers = $self->each_seqfetcher;

      my $valid_transcripts = 
	  Bio::EnsEMBL::Pipeline::Tools::GeneUtils->validate_Transcript($unchecked_transcript, 
									$self->query, 
									$GB_SIMILARITY_COVERAGE, 
									$GB_SIMILARITY_MAX_INTRON, 
									$GB_SIMILARITY_MIN_SPLIT_COVERAGE, 
									\@seqfetchers, 
									$GB_SIMILARITY_MAX_LOW_COMPLEXITY);
      if (not defined $valid_transcripts) {
	  print STDERR "   Transcript $count : REJECTED because validation failed\n";
	  next;
      }
      
      # make genes from valid transcripts
      foreach my $checked_transcript(@$valid_transcripts){
	  #print "have ".$checked_transcript." as a transcript\n";
	  
	  # add terminal stop codon if appropriate
	  $checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon($checked_transcript);

	  my $gene = new Bio::EnsEMBL::Gene;
	  $gene->type($genetype);
	  $gene->analysis($analysis_obj);
	  $gene->add_Transcript($checked_transcript);
	  
	  push (@genes, $gene);
      }
  }
  
  return \@genes;

}

=head2 remap_genes

 Title   : remap_genes
 Usage   : $self->remap_genes($runnable, @genes)
 Function: converts the coordinates of each Bio@EnsEMBL::Gene in @genes into RawContig
           coordinates for storage.
 Example : 
 Returns : array of Bio::EnsEMBL::Gene in RawContig coordinates
 Args    : @genes: array of Bio::EnsEMBL::Gene in virtual contig coordinates


=cut

sub remap_genes {
  my ($self, $genes) = @_;
  my $contig = $self->query;

  my @newf;
  my $trancount=1;
  foreach my $gene (@$genes) {
    eval {
      $gene->transform;
      # need to explicitly add back genetype and analysis.
      #$newgene->type($gene->type);
      #$gene->analysis($gene->analysis);

      foreach my $tran (@{$gene->get_all_Transcripts}) {
	foreach my $exon(@{$tran->get_all_Exons}) {
	  foreach my $sf(@{$exon->get_all_supporting_features}) {
	    # this should be sorted out by the remapping to rawcontig ... strand is fine
	    if ($sf->start > $sf->end) {
	      my $tmp = $sf->start;
	      $sf->start($sf->end);
	      $sf->end($tmp);
	    }
	  }
	}
      }
      push(@newf,$gene);

    };
    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene . " [$@]\n";
    }
    

  }

  return \@newf;
}

=head2 output

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
    
   if(@genes){
     push(@{$self->{'_output'}},@genes);
   }
   #print STDERR "have ".@{$self->{'_output'}}." as output\n";
   return @{$self->{'_output'}};
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: makes a Bio::EnsEMBL::SeqFetcher to be used for fetching protein sequences. If 
           $index is specified, then a Getseqs fetcher is made, otherwise it throws
 Example :
 Returns : Bio::EnsEMBL::SeqFetcher
 Args    : string representing path to sequence index


=cut

sub make_seqfetcher {
  my ( $self, $index, $seqfetcher_class ) = @_;
  my $seqfetcher;

  (my $class = $seqfetcher_class) =~ s/::/\//g;
  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );
    
    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );
    
  }
  else{
    $self->throw("Can't make seqfetcher\n");
  }

  return $seqfetcher;

}

=head2 each_seqfetcher

  Title   :   each_seqfetcher
  Usage   :   my @seqfetchers = $self->each_seqfetcher
  Function:   Returns an array of Bio::DB::RandomAccessI representing the various sequence indices 
              listed in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES
  Returns :   Array of Bio::DB::RandomAccessI
  Args    :   none

=cut

sub each_seqfetcher {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB array of seqfetchers
   return values(%{$self->{'_seqfetchers'}});
}

=head2 each_seqfetcher_by_type

  Title   :   each_seqfetcher_by_type
  Usage   :   my %seqfetchers_by_type = $self->each_seqfetcher_by_type
  Function:   Returns a hash of Bio::DB::RandomAccessI representing the various sequence indices 
              listed in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES keyed by type listed therein.
  Returns :   Hash of all seqfetchers linking db_type to Bio::DB::RandomAccessI
  Args    :   none

=cut

sub each_seqfetcher_by_type {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB hash of seqfetchers
   return %{$self->{'_seqfetchers'}};
}

=head2 add_seqfetcher_by_type

  Title   :   add_seqfetcher_by_type
  Usage   :   $self->add_seqfetcher_by_type('swall', $seqfetcher)
  Function:   Adds a Bio::DB::RandomAccessI into $self->{'_seqfetchers'} keyed by type
  Returns :   Nothing
  Args    :   $type - string representing db type
              $seqfetcher - Bio::DB::RandomAccesI

=cut

sub add_seqfetcher_by_type{
  my ($self, $type, $seqfetcher) = @_;
  $self->throw("no type specified\n") unless defined ($type); 
  $self->throw("no suitable seqfetcher specified: [$seqfetcher]\n") 
    unless defined ($seqfetcher) && $seqfetcher->isa("Bio::DB::RandomAccessI"); 
  $self->{'_seqfetchers'}{$type} = $seqfetcher;
}


=head2 get_seqfetcher_by_type

  Title   :   get_seqfetcher_by_type
  Usage   :   my $seqfetcher = $self->get_seqfetcher_by_type('swall')
  Function:   Fetches the seqfetcher associated with a particular db type as specified in Config::GeneBuild::Similarity::GB_SIMILARITY_DATABASES
  Returns :   Bio::DB::RandomAccessI
  Args    :   $type - string representing db type

=cut
sub get_seqfetcher_by_type{
  my ($self, $type) = @_;
  my %seqfetchers = $self->each_seqfetcher_by_type;
  foreach my $dbtype(keys %seqfetchers){
    if ($dbtype eq $type){
      return $seqfetchers{$dbtype};
    }
  }
}





=head2 fill_kill_list

 Title   : fill_kill_list
 Usage   : 
 Function: 
           
 Returns : 
 Args    : 

=cut

sub fill_kill_list {
  my ($self) = @_;
  my %kill_list;

  if (defined($GB_KILL_LIST) && $GB_KILL_LIST ne '') {
  open (KILL_LIST, "< $GB_KILL_LIST") or die "can't open $GB_KILL_LIST";
  while (<KILL_LIST>) {

    chomp;
    my @list = split;
    next unless scalar(@list); 	# blank or empty line
    $kill_list{$list[0]} = 1;
  }

  close KILL_LIST or die "error closing $GB_KILL_LIST\n";
  }
  return \%kill_list;
}


=head2 sort_hids_by_coverage

 Title   : sort_hids_by_coverage
 Usage   : my @sorted_ids = $self->sort_hids_by_coverage($database,\%hids)
 Function: This is supposed to order each hid in function of there coverage with the sequence they            have been build with
 Example :
 Returns : an array of ordered hids
 Args    : Database description(where to fetch the feature from). A hash refererence containing hid           and its features.


=cut

sub sort_hids_by_coverage{
   my ($self,$database,$hash_ref) = @_;
   my @id =  keys %{$hash_ref};
   my $seqfetcher =  $self->get_seqfetcher_by_type($database->{'type'});
   my $forward_start;
   my $forward_end;
   my $reverse_start;
   my $reverse_end;
   my $f_matches = 0;
   my $r_matches = 0;
   my $protname;
   my $seq;
   my $plength;
   my $features;
   my $best_cov;
   my %idsreturned;
   my @sorted;

 IDS: foreach my $id (@id) {

     $protname = undef;
     $seq = undef;
     $plength = 0;
     $forward_start = 0;
     $forward_end = 0;
     $reverse_start = 0;
     $reverse_end = 0;
     $best_cov = 0;
     $features = $hash_ref->{$id};
     

   FEAT: foreach my $feat(@$features) {
	  if ($feat->strand == 1) {
	      if (!defined($protname)){
		  $protname = $feat->hseqname;
	      }
	      if($protname ne $feat->hseqname){
		  warn ("$protname ne " . $feat->hseqname . "\n");
	      }
	      
	      
      
	      if((!$forward_start) || $forward_start > $feat->hstart){
		  $forward_start = $feat->hstart;
	      }
	      
	      if((!$forward_end) || $forward_end < $feat->hend){
		  $forward_end= $feat->hend;
	      }
	  }

	  if ($feat->strand == -1) {
	      if (!defined($protname)){
		  $protname = $feat->hseqname;
	      }
	      if($protname ne $feat->hseqname){
		  warn ("$protname ne " . $feat->hseqname . "\n");
	      }
	      
	      
      
	      if((!$reverse_start) || $reverse_start > $feat->hstart){
		  $reverse_start = $feat->hstart;
	      }
	      
	      if((!$reverse_end) || $reverse_end < $feat->hend){
		  $reverse_end= $feat->hend;
	      }
	  }
	      
      }

 
     $f_matches = ($forward_end - $forward_start + 1);
     $r_matches = ($reverse_end - $reverse_start + 1);

     if($self->get_length_by_id($protname)) {
	 $plength = $self->get_length_by_id($protname);
     }

     else {
     
       SEQFETCHER:
	 foreach my $seqfetcher( $self->each_seqfetcher){
	     eval{
		 $seq = $seqfetcher->get_Seq_by_acc($protname);
	     };
	 
	     if ($@) {
		 $self->warn("FPC_BMG:Error fetching sequence for [$protname] - trying next seqfetcher:[$@]\n");
	     }
	 
	     if (defined $seq) {
		 last SEQFETCHER;
	     }
	 
	 }
     
	 if(!defined $seq){
	     $self->warn("FPC_BMG:No sequence fetched for [$protname] - can't check coverage - set coverage to 1 (entry will be considered as having low coverage...sorry\n");
	 }
     
	 eval {
	     $plength = $seq->length
	     };

	 if($plength) {
	     $self->add_length_by_id($protname,$plength);
	 }
     }	 

     my $f_coverage;
     my $r_coverage;
     
#check the low complexity
     my $valid = 0;
     if ($seq) {
	 
	 eval {
	     # Ugh! 
	     my $analysis = Bio::EnsEMBL::Analysis->new(
							-db           => 'low_complexity',
							-program      => '/usr/local/ensembl/bin/seg',
							-program_file => '/usr/local/ensembl/bin/seg',
							-gff_source   => 'Seg',
							-gff_feature  => 'annot',
							-module       => 'Seg',
							-logic_name   => 'Seg'
							
						    );
	     
	     my $seg = new  Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg(    
									       -query    => $seq,
									       -analysis => $analysis,
									   );
	     
	     $seg->run;
	     
	     
	     if($seg->get_low_complexity_length > $GB_SIMILARITY_MAX_LOW_COMPLEXITY){
		 warn("discarding feature too much low complexity\n");
		 $valid = 0;
	     }
	 $valid = 1;
	 
	 
	 };
     
	 if($@){
	     print STDERR "problem running seg: \n[$@]\n";
	     $valid = 1;		# let transcript through
	 }
     }
 

#




   if(!defined($plength) || $plength == 0 || $valid == 0){
       warn("no sensible length for $protname - can't get coverage - or too much low complexity\n");
       
       #Can't check the length of the protein, set the coverage to 1 for both strand
       $f_matches = 1;
       $r_matches = 1;
       $plength = 100;
   }


   $f_coverage = $f_matches/$plength;
   $r_coverage = $r_matches/$plength;
   $f_coverage *= 100;
   $r_coverage *= 100;
     
     if ($f_coverage > $r_coverage) {
	 $best_cov = $f_coverage;
     }
     else {
	 $best_cov = $r_coverage;
     }

     $idsreturned{$protname} = $best_cov;
 }
   
   @sorted = sort { $idsreturned{$b} <=> $idsreturned{$a} }  keys %idsreturned;
   
   return @sorted;
}


=head2 prune_features

 Title   : prune_features
 Usage   : my @pruned_ids = $self->prune_features($genseq,\@sortedids,\%idhash);
 Function: This method as for goal to limit the number of hids sent to BlastMiniGenewise.
 Example :
 Returns : An array of pruned hids
 Args    : genomic sequence, an array of sorted hids, an array reference of all the hids and their             features


=cut

sub prune_features{
    my ($self,$genseq,$sortedids_array_ref,$hash_ref) = @_;
    my $forwardcountstring = '0' x $genseq->length;
    my $reversecountstring = '0' x $genseq->length;
    my %ids;
    my %hidcount;
    my %finalids;
    
    my @ident = @$sortedids_array_ref;

  IDS: foreach my $id (@ident) {
      my $features = $hash_ref->{$id};
      
      my @exons;
      
      next IDS unless (ref($features) eq "ARRAY");
      
      next IDS unless (scalar(@$features) >= 1);
      
    FEAT: foreach my $feat(@$features) {
	my $length = $feat->end - $feat->start + 1;
		
	if($feat->strand == 1) {
	    my $count = 0;
	    
	    my $str = substr($forwardcountstring,$feat->start,$length);
	    
	    foreach my $byte (split //, $str) {
		$count = $count + $byte;
	    }
	   

	    $hidcount{$id}->{'lengthforward'} = $hidcount{$id}->{'lengthforward'} + $length;
	    $hidcount{$id}->{'coverageforward'} = $hidcount{$id}->{'coverageforward'} + $count;
	   
	     
	    if ($count > 10) {
		next FEAT;
	    }
	   
	    else {
		substr($forwardcountstring,$feat->start,$length) = '1' x $length; 
		$ids{$feat->hseqname} = 1;
	    }
	}
	   
	elsif ($feat->strand == -1) {
	    my $count = 0;
	    my $length = $feat->end - $feat->start + 1;
	    
	   

	    my $str = substr($reversecountstring,$feat->start,$length);
	       
	    foreach my $byte (split //, $str) {
		$count = $count + $byte;
	    }
	    $hidcount{$id}->{'lengthreverse'} = $hidcount{$id}->{'lengthreverse'} + $length;
	    $hidcount{$id}->{'coveragereverse'} = $hidcount{$id}->{'coveragereverse'} + $count;
	    
	    if ($count > 10) {
		next FEAT;
	    }

	    else {
		substr($reversecountstring,$feat->start,$length) = '1' x $length; 
		
		$ids{$feat->hseqname} = 1;

	    }
	}
   }
  }
    
    foreach my $hidk (keys %ids) {
	my $percforward;
	my $percreverse;
	
	if ($hidcount{$hidk}->{'lengthforward'} > 0) {
	    $percforward = $hidcount{$hidk}->{'coverageforward'} * 100 / $hidcount{$hidk}->{'lengthforward'};
	}
	
	if ($hidcount{$hidk}->{'lengthreverse'} > 0) {
	    $percreverse = $hidcount{$hidk}->{'coveragereverse'} * 100 / $hidcount{$hidk}->{'lengthreverse'};
	}

	if (($percforward <= 90)&&($percreverse <= 90)) {
	    $finalids{$hidk} = 1;
	}
    }
   
    my @returnids = keys %finalids;
    return @returnids;
}


=head2 add_length_by_id

 Title   : add_length_by_id
 Usage   : $self->add_length_by_id("QE345",345)
 Function: Caches the length of a protein
 Example :
 Returns : nothing
 Args    : protein id and its length


=cut

sub add_length_by_id{
   my ($self,$id,$length) = @_;
   $self->throw("no id specified\n") unless defined ($id); 
   $self->throw("no length specified\n") unless defined ($length);
   
   if (!defined $self->{'_idlength'}{$id}) {
       $self->{'_idlength'}{$id} = $length;
   }
}

=head2 get_length_by_id

 Title   : get_length_by_id
 Usage   : $self->get_length_by_id("QGY7980")
 Function: retrieves the length of a given protein
 Example :
 Returns : Length of a protein
 Args    : Protein id


=cut

sub get_length_by_id{
   my ($self,$id) = @_;
   my %idlength = $self->{'_idlength'};
   if ($idlength{$id}) {
       return $idlength{$id};
   }
}

=head2 genewise_db

 Title   : genewise_db
 Usage   : retrieves /sets the db pointed at by host GB_GW_HOST -
           where preliminary genes are being written
 Function:

 Returns :
 Args    :

=cut

sub genewise_db {
    my( $self, $genewise_db ) = @_;

    if ($genewise_db)
    {
        $genewise_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
            || $self->throw("Input [$genewise_db] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
        $self->{_genewise_db} = $genewise_db;
    }
    return $self->{_genewise_db};
}




1;



