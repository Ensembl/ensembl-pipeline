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

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBNAME
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
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
							      
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted  qw (
							     GB_TARGETTED_GW_GENETYPE
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General    qw (
							     GB_INPUTID_REGEX
							     GB_REPEAT_MASKING
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
      my $seqfetcher =  $self->make_seqfetcher($db->{'index'}, $db->{seqfetcher});  
      $self->add_seqfetcher_by_type($db->{'type'}, $seqfetcher);
    }

    # IMPORTANT
    # SUPER creates db, which is a reference to GB_DBHOST@GB_DBNAME containing
    # features and dna
    # Here it is used as refdb only and we need to make a connection to GB_GW_DBNAME@GB_GW_DBHOST
    my $genewise_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
							 '-host'   => $GB_GW_DBHOST,
							 '-user'   => $GB_GW_DBUSER,
							 '-pass'   => $GB_GW_DBPASS,
							 '-dbname' => $GB_GW_DBNAME,
							);
    
    
    $genewise_db->dnadb($self->db);
    $self->output_db($genewise_db);
    
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

    my $gene_adaptor = $self->output_db->get_GeneAdaptor;
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
    
    print STDERR "Fetching input id : " . $self->input_id. " \n\n";

    $self->throw("No input id") unless defined($self->input_id);
    my $input_id = $self->input_id;
   
    $input_id =~ /$GB_INPUTID_REGEX/;
    my $chrid = $1;
    my $chrstart = $2;
    my $chrend   = $3;

    
    my $sla       = $self->db->get_SliceAdaptor();
    my $slice     = $sla->fetch_by_chr_start_end($chrid,$chrstart,$chrend);
    $slice->chr_name($chrid);
    if(@$GB_SIMILARITY_MASKING){
      my $seq = $slice->get_repeatmasked_seq($GB_SIMILARITY_MASKING, $GB_SIMILARITY_SOFTMASK);
      $self->query($seq);
    }else{
      $self->query($slice);
    }
    
    

    my @genes     = @{$slice->get_all_Genes_by_type($GB_TARGETTED_GW_GENETYPE)};
    
    DATABASE: foreach my $database(@{$GB_SIMILARITY_DATABASES}){
      
      print STDERR "Fetching features for " . $database->{'type'} . 
	" with score above " . $database->{'threshold'}. "\n\n";
      my $pafa = $self->db->get_ProteinAlignFeatureAdaptor();
      #print STDERR "Fetching features from db ".$self->db->dbname." on ".$self->db->host."\n";
      my @features  = @{$pafa->fetch_all_by_Slice_and_score($slice, $database->{'threshold'}, $database->{'type'})};
      print STDERR "have ".@features." \n";
      # lose version numbers - probably temporary till pfetch indices catch up
      if(!@features){
	print STDERR "have zero features for ".$database->{'type'}." exiting\n";
	next DATABASE;
      }
      foreach my $f(@features) {
	my $name = $f->hseqname;
	if ($name =~ /(\S+)\.\d+/) { 
	  $f->hseqname($1);
	}
      }
      
      #print STDERR "Number of features = " . scalar(@features) . "\n\n";
      my %redids;
      my $trancount = 1;
      
      # check which TargettedGenewise exons overlap similarity features
      foreach my $gene (@genes) {
	foreach my $tran (@{$gene->get_all_Transcripts}) {
	  foreach my $exon (@{$tran->get_all_Exons}) {
	    if ($exon->seqname eq $slice->id) {
	      
	    FEAT: foreach my $f (@features) {
		if ($exon->overlaps($f)) {
		  $redids{$f->hseqname} = 1;
		  #		print STDERR "ID " . $f->hseqname . " covered by genewise\n";
		}
	      }
	    }
	  }
	  $trancount++;
	}
      }
      
      my %idhash;
      
      # collect those features which haven't been used by Targetted GeneWise
      foreach my $f (@features) {
	      #print "Feature : " . $f->gffstring . "\n";
	
	if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
	    defined($f->hseqname) &&
	    $redids{$f->hseqname} != 1) {
	  $idhash{$f->hseqname} = 1;
	  
	}
      }
      

      my @ids = keys %idhash;

      my $seqfetcher =  $self->get_seqfetcher_by_type($database->{'type'});
      #print STDERR "Feature ids are @ids\n";
      
      my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $self->query,
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
    #print STDERR "HAVE CONVERTED OUTPUT\n";
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
    print STDERR $anaAdaptor . "\n";

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

  foreach my $tmpf (@{$results}) {
    my $unchecked_transcript = 
      Bio::EnsEMBL::Pipeline::Tools::GeneUtils->SeqFeature_to_Transcript($tmpf, 
									 $self->query, 
									 $analysis_obj, 
									 $self->output_db, 
									 0);
    
    next unless defined ($unchecked_transcript);

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
    
    # make genes from valid transcripts
    foreach my $checked_transcript(@$valid_transcripts){
      #print "have ".$checked_transcript." as a transcript\n";
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

=head2 output_db

 Title   : output_db
 Usage   : needs to be moved to a genebuild base class
 Function: 
           
 Returns : 
 Args    : 

=cut

sub output_db {
    my( $self, $output_db ) = @_;
    
    if ($output_db) 
    {
	$output_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
	    || $self->throw("Input [$output_db] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
	$self->{_output_db} = $output_db;
    }
    return $self->{_output_db};
}

1;
