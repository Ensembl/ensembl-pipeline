#
#
# Cared for by ensembl  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Contig_BlastMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Contig_BlastMiniGenewise->new(
                                             -dbobj     => $db,
                                             -input_id  => $id,
                                             -golden_path => $gp
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Contig_BlastMiniGenewise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Data::Dumper;
use Bio::EnsEMBL::Pipeline::BioperlDBConf qw (
                                              BIOPERLDB
                                              BPNAME
                                              BPUSER
                                              BP_DBI_DRIVER
                                              BP_SUPPORTING_DATABASES
                                             );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases   qw (
							       GB_DBHOST
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General     qw (
							       GB_SKIP_BMG
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Similarity  qw (
							       GB_SIMILARITY_DATABASES
							       GB_SIMILARITY_GENETYPE
							       GB_SIMILARITY_MASKING
							       GB_SIMILARITY_SOFTMASK
							      );



@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);  

    my ($path, $type, $threshold) = $self->_rearrange([qw(GOLDEN_PATH TYPE THRESHOLD)], @args);

   
    $self->db->assembly_type($path);

    $self->throw("no protein source databases defined in Config:GeneBuild::Similarity::GB_SIMILARITY_DATABASES\n") 
      unless scalar(@{$GB_SIMILARITY_DATABASES});
    
    foreach my $db(@{$GB_SIMILARITY_DATABASES}){
      my $seqfetcher = $self->make_seqfetcher($db->{'index'}, $db->{seqfetcher});
      $self->add_seqfetcher_by_type($db->{'type'}, $seqfetcher);
    }

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


=head1 RunnableDB implemented methods

=head2 db

    Title   :   db
    Usage   :   $self->db($obj);
    Function:   Gets or sets the value of db
    Returns :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
                (which extends Bio::EnsEMBL::DBSQL::DBAdaptor)
    Args    :   A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

=cut

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=head2 vcontig

 Title   : vcontig
 Usage   : $obj->vcontig($newval)
 Function: 
 Returns : value of vcontig
 Args    : newvalue (optional)

=head1 FPC_BlastMiniGenewise implemented methods

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetchs output data from a frozen perl object
                stored in file $file_name
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub fetch_output {
    my($self,$output) = @_;
    
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

    my $gene_adaptor = $self->db->get_GeneAdaptor;

  GENE: foreach my $gene ($self->output) {      
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
    
    #print STDERR "Fetching input: " . $self->input_id. " \n";
    $self->throw("No input id") unless defined($self->input_id);

    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);

    if(@$GB_SIMILARITY_MASKING){
      my $seq = $contig->get_repeatmasked_seq($GB_SIMILARITY_MASKING, $GB_SIMILARITY_SOFTMASK);
      $self->query($seq);
    }else{
      $self->query($contig);
    }

    #print STDERR "Length is " . $genseq->length . "\n";
    #print STDERR "Fetching features \n";

    #print STDERR "contig: " . $contig . " \n";

    my @features;
    if ($BIOPERLDB) {
          print STDERR "Fetching all HSPs\n";
          my @hsps = $contig->get_all_HSPs;

                

          # _select_features() to pick out the best HSPs with in a region.
      my @features = $self->_select_features (@hsps) unless (scalar(@hsps) ==0); 

          my %bdbs;
      foreach my $feat (@features){
              my $bioperldb = $feat->analysis->db;
              push (@{$bdbs{$bioperldb}},$feat);
      }

      my $bpDBAdaptor = $self->bpDBAdaptor;
      my (@bioperldbs) = split /\,/,$BP_SUPPORTING_DATABASES;

      MINIRUN:foreach my $bioperldb (@bioperldbs){
                
                print STDERR "Creating MiniGenewise for $bioperldb\n";

                $self->seqfetcher($bpDBAdaptor->fetch_BioSeqDatabase_by_name($bioperldb));

                if ($GB_SKIP_BMG) {

                        unless (defined @{$bdbs{$bioperldb}}){
                                print STDERR "Contig has no associated features in $bioperldb\n";
                                next MINIRUN;
                        }        
                                
                        my %scorehash;
                foreach my $f (@{$bdbs{$bioperldb}}) {
                        if (!defined $scorehash{$f->hseqname} || $f->score > $scorehash{$f->hseqname})  {
                                $scorehash{$f->hseqname} = $f->score;
                        }
                }

                        my @forder = sort { $scorehash{$b} <=> $scorehash{$a}} keys %scorehash;

                        my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise('-genomic'    => $self->query,
                                                                     '-features'   => \@{$bdbs{$bioperldb}},
                                                                     '-seqfetcher' => $self->seqfetcher,
                                                                     '-forder'     => \@forder,
                                                                     '-endbias'    => 0);
                $self->runnable($runnable);
                        
                }
                else{
 
                        my @ids;
                        foreach my $f (@{$bdbs{$bioperldb}}){
                                push (@ids, $f->hseqname);
                        }
 
                         my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'    => $self->query,
                                           '-ids'        => \@ids,
                                           '-seqfetcher' => $self->seqfetcher,
                                           '-trim'       => 1);
                $self->runnable($runnable);
 
                }
      }
    }
    else {
      my $alignadaptor = $self->db->get_ProteinAlignFeatureAdaptor();
      foreach my $database(@{$GB_SIMILARITY_DATABASES}){
      @features  = 
	$alignadaptor->fetch_by_Contig_and_score($contig, 
						 $database->{'threshold'}, 
						 $database->{'type'});
      
      #print STDERR "Number of features = " . scalar(@features) . "\n";
      my @filtered_features;
      
      my %idhash;
      
      foreach my $f (@features) {
        if ($f->isa("Bio::EnsEMBL::FeaturePair") && 
            defined($f->hseqname)) {
          $idhash{$f->hseqname} = 1;
        }
      }
    
      my @ids = keys %idhash;
      my $seqfetcher = $self->get_seqfetcher_by_type($database->{'type'});
      #print STDERR "Feature ids are @ids\n";
      my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'    => $self->query,
                                                                             '-ids'        => \@ids,
                                                                             '-seqfetcher' => $seqfetcher,
                                                                             '-trim'       => 1);
      
      
      $self->runnable($runnable);
    }

    # at present, we'll only ever have one ...
    $self->vcontig($contig);
}     

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

=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   calls the run method on each runnable, and then calls convert_output
    Returns :   nothing, but $self->output contains results
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    #Now there is more than one...
    foreach my $runnable ($self->runnable) {
      if ($runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise")){
        $runnable->minirun;
      }else{
        $runnable->run;
      }
    }
    
    $self->convert_output;

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
  my $genetype;
  foreach my $runnable ($self->runnable) {
     if ($runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise") || $runnable->isa("Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise")){
      $genetype = "similarity_genewise";
    }
else{
      $self->throw("I don't know what to do with $runnable");
    }

    my $anaAdaptor = $self->db->get_AnalysisAdaptor;

    #use logic name from analysis object if possible, else take $genetype;
    # if $self->analysis is undefined, this will fall over ...
    my $anal_logic_name;
    if(defined $self->analysis){
      $anal_logic_name = ($self->analysis->logic_name)  ?        $self->analysis->logic_name : $genetype        ;        
    }
    else{
      $anal_logic_name = $genetype;
    }
     
    my @analyses = $anaAdaptor->fetch_by_logic_name($anal_logic_name);
    my $analysis_obj;
    if(scalar(@analyses) > 1){
      $self->throw("panic! > 1 analysis for $genetype\n");
    }
    elsif(scalar(@analyses) == 1){
      $analysis_obj = $analyses[0];
    }
    else{
      # make a new analysis object
      $analysis_obj = new Bio::EnsEMBL::Analysis
        (-db              => 'NULL',
         -db_version      => 1,
         -program         => $genetype,
         -program_version => 1,
         -gff_source      => $genetype,
         -gff_feature     => 'gene',
         -logic_name      => $genetype,
         -module          => 'Contig_BlastMiniGenewise',
      );
    }

    my @results = $runnable->output;
    my @genes = $self->make_genes($genetype, $analysis_obj, \@results);

    $self->output(@genes);

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
  my $contig = $self->vcontig;
  my @tmpf   = @$results;
  my @genes;

  foreach my $tmpf (@tmpf) {
    my $gene       = new Bio::EnsEMBL::Gene;
    my $transcript = $self->_make_transcript($tmpf, $contig, $genetype, $analysis_obj);

    $gene->type($genetype);
    $gene->analysis($analysis_obj);
    $gene->add_Transcript($transcript);

    push (@genes, $gene);
  }

  return @genes;

}

=head2 _make_transcript

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

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}
 

  my $transcript   = new Bio::EnsEMBL::Transcript;
  my $translation  = new Bio::EnsEMBL::Translation;    
  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->contig_id($contig->dbID);
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase($exon_pred->phase);
    $exon->end_phase($exon_pred->end_phase);
    $exon->attach_seq($contig);
    
    # sort out supporting evidence for this exon prediction
    my @sf = $exon_pred->sub_SeqFeature;
    
    my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@sf); 
    
    $align->seqname($contig->dbID);
    $align->score(100);
    $align->analysis($analysis_obj);
    $exon->add_Supporting_Feature($align);
   
    
    push(@exons,$exon);
    
    $excount++;
  }
  
  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
  } 
  else {
    
    #print STDERR "num exons: " . scalar(@exons) . "\n";

    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }
    
    foreach my $exon (@exons) {
      $transcript->add_Exon($exon);
    }
    
    $translation->start_exon($exons[0]);
    $translation->end_exon  ($exons[$#exons]);
    
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }
    
    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
  }
  
  return $transcript;
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
    
   if(defined @genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
}


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

sub add_seqfetcher_by_type{
  my ($self, $type, $seqfetcher) = @_;
  $self->throw("no type specified\n") unless defined ($type); 
  $self->throw("no suitable seqfetcher specified: [$seqfetcher]\n") 
    unless defined ($seqfetcher) && $seqfetcher->isa("Bio::DB::RandomAccessI"); 
  $self->{'_seqfetchers'}{$type} = $seqfetcher;
}


sub get_seqfetcher_by_type{
  my ($self, $type) = @_;
  my %seqfetchers = $self->each_seqfetcher_by_type;
  foreach my $dbtype(keys %seqfetchers){
    if ($dbtype eq $type){
      return $seqfetchers{$dbtype};
    }
  }
}

sub each_seqfetcher_by_type {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB hash of seqfetchers
   return %{$self->{'_seqfetchers'}};
}

sub each_seqfetcher {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB array of seqfetchers
   return values(%{$self->{'_seqfetchers'}});
}


=head2 bpDBAdaptor

  Title   : bpDBAdaptor
  Usage   : $self->bpDBAdaptor($bpDBAdaptor)
  Function: get set the DBAdaptor used by SeqFetcher
  Returns : Bio::DB::SQL::DBAdaptor
  Args    : Bio::DB::SQL::DBAdaptor

=cut


sub bpDBAdaptor {
  my ($self) = @_;
  
  if (defined( $self->{'_bpDBAdaptor'})) {
    print STDERR "Returning a ".$self->{'_bpDBAdaptor'}."\n";
    return $self->{'_bpDBAdaptor'};   
  }
  
  else{
    my $bpname      = $BPNAME || undef;
    my $bpuser      = $BPUSER || undef;
    my $dbhost      = $GB_DBHOST || undef;
    my $DBI_driver  = $BP_DBI_DRIVER || undef;
    my $dbad        = Bio::DB::SQL::DBAdaptor->new(
                                                   -user => $bpuser,
                                                   -dbname => $bpname,
                                                   -host => $dbhost,
                                                   -driver => $DBI_driver,
                                                  );
    $self->{'_bpDBAdaptor'}=$dbad->get_BioDatabaseAdaptor();
    print STDERR "Creating a ".$self->{'_bpDBAdaptor'}."\n";
    
  }
  return $self->{'_bpDBAdaptor'};

}

=head2 _select_features
 
  Title   : _select_features
  Usage   : $self->_select_features(@features)
  Function: obtain the best scoring HSP within a certain area
  Returns : Array of FeaturePairs
  Args    : Array of selected hseqnames
 
=cut

sub _select_features {

        my ($self,@hsps) = @_;

        @hsps = sort {
        $a->strand <=> $b->strand
                    ||
        $a->start <=> $b->start } @hsps;


        my @clusters;
        my $prev = shift @hsps;
        my $hsp_cluster = Bio::EnsEMBL::SeqFeature->new() ;

        $hsp_cluster->add_sub_SeqFeature($prev,'EXPAND');

        push (@clusters,$hsp_cluster);

        foreach my $hsp (@hsps){
        if ($hsp->overlaps($hsp_cluster,'strong')){
                $hsp_cluster->add_sub_SeqFeature($hsp,'EXPAND');
        }
        else{
                $hsp_cluster = Bio::EnsEMBL::SeqFeature->new();
                $hsp_cluster->add_sub_SeqFeature($hsp,'EXPAND');
                push (@clusters,$hsp_cluster);
                 }
        }


        my @selected_hsps;

        foreach my $cluster (@clusters){

        my $new_cluster = Bio::EnsEMBL::SeqFeature->new() ;

        my @hsps = $cluster->sub_SeqFeature;

        @hsps = sort { $b->sub_SeqFeature_Coverage<=> $a->sub_SeqFeature_Coverage} @hsps;

                my $longest_hsp = shift @hsps;

        push (@selected_hsps,$longest_hsp);

        HSP: foreach my $hsp (@hsps){
            my $overlap =0;
            my $missing_exon =0;

        HSP_HIT: foreach my $hsp_hit ($hsp->sub_SeqFeature){
                        my $hit =0;

                 LONG:   foreach my $longest_hit ($longest_hsp->sub_SeqFeature){
                if ($hsp_hit->overlaps($longest_hit)){
                                        $hit =1;
                    my ($overlap_start,$overlap_end);
                    $overlap_start = ($longest_hit->start < $hsp_hit->start) ? $hsp_hit->start : $longest_hit->start;
                    $overlap_end = ($longest_hit->end > $hsp_hit->end) ? $hsp->end : $longest_hit->end;

                    $overlap += $overlap_end - $overlap_start;
                } 
            }
                        $missing_exon = 1 unless ($hit); 
        }
 
        if (($overlap == 0 ) || (($missing_exon)&&( int($hsp->sub_SeqFeature_Coverage/$longest_hsp->sub_SeqFeature_Coverage * 100) > 80))){
            $new_cluster->add_sub_SeqFeature($hsp,'EXPAND');
        }
 
        }
        push (@clusters,$new_cluster) unless scalar($new_cluster->sub_SeqFeature == 0);
        }
 
        my @features;
 
        foreach my $selected_hsp (@selected_hsps){
        push (@features,$selected_hsp->sub_SeqFeature);
        }
        return @features;
      }
}

1;
