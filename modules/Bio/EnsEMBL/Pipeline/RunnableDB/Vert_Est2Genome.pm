#!/usr/local/bin/perl

#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome->new(
					     -dbobj     => $db,
					     -input_id  => $id
                                             );
    $obj->fetch_input;
    $obj->run;

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::AlignFeature;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::EnsEMBL::Analysis::MSPcrunch;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::SeqIO;
use Bio::Tools::Blast;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  
  $self->{'_fplist'} = []; #create key to an array of feature pairs
  
  my( $dbobj, $input_id ) = $self->_rearrange(['DBOBJ',
					       'INPUT_ID'], @args);
  
  $self->throw("No database handle input")           
    unless defined($dbobj);
  $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI") 
    unless $dbobj->isa("Bio::EnsEMBL::DB::ObjI");
  
  $self->dbobj($dbobj);
    
  $self->throw("No input id input") 
    unless defined($input_id);
  $self->input_id($input_id);
  
  return $self; # success - we hope!
}

=head2 features

    Title   :   features
    Usage   :   $self->features(@features)
    Function:   Get/set method for features list
    Returns :   Array of Bio::Ensembl::SeqFeature
    Args    :   Array of Bio::Ensembl::SeqFeature

=cut

sub features {
    my ($self, @features) = @_;
    if (@features)
    {
        push (@{$self->{_features}}, @features);
    }
    return @{$self->{_features}};
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for est2genome from the database and 
                by rerunning blast.[Possibly needs to be moved down to run 
                level?]
    Returns :   nothing
    Args    :   none

=cut

 sub fetch_input {
   my( $self) = @_; 
   print STDERR "Fetching input \n";
   $self->throw("No input id") unless defined($self->input_id);
   my $contigid  = $self->input_id;
   my $contig    = $self->dbobj->get_Contig($contigid);
   my $genseq = $contig->get_repeatmasked_seq();
   my @features = $contig->get_all_SimilarityFeatures;
   $self->{_genseq} = $genseq;
   $self->{_features} = [];
   $self->features(@features);
  
   # process the features for blast
   my @new_features = $self->_process_features(@features);
   
   # set up the runnables
   $self->_prepare_runnables(@new_features);
 }

=head2 _process_features

    Title   :   _process_features
    Usage   :   $self->_process_features(@features)
    Function:   Screens FeaturePairs in @features for vert EST blast hits, validates and
                reblasts sequences that previously had hits, and returns the new hits as 
                FeaturePairs associated with $genseq
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   array of Bio::EnsEMBL::FeaturePair

=cut

sub _process_features {
  my ($self,@features) = @_;
  my $genseq = $self->{_genseq};
  my @mrnafeatures = ();
  my @new_features = ();
  my %idhash = ();
  
  unless (@features)
    {
      print STDERR "Contig has no associated features\n";
      return;
    }
  
  foreach my $f (@features) {
    if (defined($f->analysis)      && defined($f->score) && 
	defined($f->analysis->db)  && $f->analysis->db eq "vert") {
      
      if (!defined($idhash{$f->hseqname})) { 
	push(@mrnafeatures,$f);
	$idhash{$f->hseqname} =1;
      } 
    }
  }
  print STDERR "Number of features pre blast is " . scalar(@mrnafeatures) . "\n";
  
  unless (@mrnafeatures)
    {
      print STDERR ("Contig has no features suitable for Blast\n");
      return;
    }
  
  # run blast   
  my @seq         = $self->get_Sequences(@mrnafeatures);
  my @valid_seq   = $self->validate_sequence(@seq);
  my $blastdb     = $self->make_blast_db(@valid_seq);
  my @newfeatures = $self->run_blast($genseq,$blastdb);
  
  print STDERR "Number of features post blast is " . scalar(@newfeatures) . "\n";
  
  unless (@newfeatures)
    {
      print STDERR ("No new features returned by blast\n");
      return;
    }
  
  return @newfeatures;
}

=head2 _prepare_runnables
    Title   :   _prepare_runnables
    Usage   :   $self->_prepare_runnables(@new_features)
    Function:   Makes a Bio::EnsEMBL::Runnable::ALignFeature for plus strand blast 
                FeaturePairs, and one for minus strand hits as appropriate.
    Returns :   Nothing, but $self->{_runnables} contains the two runnables
    Args    :   array of Bio::EnsEMBL::FeaturePair

=cut

sub _prepare_runnables {
  my ($self, @new_features) = @_;
  my $genseq = $self->{_genseq};
  my @plusfeat;
  my @minusfeat;
   foreach my $f(@new_features) {
     if ($f->hstrand == -1) { push (@minusfeat, $f);  }
     else { push (@plusfeat,$f); }
   }
  
  if( scalar(@plusfeat) ) {
    print STDERR scalar(@plusfeat) . "features on plus strand\n";
    my $prunnable = new Bio::EnsEMBL::Pipeline::Runnable::AlignFeature('-genomic'  => $genseq,
								       '-features' => \@plusfeat);
   
    $self->runnable($prunnable);
  }
  else { print STDERR "no features on plus strand\n"; }

  if( scalar  (@minusfeat) ) {
    print STDERR scalar(@minusfeat) . "features on minus strand\n";
    my $mrunnable = new Bio::EnsEMBL::Pipeline::Runnable::AlignFeature('-genomic'  => $genseq,
								       '-features' => \@minusfeat);
    
    $self->runnable($mrunnable);
  }
  else { print STDERR "no features on minus strand\n"; }
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   None

=cut

sub output {
    my ($self) = @_;
   
    if (!defined($self->{_output})) {
      $self->{_output} = [];
    } 
    return @{$self->{_output}};
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   1
    Args    :   none

=cut

sub write_output {

  my($self) = @_;
  
  my @features = $self->output();
  
  my $db = $self->dbobj();
  my $gene_obj = $db->gene_Obj;
  
  my $EXON_ID_SUBSCRIPT       = "VEGE";
  my $TRANSCRIPT_ID_SUBSCRIPT = "VEGT";
  my $GENE_ID_SUBSCRIPT       = "VEGG";
  my $PROTEIN_ID_SUBSCRIPT    = "VEGP";
  
  my $sth = $db->prepare("lock tables gene write, genetype write, exon write, transcript write, exon_transcript write, translation write,dna read,contig read,clone read,feature read,analysis write");

  $sth->execute;
  
  eval {
    (my $gcount = $gene_obj->get_new_GeneID($GENE_ID_SUBSCRIPT))
      =~ s/$GENE_ID_SUBSCRIPT//;
    (my $tcount = $gene_obj->get_new_TranscriptID($TRANSCRIPT_ID_SUBSCRIPT))
      =~ s/$TRANSCRIPT_ID_SUBSCRIPT//;
    (my $pcount = $gene_obj->get_new_TranslationID($PROTEIN_ID_SUBSCRIPT))
      =~ s/$PROTEIN_ID_SUBSCRIPT//;
    (my $ecount = $gene_obj->get_new_ExonID($EXON_ID_SUBSCRIPT))
      =~ s/$EXON_ID_SUBSCRIPT//;
    
    
    foreach my $gene (@features) {
      $gene->id($GENE_ID_SUBSCRIPT . $gcount);
      $gcount++;
      
      # Convert all exon ids and save in a hash
      my %namehash;
      
      foreach my $ex ($gene->each_unique_Exon) {
	$namehash{$ex->id} = $EXON_ID_SUBSCRIPT.$ecount;
	$ex->id($EXON_ID_SUBSCRIPT.$ecount);
	$ecount++;
      }
      
      foreach my $tran ($gene->each_Transcript) {
	$tran->id($TRANSCRIPT_ID_SUBSCRIPT . $tcount);
	$tran->translation->id($PROTEIN_ID_SUBSCRIPT . $pcount);
	
	my $translation = $tran->translation;
	
	$tcount++;
	$pcount++;
	
	foreach my $ex ($tran->each_Exon) {
	  my @sf = $ex->each_Supporting_Feature;

	  if (defined($namehash{$translation->start_exon_id}) && $namehash{$translation->start_exon_id} ne "") {
	    $translation->start_exon_id($namehash{$translation->start_exon_id});
	  }
	  if (defined($namehash{$translation->end_exon_id}) &&$namehash{$translation->end_exon_id} ne "") {
	    $translation->end_exon_id  ($namehash{$translation->end_exon_id});
	  }
	}
	
      }
      
      $gene_obj->write($gene);
    }
  };

  if ($@) {
    $sth = $db->prepare("unlock tables");
    $sth->execute;
    
    $self->throw("Error writing gene for " . $self->input_id . " [$@]\n");
  } else {
    $sth = $db->prepare("unlock tables");
    $sth->execute;
  }
  return 1;
}

=head2 runnable

    Title   :   runnable
    Usage   :   $self->runnable($arg)
    Function:   Sets a runnable for this RunnableDB
    Returns :   Bio::EnsEMBL::Pipeline::RunnableI
    Args    :   Bio::EnsEMBL::Pipeline::RunnableI

=cut

sub runnable {
  my ($self,$arg) = @_;
  
  if (!defined($self->{_runnables})) {
    $self->{_runnables} = [];
  }
  
  if (defined($arg)) {
    if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
      push(@{$self->{_runnables}},$arg);
    } else {
      $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
  }
  
  return @{$self->{_runnables}};
  
}

=head2 run

    Title   :   run
    Usage   :   $self->run()
 Function:   Runs the est2genome analysis, producing Bio::EnsEMBL::Gene predictions
    Returns :   Nothing, but $self{_output} contains the predicted genes.
    Args    :   None

=cut

sub run {
  my ($self) = @_;
  
  $self->throw("Can't run - no runnable objects") unless defined($self->{_runnables});
  
  foreach my $runnable ($self->runnable) {
    $runnable->minirun;
  }  

  # sort out predicted genes
  $self->_convert_output();
}

=head2 _convert_output

    Title   :   _convert_output
    Usage   :   $self->_convert_output()
    Function:   Converts est2genome output into an array of genes remapped into genomic coordinates
    Returns :   Nothing, but $self->{_output} contains remapped genes
    Args    :   None
=cut

# get merged features into a form where they can be stored in the database.
sub _convert_output {
  my ($self) =@_;
  my $count=1;

  # make an array of genes for each runnable
  foreach my $runnable ($self->runnable) {
    my $contig = $self->dbobj->get_Contig($self->input_id);
    my @features   = $runnable->output;
    my %homols = ();
    my @genes = ();
    
    # sort the hits into a hash keyed by hid
    foreach my $f(@features) {
      push (@{$homols{$f->hseqname()}}, $f);
    }
    
    # merge the features for each hid & make one gene per hid
    foreach my $id (keys %homols) {
      my @tmpf =  @{$homols{$id}};
      my $gene = $self->_make_gene($contig, $count, @tmpf);
      if(defined($gene)) {
	push (@genes, $gene);
	$count++;
      }
    }
    
    if (!defined($self->{_output})) {
      $self->{_output} = [];
    }
    
    push(@{$self->{_output}},@genes);
  }
}

=head2 _make_gene

    Title   :   _make_gene
    Usage   :   $self->_make_gene($contig,$count,@exonfeat)
    Function:   Converts merged exon features into a gene
    Returns :   Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, integer, array of Bio::EnsEMBL::FeaturePair

=cut

sub _make_gene {
  my ($self, $contig, $count, @allfeat) = @_;
  my $excount = 1;
  my $time   = time; 
  chomp($time); 

  print STDERR "number of features to _make_gene: " . scalar(@allfeat) . "\n";

  my $gene   = new Bio::EnsEMBL::Gene;
  $gene->type('ve2g');
  $gene->id($self->input_id . ".exo.$count");
  $gene->created($time);
  $gene->modified($time);
  $gene->version(1);
  
  # make a transcript
  my $tran   = new Bio::EnsEMBL::Transcript;
  $tran->id($self->input_id . ".exo.$count");
  $tran->created($time);
  $tran->modified($time);
  $tran->version(1);
  
  my $transl = new Bio::EnsEMBL::Translation;
  $transl->id($self->input_id . ".exo.$count");
  $transl->version(1);
  
  # add transcript to gene
  $gene->add_Transcript($tran);
  $tran->translation($transl);
  
  my @exons = $self->_make_exons($contig, @allfeat);

  foreach my $exon(@exons) {
    print STDERR "debug debug\n";
    $exon->id($self->input_id . ".exo.$count.$excount");
    $exon->created($time);
    $exon->modified($time);   
    $excount++;
  }

  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
    return;
  } 
  else {
    # assemble exons
    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } 
    else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }
    
    foreach my $exon(@exons){
      $tran->add_Exon($exon);
    }
    
    $transl->start_exon_id($exons[0]->id);
    $transl->end_exon_id  ($exons[$#exons]->id);
    
    if ($exons[0]->strand == 1) {
      $transl->start($exons[0]->start);
      $transl->end  ($exons[$#exons]->end);
    } else {
      $transl->start($exons[0]->end);
      $transl->end  ($exons[$#exons]->start);
    }
    return $gene;
  }
}

=head2 _make_exons

    Title   :   _make_exons
    Usage   :   $self->_make_exons($contig, @features)
    Function:   Produces an array of exons corresponding to the input FeaturePairs
    Returns :   Array of Bio::EnsEMBL::Exon
    Args    :   Bio::EnsEMBL::Virtual::Contig, array of Bio::EnsEMBL::FeaturePair

=cut

sub _make_exons {
  my ($self,$contig,@features)  = @_;
  # should be settable? parameters?
  my $overlap       = 20;
  my $query_gap     = 15;

  print STDERR "number of features to _make_exons: " . scalar(@features) . "\n";

  my $analysis_obj = Bio::EnsEMBL::Analysis->new
    (
     -db              => 'vert',
     -dbversion       => 1,
     -program         => 'est_genome',
     -program_version => '1',
     -gff_source      => 'Vert_est2genome',
     -gff_feature     => 'predicted_exon'
    );
  
  my @exons   = ();  
  my $count = 0;
  @features = sort { $a->start <=> $b->start} @features;

  # make a Exon corresponding to the first feature
  my $exon = $self->_new_exon($count,$features[0]);
  $self->throw("Failed to make exon") unless defined($exon);
  $exon->contig_id($contig->id);
  $exon->attach_seq($contig->primary_seq);
  push(@exons,$exon);
  
  for (my $i=0; $i < $#features; $i++) {
    my $id  = $features[$i]  ->id;
    my $id2 = $features[$i+1]->id;
    
    # First case is if start of next hit is < end of previous
    if ($features[$i]->end > $features[$i+1]->start && 
	($features[$i]->end - $features[$i+1]->start) < $overlap) {
      if ($features[$i]->strand eq "1") {
	$exons[$count]-> end($features[$i+1]->end);
	$exons[$count]->hend($features[$i+1]->hend);
      } else {
	$exons[$count]-> end($features[$i+1]->end);
	$exons[$count]->hend($features[$i+1]->hstart);
      }

      $exons[$count]->add_Supporting_Feature($features[$i+1]);

      # what to do with score now?
      ## How many bases match in this feature?
      #my $bases = ($features[$i+1]->score/100) * ($features[$i+1]->end - $features[$i+1]->start + 1);
      #$newfeatures[$count]->score($newfeatures[$count]->score + $bases);
      
      if ($features[$i+1]->hstart == $features[$i+1]->hend) {
	$features[$i+1]->strand($features[$i]->strand);
      }
      
      # Allow a small gap if < $query_gap, $homol_gap
    } elsif (($features[$i]->end < $features[$i+1]->start) &&
	     abs($features[$i+1]->start - $features[$i]->end) <= $query_gap) {

      $exons[$count]->add_Supporting_Feature($features[$i+1]);     
      if ($features[$i]->strand eq $features[$i]->hstrand || $features[$i]->hstrand eq "1"){
	$exons[$count]->end($features[$i+1]->end);
      }
      else {
	# strand = 1, hstrand = -1
	$exons[$count]->end($features[$i+1]->end);
      }	  
 
      # no place for score on exon ...     
      #      my $bases = ($features[$i+1]->score/100) * ($features[$i+1]->end - $features[$i+1]->start + 1);
      #      $newfeatures[$count]->score($newfeatures[$count]->score + $bases);
        
      if ($features[$i+1]->hstart == $features[$i+1]->hend) {
	$features[$i+1]->strand($features[$i]->strand);
      }
      
    } else {
      # we can't extend the merged homologies so start a
      # new exon
      $count++;
      $i++;
      
      my $newexon = $self->_new_exon($count,$features[$i]);
      $self->throw("Failed to make exon") unless defined($newexon);
      $newexon->contig_id($contig->id);
      $newexon->attach_seq($contig->primary_seq);
      push(@exons,$newexon);

      $i--;
    }
  }
  
# need some jiggery pokery to get a new analysis for this lot. 
  foreach my $ex(@exons) {
    foreach my $fp($ex->each_Supporting_Feature) {
      $fp->analysis($analysis_obj);
      # fix those strands
      if($fp->strand == $fp->hstrand) {
	$fp->strand(1);
	$fp->hstrand(1);
      }
      else {
	$fp->strand(-1);
	$fp->hstrand(-1);
      }
      
      # hmmmm
      $ex->strand($fp->strand);
    }
  }

  return (@exons);
}

=head2 _new_exon

    Title   :   _new_exon
    Usage   :   $self->_new_exon($count,$featurepair)
    Function:   Produces an exon from a Featurepair    
    Returns :   Bio::EnsEMBL::Exon
    Args    :   count, Bio::EnsEMBL::FeaturePair

=cut

sub _new_exon {
  my ($self, $count, $fp) = @_;
  my $exon = new Bio::EnsEMBL::Exon;

  print STDERR "feature passed to _new_exon: \n";
  $self->_print_FeaturePair($fp);

  $exon->version(1);
  $exon->start($fp->start);
  $exon->end  ($fp->end);
  $exon->phase(0);
  $exon->add_Supporting_Feature($fp);
  print STDERR "_new_exon: $exon\n";
  print STDERR "_new_exon start: " . $exon->start . "\n";
  return ($exon);
}

=head2 get_Sequences

    Title   :   get_Sequences
    Usage   :   $self->get_Sequences(@features)
    Function:   Gets a Bio::Seq for each of the hit sequences in @features
    Returns :   Array of Bio::Seq
    Args    :   Array of Bio::EnsEMBL::FeaturePair

=cut
  
sub get_Sequences {
  my ($self,@pairs) = @_;
  
  my @seq;
  
  foreach my $pair (@pairs) {
    my $id = $pair->hseqname;
    if ($pair->analysis->db eq "vert") {
      
      eval {
	my $seq = $self->get_Sequence($id);
	push(@seq,$seq);
      };
      if ($@) {
	$self->warn("Couldn't fetch sequence for $id [$@]");
      } 
    }
  }
  return @seq;
}

=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id);
  Function: Fetches sequence that has ID $id. Tries a variety of methods for 
            fetching sequence. If sequence found, it is cached.
  Returns : Bio::PrimarySeq
  Args    : ID string 

=cut

sub get_Sequence {
  my ($self,$id) = @_;

  if (defined($self->{_seq_cache}{$id})) {
    return $self->{_seq_cache}{$id};
  } 
  
  my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;
  my $seq = $seqfetcher->run_efetch($id);
  
  if(!defined $seq) {
    $seq = $seqfetcher->run_getz($id, 'embl');
  }
  
  if(!defined $seq) {
    $self->throw("Couldn't find sequence for [$id]");
  }

  print (STDERR "Found sequence for $newid [" . $seq->length() . "]\n");
  $self->{_seq_cache}{$id} = $seq;
  
  return $seq;
  
}

=head2 validate_sequence

    Title   :   validate_sequence
    Usage   :   $self->validate_sequence(@seq)
    Function:   Takes an array of Seq or PrimarySeq objects and 
                returns valid ones, removing invalid characters 
                for nucleic acid sequences and rejecting sequences 
                that are not nucleic acid
    Returns :   Array of Bio::Seq
    Args    :   Array of Bio::Seq

=cut

sub validate_sequence {
  my ($self, @seq) = @_;
  my @validated;
  foreach my $seq (@seq)
    {
      print STDERR ("mrna feature $seq is not a Bio::PrimarySeq or Bio::Seq\n") 
	unless ($seq->isa("Bio::PrimarySeq") ||
		$seq->isa("Bio::Seq"));
      my $sequence = $seq->seq;
      if ($sequence !~ /[^acgtn]/i)
        {
	  push (@validated, $seq);
        }
      else 
        {
	  $_ = $sequence;
	  my $len = length ($_);
	  my $invalidCharCount = tr/mrwsykvhdbxMRWSYKVHDBX/n/;
	  #extract invalid characters
	  $sequence =~ s/[ACGTN]//ig;
	  if ($invalidCharCount / $len > 0.05)
            {
	      $self->warn("Ignoring ".$seq->display_id()
			  ." contains more than 5% ($invalidCharCount) "
			  ."odd nucleotide codes ($sequence)\n Type returns "
			  .$seq->moltype().")\n");
            }
	  else
            {
	      $self->warn ("Cleaned up ".$seq->display_id
			   ." for blast : $invalidCharCount invalid chars ($sequence)\n");
	      $seq->seq($_);
	      push (@validated, $seq);
            }
        }
    } 
  return @validated;  
}

=head2 make_blast_db

    Title   :   make_blast_db
    Usage   :   $self->make_blast_db(@seq)
    Function:   Makes a blast formatted databases from the sequences in @seq
    Returns :   Name of file containing blast database
    Args    :   Array of Bio::Seq

=cut

sub make_blast_db {
  my ($self,@seq) = @_;
  
  my $blastfile = $self->get_tmp_file('/tmp/','blast','fa');
  my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
			      -file   => ">$blastfile");
  print STDERR "seq io is " . $seqio . "\n";
  print STDERR "Blast db file is $blastfile\n";
  foreach my $seq (@seq) {
    print STDERR "Writing seq " . $seq->id ."\n";
    $seqio->write_seq($seq);
  }
  
  close($seqio->_filehandle);
  
  my $status = system("pressdb $blastfile");
  print (STDERR "Status from pressdb $status\n");
  
  return $blastfile;
}

=head2 run_blast

    Title   :   run_blast
    Usage   :   $self->run_blast($seq,$db)
    Function:   Runs blastn, querying $seq against $db
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   a Bio::Seq object (query sequence), name of file containing 
                the blast database.

=cut

sub run_blast {
  my ($self,$seq,$db) = @_;
  
  my $blastout = $self->get_tmp_file("/tmp/","blast","tblastn_vert.msptmp");
  my $seqfile  = $self->get_tmp_file("/tmp/","seq","fa");
  
  my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
			      -file   => ">$seqfile");
  
  $seqio->write_seq($seq);
  close($seqio->_filehandle);
  
  my $command  = "blastn $db $seqfile B=500 -hspmax 1000 2> /dev/null |MSPcrunch -d - >  $blastout";
  print STDERR ("Running command $command\n");
  my $status = system($command );
  
  print("Exit status of blast is $status\n");
  open (BLAST, "<$blastout") 
    or $self->throw ("Unable to open Blast output $blastout: $!");    
  if (<BLAST> =~ /BLAST ERROR: FATAL:  Nothing in the database to search!?/)
    {
      print "nothing found\n";
      return;
    }
  
  my $msp = new Bio::EnsEMBL::Analysis::MSPcrunch(-file => $blastout,
						  -type => 'DNA-DNA',
						  -source_tag => 'vert_eg',
						  -contig_id => $self->input_id,
						 );
  unlink $db.".csq";
  unlink $db.".ntb";
  unlink $db.".nhd";
  unlink $blastout;
  unlink $seqfile;
  unlink $db;
  
  my @pairs = $msp->each_Homol;
  
  foreach my $pair (@pairs) {
    $self->_print_FeaturePair($pair);
  }
  return @pairs;
}

=head2 _print_FeaturePair

    Title   :   _print_FeaturePair
    Usage   :   $self->_print_FeaturePair($pair)
    Function:   Prints attributes of a Bio::EnsEMBL::FeaturePair
    Returns :   Nothing
    Args    :   A Bio::EnsEMBL::FeaturePair

=cut

sub _print_FeaturePair {
  my ($self,$pair) = @_;
  
  print STDERR $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . 
               $pair->score . "\t" . $pair->strand . "\t" . $pair->hseqname . "\t" . 
	       $pair->hstart . "\t" . $pair->hend . "\t" . $pair->hstrand . "\n";
}

=head2 get_tmp_file

    Title   :   get_tmp_file
    Usage   :   $self->get_tmp_file($dir,$stub,$ext)
    Function:   Generates a unique file name in directory $dir with combination 
                of $stub, a random number, and $ext
    Returns :   Array of Bio::EnsEMBL::FeaturePair
    Args    :   None

=cut

sub get_tmp_file {
    my ($self,$dir,$stub,$ext) = @_;

    
    if ($dir !~ /\/$/) {
	$dir = $dir . "/";
    }

    # but just ignores the return value!!!
    $self->check_disk_space($dir);

    my $num = int(rand(10000));
    my $file = $dir . $stub . "." . $num . "." . $ext;

# this introduces a race condition; suspect it's the way it's done all over EnsEMBL; could be a problem with big pipeline runs?
    while (-e $file) {
	$num = int(rand(10000));
	$file = $stub . "." . $num . "." . $ext;
    }			
    
    return $file;
}

=head2 check_disk_space

    Title   :   check_disk_space
    Usage   :   $self->check_disk_space($dir, $minimumkb)
    Function:   Checks space in $dir
    Returns :   1 if "enough" space, 0 otherwise
    Args    :   $dir (directory name), $minimumkb (minimumkb needed)

=cut

sub check_disk_space {
    my ($self,$dir,$minimumkb) = @_;

    # what if $minimumkb is unset? check for undef?

    $self->throw("No directory entered") unless defined($dir);

    open(DF,"df -k $dir |");

    my @lines = <DF>;
    $self->throw("Wrong number of lines output from df") unless scalar(@lines) == 2;
    my @f = split(' ',$lines[1]);

    my $kbytes = $f[3];

    if ($minimumkb && $kbytes > $minimumkb) {
	return 1;
    } else {
	return 0;
    }
}

1;
