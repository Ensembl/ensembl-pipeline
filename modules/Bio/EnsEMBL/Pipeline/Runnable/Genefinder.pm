

package Bio::EnsEMBL::Pipeline::Runnable::Genefinder;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis; 
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::TranscriptFactory;
use Bio::EnsEMBL::Root;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    $self->{'_exons'}     = [];    # an array of Bio::Seqfeatures (exons)
    $self->{'_genes'}     = [];    # an array of arrays of SeqFeature 
    $self->{_transcripts} = [];
    $self->{'_tablenamefile'} = undef;
    $self->{'_intronpenalty'} = undef;
    $self->{'_exonpenalty'} = undef;
    $self->{'_query'}     = undef; # location of Bio::Seq object
    $self->{'_genefinder'}   = undef; # location of Genefinder script
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of genscan
    $self->{'_protected'} = [];    # a list of file suffixes protected from deletion
    $self->{'_parameters'} =undef; #location of parameters for genscan
    print "@args\n";
    my($query, $genefinder, $parameters, $tablenamefile, $exonpenalty, $intronpenalty) = 
        $self->_rearrange([qw(CLONE GENEFINDER PARAM TABLENAMEFILE EXONPENALTY INTRONPENALTY)], @args);

   
    $self->query($query);
    $genefinder = 'genefinder'       unless ($genefinder);
    $tablenamefile = '/usr/local/ensembl/data/nemtables/sanger.tablenamefile.cod' unless($tablenamefile);
    $intronpenalty = '/usr/local/ensembl/data/nemtables/intron_penalty.lookup' unless($intronpenalty);
    $exonpenalty = '/usr/local/ensembl/data/nemtables/exon_penalty.lookup' unless($exonpenalty);
    $self->genefinder($self->find_executable($genefinder));
    $self->tablenamefile($tablenamefile);
    $self->intronpenalty($intronpenalty);
    $self->exonpenalty($exonpenalty);
    

    if ($parameters)    
    { $self->parameters($parameters) ; }
    else                
    {$self->parameters(''); }     

    return $self;
}


###################
#get/set methods
###################

sub query {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_query'} = $seq ;
        $self->filename($self->query->id.".$$.seq");
        $self->results($self->filename.".genscan");
    }
    return $self->{'_query'};
}


sub genefinder {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Genefinder not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_genefinder'} = $location ;
    }
    return $self->{'_genefinder'};
}


sub tablenamefile{
  my ($self, $location) = @_;
  
  if ($location)
    {
      $self->{'_tablenamefile'} = $location ;
    }
  return $self->{'_tablenamefile'};
  

}

sub intronpenalty{
  my ($self, $location) = @_;
  
  if ($location)
    {
      $self->{'_intronpenalty'} = $location ;
    }
  return $self->{'_intronpenalty'};
  

}

sub exonpenalty{
  my ($self, $location) = @_;
  
  if ($location)
    {
      $self->{'_exonpenalty'} = $location ;
    }
  return $self->{'_exonpenalty'};
  

}

sub parameters {
  my ($self, $param) = @_;
  if ($param)
    {
      $self->{'_parameters'} = $param;
    }
  return $self->{'_parameters'};
}

sub exons {
  my ($self, $exon) =@_;
  if ($exon)
    {
      $exon->isa("Bio::EnsEMBL::SeqFeature") 
	|| $self->throw("Input isn't a Bio::EnsEMBL::SeqFeature");
      push(@{$self->{'_exons'}}, $exon);
    }
  return @{$self->{'_exons'}};
}

#empties exon array after they have been converted to genes
sub clear_exons {
    my ($self) = @_;
    $self->{'_exons'} = [];
  }

sub genefinder_genes {
    my ($self, $gene) =@_;
    if ($gene)
    {
        $gene->isa("Bio::EnsEMBL::SeqFeature") 
                || $self->throw("Input isn't a Bio::EnsEMBL::SeqFeature");
        push(@{$self->{'_genes'}}, $gene);
        @{$self->{'_genes'}} = sort { $a->seqname <=> $b->seqname } @{$self->{'_genes'}};
    }
    return @{$self->{'_genes'}};
  }


sub add_transcript{
  my ($self, $transcript) =@_;
  if ($transcript)
    {
      $transcript->isa("Bio::EnsEMBL::Transcript") 
	|| $self->throw("Input isn't a Bio::EnsEMBL::Transcript");
      push(@{$self->{'_transcripts'}}, $transcript);
    }
  return @{$self->{'_transcripts'}};
}

sub each_transcript{
  my ($self) = @_;

  if (!defined($self->{'_transcripts'})) {
    $self->{_transcripts} = [];
  }
  return @{$self->{'_transcripts'}};
}

###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs genefinder and creates array of sub-seqfeatures
    Returns :   none
    Args    :   none

=cut

sub run {

    my ($self, $dir) = @_;
    #check seq
    my $seq = $self->query() || $self->throw("Seq required for Genefinder\n");
    #set directory if provided
    if($dir){
      $self->workdir($dir);
    }else{
      $self->workdir('/tmp');
    }
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run genefinder       
    $self->run_genefinder();
    my $size = -s $self->results;
    #print "2 ".$size."\n";
    #parse output and create features
    $self->parse_results();
    $size = -s $self->results;
    #print "3 ".$size."\n";
    $self->deletefiles();
  }


sub run_genefinder {
    my ($self) = @_;
    #print "Running genefinder on ".$self->filename."\n";
    #print STDERR `pwd`;
    my $command = $self->genefinder.' -tablenamefile '.$self->tablenamefile.' -intronPenaltyFile '.$self->intronpenalty.' -exonPenaltyFile '.$self->exonpenalty.' /tmp/'.$self->filename.' > '.$self->results; 
    #print STDERR $command."\n";
    system ($command);
    my $size = -s $self->results;
    #print $size."\n";
    $self->throw($self->results." not created by genefinder\n") unless (-s $self->results);
}



sub parse_results {
  my ($self) = @_;
  #print STDERR "parsing results\n";
  my %phaseMap = ( 0 => 0, 1 => 2, 2 => 1 ) ; # bizarre rule
  
  my ($seq, $prefix, $seqnum, $exon) ;
  my ($strand, $score, $phase, $tsl_offset, $tsl_score, $start, $end, @elts) ;
  
  my $resultsfile = $self->results;
  my $size = -s $self->results;
  #print "size of file = ".$size."\n";
  #print STDERR "opening ".$self->results."\n";
  open(FH, "<".$self->results) || die("couldn't open file $! \n");
  
  while (<FH>) {
    #print "running through file\n";
    #print;
    chomp;
    if (/Sequence: (\S+)/) {
      $seq = $1 ;
      $prefix = $seq ; $prefix =~ s/run20/gf/ ;
      $seqnum = 0 ;
    }
    
    next unless /\.\./ ;
    
    # Typical lines are:
    #
    #    [ 11.52 ]   51958..52077   52126..52399   52450..52676 end
    #    [ 12.58 ]  [TSL: -17  3.2]   37537..37698   37746..37877 end
    # *  [ 29.80 ]   47920..48144   48192..48389   48440..49238   50304..50326 start [TSL: -7  1.9] 
    #    [ 18.14 ]   48..975 U1
    # *U0 [ 9.01 ]   342..466   841..865   971..1055   1108..1176 U2
    #
    # The tricky cases are the ones with U\d, which means the gene is continues beyond the sequence.
    # The digit determines the phase.  This can come at the start or end of the sequence.
    
    ++$seqnum ;
    $exon = 0 ;
    $phase = 0 ;
    
    # extract strand, score and TSL info
    
    if (s/^\*//) { $strand = '-1' ; } else { $strand = '1' ; }
    if (s/\[ -?(\d+\.\d+) \]//) { $score = $1 ; } else { die "can't extract score: $_" ; }
    if (s/\[TSL: (\S+)\s+(\S+)\]//) { $tsl_offset = $1 ; $tsl_score = $2 ; } else { undef $tsl_offset, $tsl_score ; }
    
    # make array of remaining material, reverse if '-' strand, and process

    @elts = split ;
    @elts = reverse @elts if $strand eq '-1' ;
    for (@elts) {
      if (/U(\d)/) {
	$phase = $phaseMap{$1} ;
      } elsif (($start, $end) = /(\d+)\.\.(\d+)/) {
	++$exon ;
	
	#print join ("\t", $seq,"Genefinder","CDS",$start,$end,$score,$strand,$phase,
	#              "Sequence $prefix.$seqnum ; Exon $prefix.$seqnum.$exon\n") ;
	#print STDERR "seqname = $prefix.$seqnum.$exon\n";
	my $exonname = $prefix.".".$seqnum;
	#print STDERR "exonname = ".$exonname."\n";
	$self->create_feature($start, $end, $score, $strand, $phase, $exonname, 'genefinder', '1', 'genefinder', 'prediction');
	# phase update is tricky
	# next exon phase is number of bases remaining from partial codon of previous exon
	
	$phase = (3 - ($end - $start + 4 - $phase) % 3) % 3 ;
	
	#                          \------------------------/
	#                             3 + number of bp used from start of first full codon
	#                          \-----------------------------/
	#                             number of bases started into partial codon
	  
	if ($exon == 1 && defined $tsl_offset) {
	  if ($strand eq '1') {
	    print join ("\t", $seq,"Genefinder","TSL",
			($start+$tsl_offset-1),($start+$tsl_offset),
			$tsl_score,$strand,'.',"Sequence $prefix.$seqnum\n") ;
	  } else {
	    print join ("\t", $seq,"Genefinder","TSL",
			($end-$tsl_offset),($end-$tsl_offset+1),
			$tsl_score,$strand,'.',"Sequence $prefix.$seqnum\n") ;
	  }
	}
      } elsif (!/start/ && !/end/) {
	die "unknown field $_\n" ;
      }
    }
  }
  $self->create_genes;

}



sub add_Transcript {
  my ($self,$transcript) = @_;

  if (defined($transcript)) {
    if (!defined($self->{_transcripts})) {
      $self->{_transcripts} = [];
    }
    push(@{$self->{_transcripts}},$transcript);
  }
}

sub each_Transcript {
  my ($self) = @_;

  if (!defined($self->{_transcripts})) {
    $self->{_transcripts} = [];
  }
  return @{$self->{_transcripts}};
}


sub create_feature {
    my ($self, $start, $end, $score, $strand, $phase, $seqname, $program, $program_version, $source, $primary) = @_;
    #$self->create_gene();
    #print "seqname ".$seqname."\n";
    #create analysis object
    my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => $program,
                            -program_version => $program_version,
                            -gff_source      => $source,
                            -gff_feature     => $primary);

    #create and fill Bio::EnsEMBL::Seqfeature objects   
    my $exon = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $seqname,
                            -start   => $start,
                            -end     => $end,
                            -strand  => $strand,
                            -score   => $score,
                            -phase   => $phase,
                            -analysis => $analysis_obj);  
    $self->exons($exon);
}

#creates groups of exons as subseqfeatures.
#relies on seqname of exons being in genefinder format 3.01, 3.02 etc
sub create_genes {
    my ($self) = @_;
    print "creating genes \n";
    my (%genes, %gene_start, %gene_end, %gene_score, %gene_p,
        %gene_strand, %gene_source, %gene_primary, %gene_analysis);

    my @ordered_exons = sort { $a->start <=> $b->start } $self->exons();

    #sort exons into hash by initial numbers of seqname (genes)
    foreach my $exon (@ordered_exons)
    {
      #print $exon->seqname."\n";
        my ($group_number) = $exon->seqname;
	#print "seqname =  ".$exon->seqname."\n";
	print $group_number."\n";
        #intialise values for new gene
        unless (defined ($genes {$group_number}))
        {
	  print "creating a new trancripts for ".$group_number."\n";
            $genes          {$group_number} = [];
            $gene_start     {$group_number} = $exon->start;
            $gene_end       {$group_number} = $exon->end;
            $gene_score     {$group_number} = 0 ;
            $gene_strand    {$group_number} = $exon->strand;
            $gene_analysis  {$group_number} = $exon->analysis;
        }
        #fill array of exons
	print "adding an exon to ".$group_number."\n";
        push (@{$genes {$group_number}}, $exon);
        #calculate gene boundaries and total score
        $gene_start {$group_number} = $exon->start() 
            if ($exon->start() < $gene_start{$group_number});
        $gene_end   {$group_number} = $exon->end() 
            if ($exon->end() > $gene_end{$group_number});
        $gene_score {$group_number} += $exon->score();
        
    }

    #create Bio::SeqFeature objects (genes) with SubSeqFeatures (exons)
    foreach my $gene_number (keys(%genes))
    {
        my $gene = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname     => $gene_number,
                            -strand      => $gene_strand   {$gene_number},
                            -score       => $gene_score    {$gene_number}
                                            /(scalar @{$genes {$gene_number}}),
                            -start       => $gene_start    {$gene_number},
                            -end         => $gene_end      {$gene_number},
                            -analysis    => $gene_analysis {$gene_number}, )
                    or $self->throw("Couldn't create Bio::EnsEMBL::SeqFeature object");

        foreach my $exon (@{$genes {$gene_number}})
        {
          $gene->add_sub_SeqFeature($exon, '');
        }
        $self->genefinder_genes($gene); #add gene to main object
	my $tran = Bio::EnsEMBL::TranscriptFactory::fset2transcript_with_seq($gene, $self->query);
	#print "have ".$tran."\n";
	$self->add_transcript($tran);
      }
}

##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of SeqFeatures representing predicted genes 
                with exons stored as SubSeqFeatures.
    Returns :   An array of SeqFeatures (genes) containing sub-seqfeatures (exons)
    Args    :   none

=cut

sub output {

  my ($self) = @_;

  my @feat;

  my $analysis = Bio::EnsEMBL::Analysis->new(   -db              => undef,
                                                -db_version      => undef,
                                                -program         => 'genefinder',
                                                -program_version => 1,
                                                -gff_source      => 'genefinder',
                                                -gff_feature     => 'prediction',
                                                -logic_name      => 'genefinder',
                                                );

  
  foreach my $transcript ($self->each_transcript) {
    my @exons = $transcript->get_all_Exons;

    if ($exons[0]->strand == 1) {
      @exons = sort {$a->start <=> $b->start } @exons;
    } else {
      @exons = sort {$b->start <=> $a->start } @exons;
    }
    
    #print "\ntranscript  translates to ".$transcript->translate->seq."\n\n";
    foreach my $exon (@exons) {
      my $f = new Bio::EnsEMBL::SeqFeature(-seqname => $exon->seqname,
                                           -start   => $exon->start,
                                           -end     => $exon->end,
                                           -strand  => $exon->strand,
                                           -phase   => $exon->phase,
                                           -end_phase => $exon->end_phase,
                                           -score   => $exon->score,
                                           -analysis     => $analysis);
      my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $exon->seqname,
                                            -start   => $exon->start,
                                            -end     => $exon->end,
                                            -strand  => $exon->strand,
                                            -phase   => $exon->phase,
                                            -end_phase => $exon->end_phase,
                                            -score   => $exon->score,
                                            -p_value   => $exon->p_value,
                                            -analysis     => $analysis);

#      print STDERR $exon->start . " " . $exon->end . " " . $exon->phase . " " . $exon->strand . "\n";

      $f->analysis($analysis);
      $f2->analysis($analysis);

      my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
                                             -feature2 => $f2);

      push(@feat,$fp);
    }

  }
  return @feat;
}

=head2 output_exons

    Title   :   output_exons
    Usage   :   obj->output_exons()
    Function:   Returns an array of SeqFeatures corresponding to exons
    Returns :   An array of SeqFeatures corresponding to exons
    Args    :   none

=cut

sub output_exons {
    my ($self) = @_;
    my @exons;

    foreach my $gene ($self->genefinder_genes)
    {
        push (@exons, $gene->sub_SeqFeature);
    }
    print STDERR "No exons predicted\n" unless (@exons);
    @exons = sort { $a->seqname <=> $b->seqname } @exons;
    return @exons;
}

=head2 output_singlefeature

    Title   :   output_singlefeature
    Usage   :   obj->output_singlefeature()
    Function:   Returns a single SeqFeature with exons as sub-SeqFeatures
    Returns :   A single SeqFeature with exons as sub-SeqFeatures
    Args    :   none

=cut

sub output_singlefeature {
    my ($self) = @_;
    my ($start, $end, $score, $analysis, $primary, $source, $p_value);

    my (@genes) = $self->genes();
    print STDERR "No exons predicted\n" unless (@genes);
    #calculate boundaries and aggregate values
    foreach my $gene (@genes)
    {
        $start      =  $gene->start()  if (!defined($start) || $gene->start() < $start);
        $end        =  $gene->end()    if (!defined($end)   || $gene->end()   > $end);
        $score      += $gene->score();
        $p_value    += $gene->p_value();
        $analysis   =  $gene->analysis();
        $primary    =  $gene->primary_tag();
        $source     =  $gene->source_tag();
    }
    $score  = $score/scalar(@genes); #average score

    my $single = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname        => 'genefinder',
                            -strand         => 1,
                            -score          => $score,
                            -start          => $start,
                            -end            => $end,
                            -p_value        => $p_value,
                            -analysis       => $analysis )
                    or $self->throw("Couldn't create Bio::EnsEMBL::SeqFeature object");

    foreach my $gene (@genes)
    {
        foreach my $exon ($gene->sub_SeqFeature())
        {
            $single->add_sub_SeqFeature($exon, '');
        }    
    }
    return $single;
}

1;
