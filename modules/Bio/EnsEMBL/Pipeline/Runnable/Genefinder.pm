#
# Cared for by esnembl-dev  <ensembl-dev@ebi.ac.uk>
#
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Genefinder

=head1 SYNOPSIS

    my $genefinder = Bio::EnsEMBL::Pipeline::Runnable::Genefinder->new(-query => $seq) #should be a Bio::PrimarySeq obj

  will need to alter the hard coded paths to config files from genefinder for this to work properly
 
=head1 DESCRIPTION

  This package is based on Bio::EnsEMBL::Pipeline::Runnable::Genscan


=head2 Methods:

=over 4

=item new($seq_obj)

=item    genefinder($path_to_Genefinder)

=item    workdir($directory_name)

=item    run()

=item    output()

=back

=head1 SEE ALSO

=over 4

=item B<Bio::EnsEMBL::Pipeline::RunnableI>

=item B<Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder> 

=back

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

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
#    print STDERR "ARGS: @args\n";
    my($query, $genefinder, $parameters, $tablenamefile, $exonpenalty, $intronpenalty) = 
        $self->_rearrange([qw(QUERY GENEFINDER PARAM TABLENAMEFILE EXONPENALTY INTRONPENALTY)], @args);

   
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


=head2 accessor methods

  Arg [1]   : scalar/object
  Function  : get/set value to argument passed in
  Returntype: set value
 Exceptions: query will throw if not passed a Bio::PrimarySeq and genefinder will throw if location of executable doesn't exist'
  Caller    : $self
  Example   : $self->query($seq);

=cut



sub query {
    my ($self, $seq) = @_;
  
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input ".$seq." isn't a Bio::Seq or Bio::PrimarySeq");
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


=head2 genefinder_genes

  Arg [1]   : Bio::EnsEMBL::SeqFeature
  Function  : add arg to the array of genes
  Returntype: array of genes
  Exceptions: if arg isn't a Bio::EnsEMBL::SeqFeature'
  Caller    : $self
  Example   : $self->genefinder_genes($gene);

=cut


sub genefinder_genes {
    my ($self, $gene) =@_;
    if ($gene)
    {
        $gene->isa("Bio::EnsEMBL::SeqFeature") 
                || $self->throw("Input isn't a Bio::EnsEMBL::SeqFeature");
        push(@{$self->{'_genes'}}, $gene);
        @{$self->{'_genes'}} = sort { $a->seqname cmp $b->seqname } @{$self->{'_genes'}};
    }
    return @{$self->{'_genes'}};
  }


=head2 add_transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : adds transcript object to array
  Returntype: array of transcripts
 Exceptions: throws if arg isnt a Bio::EnsEMBL::Transcript
  Caller    : $self
  Example   : $self->add_transcript($transcript);

=cut


sub add_Transcript{
  my ($self, $transcript) =@_;
  if ($transcript)
    {
      $transcript->isa("Bio::EnsEMBL::Transcript") 
	|| $self->throw("Input isn't a Bio::EnsEMBL::Transcript");
      push(@{$self->{'_transcripts'}}, $transcript);
    }
  return @{$self->{'_transcripts'}};
}


=head2 each_Transcript

  Arg [1]   : none
  Function  : returns array of transcripts
 Returntype: array of transcripts
 Exceptions: throws if arg isnt a Bio::EnsEMBL::Transcript
  Caller    : $self
  Example   : $self->add_Transcript($transcript);

=cut

sub each_Transcript{
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

  Arg [1]   : working directory
  Function  : run the genefinder analysis
  Returntype: none
  Exceptions: throws if no  sequence object is found
  Caller    : runnableDB normally
  Example   : $runnable->run (see Bio::EnsEMBL::Pipeline::RunnableDB::Genefinder)

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


=head2 run_genefinder

  Arg [1]   : none
  Function  : constructs the command line and runs genefinder
  Returntype: non
  Exceptions: throws if no result file is produces
  Caller    : $self
  Example   : $self->run_genefinder

=cut



sub run_genefinder {
    my ($self) = @_;
    #print "Running genefinder on ".$self->filename."\n";
    #print STDERR `pwd`;
    my $command = $self->genefinder.' -tablenamefile '.$self->tablenamefile.' -intronPenaltyFile '.$self->intronpenalty.' -exonPenaltyFile '.$self->exonpenalty.' /tmp/'.$self->filename.' > '.$self->results; 
    print STDERR $command."\n";
    system ($command);
    my $size = -s $self->results;
    #print $size."\n";
    $self->throw($self->results." not created by genefinder\n") unless (-s $self->results);
}



=head2 parse_results

  Arg [1]   : none
  Function  : parse genefinders evil output
  Returntype: none
  Exceptions: throws if couldn't open outputfile
  Caller    : $self'
  Example   : $self->parse_results

=cut



sub parse_results{
  my ($self) = @_;
  
  # NOTE genefinder output is evil
  # *  [ 7.08 ]   1972..2061   2368..2553	3397..3617   10379..10703 start 
  
  #   [ 9.56 ]   13592..14005 end
  
  #   [ 15.99 ]  [TSL: -1	3.0]   14096..14145   15781..15921   15989..16093   16188..16292   16340..16493	  16761..17021 end
  
  #   [ 43.18 ]  [TSL: -4	3.1]   18068..18139   18193..18392   18438..18607   19152..19350   20217..20481	  20884..20968	 21508..21795	22647..22978 end

  #*  [ 26.51 ]   23612..24130 U0  not this is an incomplete gene in this case the completion is of the end of the sequence 
  #but it can just be interupted by another gene in output file this is the 3' end of a reverse strand gene the U0 indicates the exon start in phase 0
  my $resultsfile = $self->results;
  open(FH, "<".$self->results) || die("couldn't open file ".$self->results." $! \n");
  my $prefix;
  my @forward_lines;
  my @reverse_lines;
  while(<FH>){ # this loop sorts the lines in to foward and reverse genes
    #print;
    chomp;
    if (/Sequence: (\S+)/) {
      $prefix = $1; 
    }
    if($_ =~ /\.\./){
      if($_ =~ /\*/){
	push(@reverse_lines, $_);
      }else{
	push(@forward_lines, $_);
      }
      
    }  
  }
  my $gene_count = 0;
  my $phase = 0;
  my $strand = 1;
  #print STDERR "forward strand\n";
  foreach my $line(@forward_lines){
    #print "parsing gene ".$gene_count."\n";
    #print STDERR "line ".$line,"\n";
    $line =~ s/\[TSL: \S+\s+\S+\]//; #TSL is striped iout, not sure what it is
    $line =~ s/end//; #end is stripped off
    $line =~ s/\[ -?(\d+\.\d+) \]//; #score is parsed out
    my $score = $1;
    #print STDERR "have line ".$line."\n"; 
    my @values = split /\s+/, $line; #remaining line is split into pairs of exon coords
    #print STDERR "have values @values\n";
    if(!$values[0]){
      my $first = shift @values;  
    }
    if($values[0] =~ /U(\d)/){
      $phase = $1;
      my $first = shift @values; 
    }
    #print STDERR "have values @values\n";
    my $count = 0;
    my $exonname = $prefix.".".$gene_count;
    foreach my $coord(@values){ #for each pair of values a feature is created
      #print STDERR "start phase of this gene is ".$phase."\n";
      if($coord =~ /U\d/){
	next;
      }
      #print STDERR "have coord ".$coord."\n";
      my ($start, $end) = split /\.\./, $coord;
      #print STDERR "have start ".$start," end ".$end."\n";
      my $end_phase = ($phase + ($end-$start) + 1)%3;
      $self->create_feature($start, $end, $score, $strand, $phase, $end_phase, $exonname);
      $phase =  $end_phase;
      my $count++;
      if($values[$count] =~ /U(\d)/){
	if($phase != $1){
	  $self->warn(" phase ".$phase." and continuation phase ".$1." aren't the same may be issues with translation\n");
	}
      }
    }
    $phase = 0;
    $gene_count++;
  }
  #print STDERR "reverse strand\n";
  $phase = 0;
  $strand = -1;
  foreach my $line(@reverse_lines){ # this is pretty much as for forward strand genes but the line of coordinate pairs has to be reversed before processing
    #print "parsing gene ".$gene_count."\n";
    #print "phase = ".$phase."\n";
    #print STDERR "line ".$line,"\n";
    $line =~ s/\[TSL: \S+\s+\S+\]//;
    $line =~ s/start//;
    $line =~ s/\[ -?(\d+\.\d+) \]//;
    my $score = $1;
    if(!$line =~ /\*/){
      $self->throw("this gene $line isn't on the reverse strand why is it here");
    }else{
      $line =~ s/\*//;
    }
    my @elements = split /\s+/, $line;
    $line = '';
    foreach my $v (reverse(@elements)){
      $line .= $v." ";
    } 
    #print STDERR $line."\n";
    my @values = split /\s+/, $line;
    #print STDERR "@values\n";
    if(!$values[0]){
      my $first = shift @values;  
    }
    if($values[0] =~ /U(\d)/){
     # print "value is ".$values[0]."\n";
      $phase = $1;
      my $phase_variable = shift @values; 
     # print STDERR "phase has been set to ".$phase."\n";
    }
    #print STDERR "@values\n";
    my $count = 0;
    my $exonname = $prefix.".".$gene_count;
    #print "phase = ".$phase."\n";
    foreach my $coord(@values){
      #print STDERR "start phase of this gene is ".$phase."\n";
      if($coord =~ /U\d/){
	next;
      }
     # print STDERR "phase is ".$phase."\n";
#      print STDERR "have coord ".$coord."\n";
      my ($start, $end) = split /\.\./, $coord;
#      print STDERR "have start ".$start," end ".$end."\n";
      my $end_phase = ($phase + ($end-$start) + 1)%3;
      $self->create_feature($start, $end, $score, $strand, $phase, $end_phase, $exonname,);
      $phase =  $end_phase;
      my $count++;
      if($values[$count] =~ /U(\d)/){
	if($phase != $1){
	  $self->warn(" phase ".$phase." and continuation phase ".$1." aren't the same may be issues with translation\n");
	}
      }
    }
    $phase = 0;
    $gene_count++;
  }
  $self->create_genes;
}



=head2 create_feature

  Arg [1]   : int
  Arg [2]   : int
  Arg [3]   : float
  Arg [4]   : int (strand should only be 1 or -1)
  Arg [5]   : int
  Arg [6]   : int
  Arg [7]   : int 
  Arg [8]   : string
  Function  : create a seqfeature from passed in values and add to array of exons
  Returntype: none
  Exceptions: none
  Caller    : $self
  Example   : (line 452)

=cut




sub create_feature {
    my ($self, $start, $end, $score, $strand, $phase, $end_phase, $seqname) = @_;
    #$self->create_gene();
    #print STDERR " Creating exon seqfeature seqname ".$seqname." start ".$start. " end ".$end." strand ".$strand." phase ".$phase." end phase ".$end_phase." score ".$score."\n" if($strand == -1);
  #  $self->throw("create feature");
    #create and fill Bio::EnsEMBL::Seqfeature objects   
    my $exon = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $seqname,
                            -start   => $start,
                            -end     => $end,
                            -strand  => $strand,
                            -score   => $score,
                            -phase   => $phase,
                            -end_phase => $end_phase);  
    $self->exons($exon);
}


=head2 create_genes

  Arg [1]   : none
  Function  : groups exons on seqname and turns into transcripts and genes
  Returntype: none
  Exceptions: throws if problem creating seqfeature
  Caller    : $self
  Example   : $self->create_genes

=cut


#creates groups of exons as subseqfeatures.
#relies on seqname of exons being in genefinder format 3.01, 3.02 etc
sub create_genes {
    my ($self) = @_;
    #print "creating genes \n";
    my (%genes, %gene_start, %gene_end, %gene_score, %gene_p,
        %gene_strand, %gene_source, %gene_primary, %gene_analysis);

    my @ordered_exons = sort { $a->start <=> $b->start } $self->exons();

    #sort exons into hash by initial numbers of seqname (genes)
    foreach my $exon (@ordered_exons)
    {
      #print $exon->seqname."\n";
        my ($group_number) = $exon->seqname;
	#print "seqname =  ".$exon->seqname."\n";
	#print $group_number."\n";
        #intialise values for new gene
        unless (defined ($genes {$group_number}))
        {
	  #print "creating a new trancripts for ".$group_number."\n";
            $genes          {$group_number} = [];
            $gene_start     {$group_number} = $exon->start;
            $gene_end       {$group_number} = $exon->end;
            $gene_score     {$group_number} = 0 ;
            $gene_strand    {$group_number} = $exon->strand;
            $gene_analysis  {$group_number} = $exon->analysis;
        }
        #fill array of exons
	#print "adding an exon to ".$group_number."\n";
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
	$gene->contig($self->query);
        foreach my $exon (@{$genes {$gene_number}})
        {
          $gene->add_sub_SeqFeature($exon, '');
        }
        $self->genefinder_genes($gene); #add gene to main object
	#print STDERR "contig ".$gene->contig."\n";
	$self->query->id($gene_number);
	my $tran = Bio::EnsEMBL::TranscriptFactory::fset2transcript_with_seq($gene, $self->query);
	#print "have ".$tran."\n";
	$self->add_Transcript($tran);
      }
}

##############
# input/output methods
#############


=head2 output

  Arg [1]   : none
  Function  : turns exons into feature pairs
  Returntype: array of feature pairs
  Exceptions: none
  Caller    : 
  Example   : 

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

  
  foreach my $transcript ($self->each_Transcript) {
    my @exons = @{$transcript->get_all_Exons};

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


1;
