#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Snap

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq'; 
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::Snap object
  my $Snap = Bio::EnsEMBL::Pipeline::Runnable::Snap->new (-CLONE => $seq);
  $Snap->workdir($workdir);
  $Snap->run();
  my @genes = $snap->output();
  my @exons = $snap->output_exons();
  my $seqfeature = $snap->output_singlefeature();

=head1 DESCRIPTION

Snap takes a Bio::Seq (or Bio::PrimarySeq) object and runs Snap on it. The
resulting output is parsed to produce a set of Bio::SeqFeatures. 

snap is a gene predictor written by Ian Korf (ik1@sanger.ac.uk) part the Zoe software library.

=head2 Methods:

=over 4

=item new($seq_obj)

=item    Snap($path_to_Snap)

=item    workdir($directory_name)

=item    run()

=item    output()

=item    output_exons()

=item    output_genes()

=back

=head1 SEE ALSO

=over 4

=item B<Bio::EnsEMBL::Pipeline::RunnableI>

=item B<Bio::EnsEMBL::Pipeline::RunnableDB::Snap> 

=back

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Snap;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis; 
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::TranscriptFactory;
use Bio::EnsEMBL::Root;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Snap->new (-CLONE => $seq);
    Function:   Initialises Snap object
    Returns :   a Snap Object
    Args    :   A Bio::Seq object 
                (Snap location and matrix file location optional)

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    

    $self->{'_exons'}     = [];    # an array of Bio::Seqfeatures (exons)
    $self->{'_genes'}     = [];    # an array of arrays of SeqFeatures
    $self->{'_peptides'}  = [];    # Snap predicted peptide (used for phase)
    $self->{'_snap'}   = undef; # location of Snap script
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of Snap
    $self->{'_protected'} = [];    # a list of file suffixes protected from deletion
    $self->{'_parameters'} =undef; #location of parameters for Snap
    $self->{'_hmmfile'} = undef;

    my($query, $snap, $parameters, $matrix) = 
        $self->_rearrange([qw(QUERY SNAP PARAM HMMFILE)], @args);


    $self->query($query);
#    $snap = 'snap'       unless ($snap);

    if ($snap) {
	$self->snap($snap);
    }
    else {
	$self->snap($self->find_executable($snap)) unless ($snap);
    }

    print STDERR "HMMFILE1: $matrix\n";
    
    $self->hmmfile($matrix);


    print STDERR "HMMFILE2: ".$self->hmmfile."\n";

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
#Snap takes all of the fasta header, want to make sure that the fasta header contains only one field

#	my $id = $seq->id;
#	my $s = $seq->seq;

#	my $new_seq = Bio::Seq->new(-seq=>$s,-id=>$id);

        $self->{'_query'} = $seq ;
        $self->filename($seq->id.".$$.seq");
        $self->results($self->filename.".snap");
	$self->protfile($self->filename.".snapprot");
	
	$self->file($self->results);
    }
    return $self->{'_query'};
}


=head2 snap

    Title   :   snap
    Usage   :   $obj->snap('');
    Function:   Get/set method for the location of snap
    Args    :   File path (optional)

=cut

sub snap {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Snap not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_snap'} = $location ;
    }
    return $self->{'_snap'};
}

=head2 hmmfile

    Title   :   hmmfile
    Usage   :   $obj->hmmfile('fly');
    Function:   Get/set method for the location of genscan hmmfile
    Args    :   File path (optional)

=cut

sub hmmfile {
    my ($self, $location) = @_;
    if ($location) {
	$self->{'_hmmfile'} = $location ;
    }

    return $self->{'_hmmfile'};
}


=head2 parameters

    Title   :   parameters
    Usage   :   $obj->parameters('parameters');
    Function:   Get/set method for the location of genscan parameters
    Args    :   File path (optional)

=cut

sub parameters {
    my ($self, $param) = @_;
    if ($param)
    {
        $self->{'_parameters'} = $param;
    }
    return $self->{'_parameters'};
}

=head2 protfile

 Title   : protfile
 Usage   : $obj->protfile('protein_file')
 Function: Get/Set method for the location of snap output file containing the protein files
 Returns : value of protfile
 Args    : newvalue (optional)


=cut

sub protfile{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'protfile'} = $value;
    }
    return $obj->{'protfile'};

}


sub exons {
    my ($self, $exon) =@_;
    if ($exon)
    {
      #print "adding ".$exon->seqname."\n";
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

sub snap_genes {
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

sub snap_peptides {
    my ($self, $peptide) = @_;
    push (@{$self->{'_peptides'}}, $peptide) if ($peptide);
    return @{$self->{'_peptides'}};
}

###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs snap and creates array of sub-seqfeatures
    Returns :   none
    Args    :   none

=cut

sub run {

    my ($self) = @_;
    #check seq
    my $seq = $self->query() || $self->throw("Seq required for Snap\n");
    #set directory if provided
    $self->workdir('/tmp') unless $self->workdir();
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run snap       
    $self->run_snap();
    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
    unlink ($self->protfile) or $self->throw ("Couldn't delete ".$self->protfile.":$!");
    return 1;
}

sub run_snap {
    my ($self) = @_;
#    print STDERR "Running snap on ".$self->snap.' '.$self->hmmfile.' '.$self->filename.' -gff -aa '.$self->protfile.'> '.$self->results."\n";
    
    system ($self->snap.' '.$self->hmmfile.' '.$self->filename.' -gff -aa '.$self->protfile.'> '.$self->results);
   
   $self->throw($self->results." (GFF output) not created by Snap\n") unless (-e $self->results);
    $self->throw($self->protfile." (protein file) not created by Snap\n") unless (-e $self->results);
}

=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses Snap output to give a set of seqfeatures
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
                NOTE: snap can not assign phases to exons from the output
                file unless the sequence is supplied as a Bio::Seq object.
    Returns :   none
    Args    :   optional filename

=cut

sub parse_results {
    my ($self) = @_;

    my $in  = Bio::SeqIO->new ( '-format' => 'Fasta' , -file => $self->protfile);

    while (my $seq = $in->next_seq) {
	$self->snap_peptides($seq->seq); #final peptide 
    }
    
    
    my %exon_type = ('Esngl', 'Single Exon',
                     'Einit', 'Initial Exon',
                     'Exon', 'Internal Exon',
                     'Eterm', 'Terminal Exon');

    my $filehandle;
     
    if (ref ($self->results) !~ /GLOB/)
    {
        open (SNAP, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*SNAP;
    }
    else
    {
        $filehandle = $self->results;
    }

    my @element;
    
    my $excount;
    my $transcount = 0;

    my $current_trans;

    #The big parsing loop - parses exons and predicted peptides
    while (<$filehandle>) 
    {

	
	print STDERR "$_";
	my %feature; 
	
	@element = split;
	
	$self->throw("Unable to parse Snap ouput (".scalar(@element).") Line: $_\n") unless (scalar(@element) == 9); 
	
	
	if ($current_trans ne $element[8]) {
	    $excount = 0;
	    $transcount++;
	    $current_trans = $element[8];
	}
	
	$excount++;

	my $name = $transcount.".".$excount;


	if ($element[6] eq '+')
	{
	    $feature {'start'} = $element[3];
	    $feature {'end'} = $element[4];
	    $feature {'strand'} = 1;
	}
	elsif ($element[6] eq '-')
	{
	    $feature {'start'} = $element[3];
	    $feature {'end'} = $element[4];
	    $feature {'strand'} = -1;
	}
	else {
	    die ("output file wronlgy formated\n");
	}
	
	$feature{'name'} = $name;
	$feature {'score'} = $element[5];
	$feature {'type'} = $exon_type{$element[2]};
	$feature {'program'} = 'snap';
	$feature {'program_version'} = '1.0';
	$feature {'primary'} = 'prediction';
	$feature {'source'} = 'snap';
	
	$self->create_feature(\%feature);
	
#	if (($element[2] eq "Eterm")||( $element[2] eq "Esngl")) {
	#    $transcount++;
	#    $excount = 0;
	#}
    }

    $self->create_genes();    
    print STDERR "there are ".scalar($self->snap_genes)." genes\n";

    $self->calculate_and_set_phases_new();

    $self->clear_exons(); #free up unecessary storage 
}


sub calculate_and_set_phases_new {
    my ($self) = @_;
    my $min_gene_length = 25;

    my @genes       = $self->snap_genes();
    my @peptides    = $self->snap_peptides();

    $self->throw("Mismatch in number of genes (".scalar(@genes).
                 ") and peptides ("             .scalar(@peptides).
                 ") parsed from file") unless (scalar(@genes) == scalar (@peptides));

    my $i = 0;
    my $count = 1;
  GENE: while ($i < scalar(@genes)) {

      print STDERR "Gene " . $genes[$i]->seqname . "\n";

      my $peptide = $peptides[$i];

      if (length $peptide < $min_gene_length) {
          print STDERR "peptide $i too short (", length($peptide), "bp)\n";
          $i++;
          next GENE;
      }
      
      if ($peptide =~ /X{5,}?/) {
	  print STDERR "peptide $i ($peptide) has too much low complexity\n";
	  $i++;
	  next GENE;
      }

      my @exons   = $genes[$i]->sub_SeqFeature();
      print STDERR "Exons are $#exons\n";
      my @newtran = Bio::EnsEMBL::TranscriptFactory::fset2transcript_3frame($genes[$i],$self->query);

      print STDERR "\nPeptide is " . $peptides[$i] . "\n";

      my $translation_found = 0;

      # remove any initial X's from the peptide
      $peptide =~ s/^x//i;
        
      # remove any terminal X's from the peptide
      $peptide =~ s/x$//i;

      foreach my $tran (@newtran) {

        my $temp_tran = $tran->translate->seq;

       print STDERR "Translation is " . $temp_tran . "\n";

        # clean the translated sequence

        #snap translated partial genes correctly whilst exon translation begin with M
        #$temp_tran =~ s/^M//i; #remove initial M from exon

        # remove any initial X's from the translation
        $temp_tran =~ s/^x//i;
        
        # remove any terminal X's from the translation
        $temp_tran =~ s/x$//i;

       

        my $x = 0;
        while (($x = index($temp_tran, 'X', $x)) != -1) {
          #print STDERR "Found an 'X' at ", $i + 1, "\n";
          substr($peptide, $x, 1) = 'X'
           if length($peptide) >= length($temp_tran);
          $x++;
        }

        $x = 0;
        while (($x = index($peptide, 'X', $x)) != -1) {
          #print STDERR "Found an 'X' at ", $i + 1, "\n";
          substr($temp_tran, $x, 1) = 'X'
           if length($temp_tran) >= length($peptide);
          $x++;
        }
    
        #print "\nafter\nsnap: $peptide\nensembl: $ensembl\n";#
	
        # The next line is a hack to fix the inconsistency generated by the translation when  
        # 2 bases unambiguously determined an aminoacid. 
        # The translation might have an extra aminoacid that must be removed! mc2
 
	if (length($temp_tran)==length($peptide)+1){
	    chop($temp_tran);
	}
        if (index($peptide ,$temp_tran) >= 0) {
#         print STDERR $tran->temporary_id . " " . $tran->translate->seq . "\n";

          $translation_found = 1;
          foreach my $exon (@{$tran->get_all_Exons()}) {
           #print "exon ".$exon->seqname." ".$exon->start . " " . $exon->end . " " . $exon->phase . " " . $exon->end_phase . " " .$exon->strand . "\n";
          }
          $tran->temporary_id($self->query->id . "." . $count);
          $count++;
          #print "translation = ".$temp_tran."\n";
          if($temp_tran =~/\*/){
            print "translation ".$tran->temporary_id." = ".$temp_tran."\n";
            print "transcript contains stop codons!\n";
          } else {
           #print "adding translation ".$tran->temporary_id." \n";
            $self->add_Transcript($tran);
          }
          $i++;
          next GENE;
        }
      }
      
      unless ($translation_found) {
        $self->throw("[Snap.pm] Unable to match Snap peptide ".$peptides[$i].
                     " in a translation\n");
      }

#      print "\n";
      $i++;
    }
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

sub get_all_Transcripts {
  my ( $self ) = @_;

  if (!defined($self->{_transcripts})) {
    $self->{_transcripts} = [];
  }
  return $self->{_transcripts};
}


sub each_Transcript {
  my ($self) = @_;

  my $transcripts = $self->get_all_Transcripts();
  $self->warn( "each_Transcript deprecated, use get_all_Transcripts()" );
  return @{$self->get_all_Transcripts()};
}


sub create_feature {
    my ($self, $feat) = @_;

    #create analysis object
    my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => $feat->{'program'},
                            -program_version => $feat->{'program_version'},
                            -gff_source      => $feat->{'source'},
                            -gff_feature     => $feat->{'primary'});

    #create and fill Bio::EnsEMBL::Seqfeature objects 

    my $exon = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $feat->{'name'},
                            -start   => $feat->{'start'},
                            -end     => $feat->{'end'},
                            -strand  => $feat->{'strand'},
                            -score   => $feat->{'score'},
                            -frame   => $feat->{'frame'},
                            -p_value => $feat->{'p'},
                            -analysis => $analysis_obj);

    $self->exons($exon);

}

#creates groups of exons as subseqfeatures.
#relies on seqname of exons
sub create_genes {
    my ($self) = @_;
    #print "creating genes \n";
    my (%genes, %gene_start, %gene_end, %gene_score, %gene_p,
        %gene_strand, %gene_source, %gene_primary, %gene_analysis);

    my @ordered_exons = sort { $a->seqname <=> $b->seqname } $self->exons();
#    foreach my $exon(@ordered_exons){
#      print STDERR " EXON ".$exon."\n";
#    }
    #sort exons into hash by initial numbers of seqname (genes)
    foreach my $exon (@ordered_exons)
      {
	#print "seqname ".$exon->seqname."\n"; 
        my ($group_number) = $exon->seqname =~ /(\d*)./;
	
	print STDERR "group number".$group_number."\n";
        #print "seqname =  ".$exon->seqname."\n";
        #intialise values for new gene
        unless (defined ($genes {$group_number}))
        {
            $genes          {$group_number} = [];
            $gene_start     {$group_number} = $exon->start;
            $gene_end       {$group_number} = $exon->end;
            $gene_score     {$group_number} = 0 ;
            $gene_strand    {$group_number} = $exon->strand;
            $gene_analysis  {$group_number} = $exon->analysis;
            $gene_p         {$group_number} = 0 ;
        }
        #fill array of exons
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
	#print "gene number ".$gene_number."\n";
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
        $self->snap_genes($gene); #add gene to main object
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
    my @pred;

    my $analysis = Bio::EnsEMBL::Analysis->new(
        -db              => undef,
        -db_version      => undef,
        -program         => 'snap',
        -program_version => 1,
        -gff_source      => 'snap',
        -gff_feature     => 'prediction',
        -logic_name      => 'snap',
    );

    foreach my $transcript (@{$self->get_all_Transcripts}) {

	
        my $exons = $transcript->get_all_Exons();
        

	my @exons;

        if ($exons->[0]->strand == 1) {
            @exons = sort {$a->start <=> $b->start } @{$exons};
        } else {
            @exons = sort {$b->start <=> $a->start } @{$exons};
        }

        push @pred, Bio::EnsEMBL::PredictionTranscript->new(@exons);
    }
    return @pred;
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

    foreach my $gene ($self->snap_genes)
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
                        (   -seqname        => 'snap',
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
