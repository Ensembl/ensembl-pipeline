package Bio::EnsEMBL::Pipeline::Runnable::Fgenesh;

#
#
# Cared for by Laura Clarke  <lec@sanger.ac.uk>
#
# Copyright Laura Clarke
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Fgenesh

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $seqfile = '/nfs/disk65/mq2/temp/bA151E14.seq'; 
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $seqfile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::Fgenesh object
  my $fgenesh = Bio::EnsEMBL::Pipeline::Runnable::Fgenesh->new (-QUERY => $seq);
  $fgenesh->workdir($workdir);
  $fgenesh->run();
  my @genes = $fgenesh->output();
  my @exons = $fgenesh->output_exons();
  my $seqfeature = $fgenesh->output_singlefeature();

=head1 DESCRIPTION

This package is based on the genscan runnable.
Fgenesh takes a Bio::Seq (or Bio::PrimarySeq) object and runs Fgenesh on it. The
resulting output is parsed to produce a set of Bio::SeqFeatures. 

=head2 Methods:

=over 4

=item new($seq_obj)

=item    Fgenesh($path_to_Fgenesh)

=item    workdir($directory_name)

=item    run()

=item    output()

=back

=head1 CONTACT

ensembl-dev.ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 


=cut

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis; 
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::TranscriptFactory;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::PredictionExon;
use Bio::Seq;
use Bio::EnsEMBL::Root;
use Bio::SeqIO;



@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Fgenesh->new (-QUERY => $seq);
    Function:   Initialises Fgeneshobject
    Returns :   a FgeneshObject
    Args    :   A Bio::Seq object 
                (Fgeneshlocation and matrix file location optional)

=cut

sub new{
     my ($class,@args) = @_;
     my $self = $class->SUPER::new(@args);  

     $self->{'_genes'} = []; #an array of arrays of Bio::Seqfeatures
     $self->{'_exons'} = []; #an array of seqfeatures
     $self->{'_transcripts'} = []; #an array of Bio::Ensembl::Transcript
     $self->{'_proteins'} = []; #fgenesh predicted proteins
     $self->{'_query'} = undef; #location of Bio::Seq object
     $self->{'_workdir'} = undef; #location of temp dir
     $self->{'_fgenesh'} = undef; #location of fgenesh binary
     $self->{'_filename'} = undef; #file to store Bio::Seq object
     $self->{'_results'} = undef; #file to store results of fgenesh
     $self->{'_protected'} = [];    # a list of file suffixes protected from deletion
     $self->{'_parameters'} = undef; #location of parameters for fgenesh
     $self->{'_matirx'} = undef; #location of matrix used by fgenesh
     #print "@args\n";
     my($query, $fgenesh, $parameters, $matrix) = 
	  $self->_rearrange([qw(QUERY FGENESH PARAM MATRIX)], @args);

      $self->query($query) if ($query);


     $fgenesh = 'fgenesh' unless ($fgenesh);
     $matrix  = 'hum.dat' unless ($matrix);

     
     $self->fgenesh($self->find_executable($fgenesh));
     $self->matrix ($self->find_file($matrix));


     if ($parameters){ 
       $self->parameters($parameters); 
     }else{
       $self->parameters(''); 
     }

     return $self;

 }


###################
#get/set methods
###################

=head2 query

    Title   :   query
    Usage   :    $Fgenesh->query($seq);
    Function:   sets the sequence the fgenesh object will run on
  and checks it is a Bio::Seq
    Returns :   a seq
    Args    :   A Bio::Seq object 
                

=cut

sub query {
    my ($self, $seq) = @_;
    if ($seq)
    {
      #print "seq = ".$seq."\n";
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_query'} = $seq ;
        $self->filename($self->query->id.".$$.seq");
        $self->results($self->filename.".fgenesh");
    }
    return $self->{'_query'};
}

=head2 fgenesh

    Title   :   fgenesh
    Usage   :   $obj->fgenesh('/nfs/disk100/humpub/OSFbin/fgenesh');
    Function:   Get/set method for the location of fgenesh
    Args    :   File path (optional)

=cut


sub fgenesh {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Fgenesh not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_fgenesh'} = $location ;
    }
    return $self->{'_fgenesh'};
}

=head2 matrix

    Title   :   matrix
    Usage   :   $obj->matrix('/nfs/disk100/humpub/OSFbin/HumanIso.smat');
    Function:   Get/set method for the location of fgenesh matrix
    Args    :   File path (optional)

=cut

sub matrix {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Fgenesh matrix not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_matrix'} = $location ;
    }
    return $self->{'_matrix'};
}

=head2 parameters

    Title   :   parameters
    Usage   :   $obj->parameters('parameters');
    Function:   Get/set method for the location of fgenesh parameters
    Args    :   File path (optional)

=cut

sub parameters {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->{'_parameters'} = $location ;
    }
    return $self->{'_parameters'};
}

=head2 exons

  Title : exons
  Usage : $obj->exons($seqfeature);
  Function : adds to the objects array of exons and returns the array of exons
  Args   : Seqfeature that is an exon
  

=cut

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
=head2 clear_exons

  Title : exons
  Usage : $obj->clear_exons();
  Function : emptys the exon array.
  Args   : 
  

=cut

sub clear_exons {
    my ($self) = @_;
    $self->{'_exons'} = [];
}

=head2 add_Fgenesh_Gene

  Title : add_Fgenesh_Gene
  Usage : $obj->add_Fgenesh_Gene($seqfeature);
  Function : adds to the objects array of genes and returns the array of genes
  Args   : Seqfeature that is a gene
  

=cut

sub add_Fgenesh_Gene {
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

=head2 each_Fgenesh_Gene

  Title : each_Fgenesh_Gene
  Usage : $obj->each_Fgenesh_Gene();
  Function : returns an array of fgenesh genes
  Args   : 
  

=cut


sub each_Fgenesh_Gene {
    my ($self, $gene) =@_;
   
    if (!defined($self->{'_genes'})) {
	$self->{'_genes'} = [];
    }
    
    return @{$self->{'_genes'}};
}
=head2 add_Fgenesh_Protein

  Title : add_Fgenesh_Protein
  Usage : $obj->add_Fgenesh_Protein($seqfeature);
  Function : adds to the objects array of proteins and returns the array of proteins
  Args   : a protein
  

=cut
sub add_Fgenesh_Protein {
    my ($self, $protein) =@_;
    if ($protein)
    {
        push(@{$self->{'_proteins'}}, $protein);
    }
    return @{$self->{'_proteins'}};
}
=head2 each_Fgenesh_Protein

  Title : each_Fgenesh_Protein
  Usage : $obj->each_Fgenesh_Protein();
  Function : returns an array of fgenesh Proteins
  Args   : 
  

=cut

sub each_Fgenesh_Protein {
    my ($self, $protien) =@_;
   
    if (!defined($self->{'_proteins'})) {
	$self->{'_proteins'} = [];
    }
    
    return @{$self->{'_proteins'}};
}

=head2 add_Fgenesh_Transcript

  Title : add_Fgenesh_Transcript
  Usage : $obj->add_Fgenesh_Transcript($seqfeature);
  Function : adds to the objects array of transcripts and returns the array of transcripts
  Args   : a Bio::Ensembl::Transcript
  

=cut

sub add_Fgenesh_Transcript{
    my ($self, $transcript) =@_;
    if ($transcript)
    {
        $transcript->isa("Bio::EnsEMBL::Transcript") 
                || $self->throw("Input isn't a Bio::EnsEMBL::Transcript");
        push(@{$self->{'_transcripts'}}, $transcript);
    }
    return @{$self->{'_transcripts'}};
}
=head2 each_Fgenesh_Transcript

  Title : each_Fgenesh_Transcript
  Usage : $obj->each_Fgenesh_Transcript();
  Function : returns an array of fgenesh Transcripts
  Args   : 
  

=cut
sub each_Fgenesh_Transcript {
  my ($self) = @_;

  if (!defined($self->{'_transcripts'})) {
    $self->{_transcripts} = [];
  }
  return @{$self->{'_transcripts'}};
}

sub get_all_Transcripts {
  my ( $self ) = @_;

  if (!defined($self->{_transcripts})) {
    $self->{_transcripts} = [];
  }
  return $self->{_transcripts};
}


##################
#Analysis methods#
##################

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs fgenesh and creates array of sub-seqfeatures
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir, $filename) = @_;
    #print "check query ".$self->query."\n";
    my $seq = $self->query() || $self->throw("Seq required for Fgenesh\n");
    #print "set directory if provided\n";
    if($dir){
      $self->workdir($dir);
    }
    $self->workdir('/tmp') unless $self->workdir();
    $self->checkdir();
    #write sequence to file
    #print "have checked directory writing file next\n";
    #print "filename to be set to ".$filename."\n";
    if($filename){
      $self->filename($filename);
    }else{
      $self->writefile();
    } 
    print STDERR "about to run Fgenesh\n";
#run fgenesh       
    $self->run_fgenesh();
    print STDERR "have run fgenesh\n";
    #parse output and create features
    $self->parse_results();
    print STDERR "have parsed output\n";
    $self->deletefiles();

    1;
}

=head2 run_fgenesh

    Title   :  run
    Usage   :   $obj->run_fgenesh()
    Function:   Runs fgenesh
    Returns :   none
    Args    :   none

=cut

sub run_fgenesh {
    my ($self) = @_;
    #print STDERR "Running fgenesh on ".$self->filename."\n";
    my $command = $self->fgenesh.' '.$self->matrix.' '.$self->filename.' > '.$self->results;
    print STDERR $command."\n";
    system ($command);
    $self->throw($self->results." not created by fgenesh\n") unless (-e $self->results);
    #print "leaving run_fgenesh\n";
}

=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses fgenesh output to give a set of seqfeatures
                parse_results can accept filenames, filehandles or pipes (\*STDIN)
                NOTE: fgenesh can not assign phases to exons from the output
                file unless the sequence is supplied as a Bio::Seq object.
    Returns :   none
    Args    :   optional filename

=cut

sub parse_results {
    my ($self) = @_;
    #print "parsing data\n";
    my %exon_type = ('CDSo', 'Single Exon',
                     'CDSf', 'Initial Exon',
                     'CDSi', 'Internal Exon',
                     'CDSl', 'Terminal Exon');

    my $filehandle;
    if (ref ($self->results) !~ /GLOB/)
    {
        open (FGENESH, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*FGENESH;
    }
    else
    {
        $filehandle = $self->results;
    }
    #print "have opened results file\n";
    if (<$filehandle> =~ m|no reliable predictions|i )
    {
        print STDERR "No genes predicted\n";
        return;
    }

    #The big parsing loop - parses exons and predicted peptides
    while (<$filehandle>)
    {
        # Last line before predictions contains nothing but spaces and dashes
	
	 
	
	    #print "looking at ".$_."\n";
        
	my $flag = 0;
	if (/^\s*-[-\s]+$/)  
        { 
	    #print "have found data\n";
	    GENES : while (<$filehandle>) 
	    {
		my @lines;
		
		until (/^$/)
		{
		    #print "searching for exons\n";
		    #print "line = " . $_."\n";
		    if (/CDSl|CDSi|CDSf|CDSo/i)
		    {
			#print "parsing data from ". $_ ."\n";
			my @element = split;
			push @lines, \@element; 
			$self->throw("Unable to parse fgenesh ouput (".scalar(@element).") Line: $_\n") 
			    unless (scalar(@element) == 11); 	    
			
		    } elsif (/Predicted protein/i)
		    {
                        #print "Doing PP elsif\n";
			$flag = 1;
			last GENES ;
		    }
                    #print "Getting next line\n";
                    $_ = <$filehandle>;    
		}
		
		#print "this gene has ". scalar @lines." exons\n";
		if($lines[0]->[1] eq '+')
		{
		    @lines = sort {$a->[3] <=> $b->[3]} @lines;
		}
		elsif($lines[0]->[1] eq '-')
		{
		    @lines = reverse sort {$a->[3] <=> $b->[3]} @lines;
		}
		my $exon_num=1;
		foreach my $line(@lines)
		{
		    my %feature;
		    #print "parsing data\n";
		    my $number = $line->[0]+($exon_num/100);
		    $feature {'name'} = $number;
		    if($line->[1] eq '+')
		    {
			$feature {'start'} = $line->[3];
			$feature {'end'} = $line->[5];
			$feature {'strand'} = 1;
			$feature {'phase'} = (3-($line->[7]-$line->[3]))% 3;
		    }
		    elsif($line->[1] eq '-')
		    {
			$feature {'start'} = $line->[3];
			$feature {'end'} = $line->[5];
			$feature {'strand'} = -1;
			$feature {'phase'} = (3-($line->[5]-$line->[9]))% 3;
		    }
		    
		    $feature {'score'} = $line->[6];
		    $feature {'type'} = $exon_type{$line->[2]};
		    $feature {'program'} = 'Fgenesh';
		    $feature {'program_version'} = '1.0';
		    $feature {'primary'} = 'prediction';
		    $feature {'source'} = 'fgenesh';
		    #print "creating an exon". $feature {'name'}."\n";
		    $self->create_feature(\%feature);
		    $exon_num++;
		}
		if($flag==1)
		{
		    last;
		}
	    } 
	    my $protein_string = undef;
	    my $reading_protein =0;
	    while (<$filehandle>)
	    {

#		if (/\A[\W]\Z/)
#		{
#		    $self->add_Fgenesh_Protein($protein_string);
#		    $protein_string = undef;
#		    $reading_protein = 0;
#		}

		if (/^>/)
		{
		    $reading_protein = 1;
                    if (defined $protein_string) {
                        $self->add_Fgenesh_Protein($protein_string);
                    }
		}
		elsif ($reading_protein)
		{
		    $_ =~ s/\n//g;
		    $protein_string .= $_;
		}
	    }
            if (defined $protein_string) {
	        $self->add_Fgenesh_Protein($protein_string); #final peptide
            } 
	}
    }

    $self->create_genes();
    #$self->test_phases();
    $self->clear_exons(); #free up unecessary storage 
    
}

=head2 create_genes
  
 Title: create_genes
 Usage: $obj->create_genes;
 Function:makes genes from the exons produced by running fgenesh
 Returns:
 Args:
   
=cut

sub create_genes {
    my ($self) = @_;
    my (%genes, %gene_start, %gene_end, %gene_score,
        %gene_strand, %gene_source, %gene_primary, %gene_analysis);

    my @ordered_exons = sort { $a->seqname cmp $b->seqname } $self->exons();

    #sort exons into hash by initial numbers of seqname (genes)
    foreach my $exon (@ordered_exons)
      {
	#print $exon->seqname."\n";
        my ($group_number) = ($exon->seqname =~ /(\d+)\./);
	#print $group_number."\n";
        #intialise values for new gene
        unless (defined ($genes {$group_number}))
        {
            $genes          {$group_number} = [];#
            $gene_start     {$group_number} = $exon->start;
            $gene_end       {$group_number} = $exon->end;
            $gene_score     {$group_number} = 0 ;
            $gene_strand    {$group_number} = $exon->strand;
            $gene_source    {$group_number} = $exon->source_tag ;
            $gene_primary   {$group_number} = $exon->primary_tag;
            $gene_analysis  {$group_number} = $exon->analysis;
            
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
	#print "start = ".$gene_start{$gene_number}." gene number = ".$gene_number." \n";
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
        $self->add_Fgenesh_Gene($gene); #add gene to main object
	#print STDERR "seq = ".$self->query."\n";
	my $tran = Bio::EnsEMBL::TranscriptFactory::fset2transcript_with_seq($gene, $self->query);
	#print "have ".$tran."\n";
	$self->add_Fgenesh_Transcript($tran);
    }
    #print "created the genes\n";
}

=head2

 Title: create_feature
 Usage: $obj->create_feature(\%features)
 Function: turns a hash of feature information into a seqfeatures and               adds it to exons
 Returns: Nothing
 Args: Hash containing info about exon;

=cut  
    
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
                            -phase   => $feat->{'phase'},
                            -analysis => $analysis_obj);  
    $self->exons($exon);
    #print "created the exon\n";
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
    print STDERR "Getting output\n";

    foreach my $transcript (@{$self->get_all_Transcripts}) {

        my $exons = $transcript->get_all_Exons();
        print STDERR "Have ".@$exons." exons\n";
        my @exons;

        if ($exons->[0]->strand == 1) {
            @exons = sort {$a->start <=> $b->start } @{$exons};
        } else {
            @exons = sort {$b->start <=> $a->start } @{$exons};
        }
        foreach my $e(@exons){
          bless($e, 'Bio::EnsEMBL::PredictionExon');
        }
        push @pred, Bio::EnsEMBL::PredictionTranscript->new(-exons => \@exons);
    }
    return @pred;
}


1;



















