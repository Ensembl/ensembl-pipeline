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

Bio::EnsEMBL::Pipeline::Runnable::Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new(
                                             -genomic => $genseq,
                                             -est     => $estseq 
                                             );
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new()

=head1 DESCRIPTION

Object to store the details of an est2genome run.
Stores the est2genome matches as an array of Bio::EnsEMBL::FeaturePair

=head2 Methods:

 new,
 genomic_sequence,
 est_sequence,
 run,
 output.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
#compile time check for executable
use Bio::EnsEMBL::Analysis::Programs qw(est_genome); 
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::Root::RootI;

#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-GENOMIC       => $genomicseq,
			   -EST           => $estseq,
                           -E2G           => $executable,
                           -ARGS          => $args);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::Runnable::Est2Genome object
    Returns :   A Bio::EnsEMBL::Pipeline::Runnable::Est2Genome object
    Args    :   -genomic:    Bio::PrimarySeqI object (genomic sequence)
                -est:        Bio::PrimarySeqI object (est sequence), 
                -e2g:        Path to Est2Genome executable
                -args:       Arguments when running est2genome 
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
           
    $self->{'_fplist'}     = [];    # create key to an array of feature pairs
    $self->{'_clone'}      = undef; # location of Bio::Seq object
    $self->{'_est_genome'} = undef; # location of est2genome
    $self->{'_workdir'}    = undef; # location of temp directory
    $self->{'_filename'}   = undef; # file to store Bio::Seq object
    $self->{'_estfilename'}= undef; # file to store EST Bio::Seq object
    $self->{'_results'}    = undef; # file to store results of analysis
    $self->{'_protected'}  = [];    # a list of files protected from deletion
    $self->{'_arguments'}  = undef; # arguments for est2genome
    
    my( $genomic, $est, $est_genome, $arguments ) = 
        $self->_rearrange([qw(GENOMIC EST E2G ARGS)], @args);
    print STDERR "Genomic $genomic\n"; 
    $self->genomic_sequence($genomic) if $genomic; #create & fill key to Bio::Seq
    $self->est_sequence($est) if $est; #create & fill key to Bio::Seq
    if ($est_genome) 
    {   $self->est_genome($est_genome) ;}
    else
    {   
        eval 
        { $self->est_genome($self->locate_executable('est_genome')); };
        if ($@)
        { $self->est_genome('/usr/local/pubseq/bin/est_genome'); }
    }
    if ($arguments) 
    {   $self->arguments($arguments) ;}
    else
    { $self->arguments(' -reverse ') ;      }
    
    return $self;
}

#################
# get/set methods 
#################

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
        $self->filename($value->id.".$$.seq");
        $self->results($self->filename.".est_genome.out");
    }
    return $self->{'_genomic_sequence'};
}

=head2 est_sequence

    Title   :   est_sequence
    Usage   :   $self->est_sequence($seq)
    Function:   Get/set method for est sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub est_sequence {
    my( $self, $value ) = @_;
    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_est_sequence'} = $value;
        $self->estfilename($value->id.".$$.est.seq");
    }
    return $self->{'_est_sequence'};
}

sub estfilename {
    my ($self, $estfilename) = @_;
    $self->{'_estfilename'} = $estfilename if ($estfilename);
    return $self->{'_estfilename'};
}

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

=head2 workdir

    Title   :   workdir
    Usage   :   $obj->wordir('~humpub/temp');
    Function:   Get/set method for the location of a directory to contain temp files
    Args    :   File path (optional)

=cut

=head2 arguments

    Title   :   arguments
    Usage   :   $obj->arguments(' -reverse ');
    Function:   Get/set method for est_genome arguments
    Args    :   File path (optional)

=cut

sub arguments {
    my ($self, $args) = @_;
    if ($args)
    {
        $self->{'_arguments'} = $args ;
    }
    return $self->{'_arguments'};
}

###########
# Analysis methods
##########

=head2 run

  Title   : run
  Usage   : $self->run()
            or
            $self->run("genomic.seq", "est.seq")
  Function: Runs est2genome and stores results as FeaturePairs
  Returns : TRUE on success, FALSE on failure.
  Args    : Temporary filenames for genomic and est sequences

=cut

# this is a MESS
sub run {
    my ($self, @args) = @_;
    
    # some constant strings
    my $source_tag  = "est2genome";
    my $dirname     = "/tmp/";
    
    #flag for est strand orientation
    my $estOrientation; 
    
    #check inputs
    my $genomicseq = $self->genomic_sequence ||
        $self->throw("Genomic sequence not provided");
    my $estseq = $self->est_sequence ||
        $self->throw("EST sequence not provided");
    
    #extract filenames from args and check/create files and directory
    my ($genname, $estname) = $self->_rearrange(['genomic', 'est'], @args);
    my ($genfile, $estfile) = $self->_createfiles($genname, $estname, $dirname);
        
    #use appropriate Bio::Seq method to write fasta format files
    {
        my $genOutput = Bio::SeqIO->new(-file => ">$genfile" , '-format' => 'Fasta')
                    or $self->throw("Can't create new Bio::SeqIO from $genfile '$' : $!");
        my $estOutput = Bio::SeqIO->new(-file => ">$estfile" , '-format' => 'Fasta')
                    or $self->throw("Can't create new Bio::SeqIO from $estfile '$' : $!");

        #fill inputs
        $genOutput->write_seq($self->{'_genomic_sequence'}); 
        $estOutput->write_seq($self->{'_est_sequence'});

    }
        
    #The -reverse switch ensures correct numbering on EST seq in either orientation
    my $est_genome_command = "est_genome  -reverse -genome $genfile -est $estfile |";
    print STDERR "running for " . $estseq->display_id . "\n";
    eval {
      print (STDERR "Running command $est_genome_command\n");
      open (ESTGENOME, $est_genome_command) 
	or $self->throw("Can't open pipe from '$est_genome_command' : $!");
      
      #Use the first line to get gene orientation
      my $firstline = <ESTGENOME>;
      print STDERR "firstline: \t$firstline\n";
      # put the gene on the minus strand iff splice sites imply reversed gene
      if ($firstline =~ /REVERSE/) { print STDERR "***reversed gene***\n"; $estOrientation = -1; }
      else {$estOrientation = 1}
     
      #read output
      while (<ESTGENOME>) {
	if ($_ =~ /Segmentation fault/) {
	  $self->warn("Segmentation fault from est_genome\n");
	  close (ESTGENOME) or $self->warn("problem running est_genome: $!\n");
	  return(0);
	}
	elsif ($_ =~ /^(Segment|Exon|Span)/) {
	  
	  #split on whitespace
	  my @elements = split;
	  
	  #extract values from output line
	  my $f1score  = $elements[2];
	  my $f1start  = $elements[3];
	  my $f1end    = $elements[4];
	  my $f1id     = $elements[5];
	  my ($f2start, $f2end);
	  my $f2id     = $elements[8];
	  my $f1source = $source_tag;
	  my $f2source = $source_tag;
	  #my $f1strand = 1;
	  my $f1strand = $estOrientation; # otherwise this is going to get lost later on ....
	  #my $f2strand = $estOrientation;
	  my $f2strand = 1;
	  my $f1primary = $elements[0];
	  my $f2primary = $f1primary;
	  #ensure start is always less than end
	  if ($elements[6] < $elements[7])
	    {
	      $f2start =  $elements[6]; 
	      $f2end = $elements[7];
	    }
	  else
	    {
	      $f2start =  $elements[7]; 
	      $f2end = $elements[6];
	    }              
	  #create array of featurepairs              
	  $self->_createfeatures ($f1score, $f1start, $f1end, $f1id, 
				  $f2start, $f2end, $f2id, $f1source, 
				  $f2source, $f1strand, $f2strand, 
				  $f1primary, $f2primary);
        }    

      }
      if(!close(ESTGENOME)){
	$self->warn("problem running est_genome: $!\n");
	return(0);
      }

      $self->convert_output;

    };

    $self->_deletefiles($genfile, $estfile);
    if ($@) {
        $self->throw("Error running est_genome [$@]\n");
    } else {
        return 1;
    }
}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
  my ($self) = @_;
  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }
  return @{$self->{'_output'}};
}

sub convert_output {
  my ($self) = @_;
  my @genes;
  my @exons;
  my @supp_feat;

  # split the different features up
  foreach my $f(@{$self->{'_fplist'}}){
    if ($f->primary_tag eq 'Span'){
      push(@genes, $f);
    }
    elsif($f->primary_tag eq 'Exon'){
      push(@exons, $f);
    }
    elsif($f->primary_tag eq 'Segment'){
      push(@supp_feat, $f);
    }
  }
  
  # now reassemble them
  # add exons to genes
  foreach my $ex(@exons){
    my $added = 0;
    
    foreach my $g(@genes){
      if($ex->start >= $g->start  && $ex->end <= $g->end
	 && $ex->strand == $g->strand && !$added){
	$g->feature1->add_sub_SeqFeature($ex,'');
	$added = 1;
      }
    }
    $self->warn("Exon $ex could not be added to a gene ...\n") unless $added;     
  }

  # add supporting features to exons
  foreach my $sf(@supp_feat){
    my $added = 0;
    
    foreach my $e(@exons){
      if($sf->start >= $e->start  && $sf->end <= $e->end
	 && $sf->strand == $e->strand && !$added){
	$e->feature1->add_sub_SeqFeature($sf,'');
	$added = 1;
      }
    }
    $self->warn("Feature $sf could not be added to an exon ...\n") unless $added;     
  }
  
  push(@{$self->{'_output'}},@genes);

}

sub _createfeatures {
    my ($self, $f1score, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id,
        $f1source, $f2source, $f1strand, $f2strand, $f1primary, $f2primary) = @_;
    

    my $analysis_obj    = new Bio::EnsEMBL::Analysis
                                (-db              => "none",
                                 -db_version      => "none",
                                 -program         => "est_genome",
                                 -program_version => "none",
                                 -gff_source      => $f1source,
                                 -gff_feature     => $f1primary,);
    #create features
    my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start      =>   $f1start,
                                              -end         =>   $f1end,
                                              -seqname     =>   $f1id,
                                              -strand      =>   $f1strand,
                                              -score       =>   $f1score,
					      -percent_id  =>   $f1score, 
                                              -analysis    =>   $analysis_obj );
     
   my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $f2start,
                                              -end         =>   $f2end,
                                              -seqname     =>   $f2id,
                                              -strand      =>   $f2strand,
                                              -score       =>   $f1score,
					      -percent_id  =>   $f1score, 
                                              -analysis    =>   $analysis_obj );

    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
                                             -feature2 => $feat2) ;

        push(@{$self->{'_fplist'}}, $fp);
}

sub _createfiles {
    my ($self, $genfile, $estfile, $dirname)= @_;
    
    #check for diskspace
    my $spacelimit = 0.01; # 0.01Gb or about 10 MB
    #my $dir ="./";
    my $dir =$dirname;
    unless ($self->_diskspace($dir, $spacelimit)) 
    {
        $self->throw("Not enough disk space ($spacelimit Gb required)");
    }
            
    #if names not provided create unique names based on process ID    
    $genfile = $self->_getname("genfile") unless ($genfile);
    $estfile = $self->_getname("estfile") unless ($estfile);    
    $genfile = $dirname . $genfile; 
    $estfile = $dirname . $estfile; 
    # Should check we can write to this directory 
    $self->throw("No directory $dirname") unless -e $dirname;

    #chdir ($dirname) or $self->throw ("Cannot change to directory '$dirname' ($?)"); 
    return ($genfile, $estfile);
}
    
sub _getname {
    my ($self, $typename) = @_;
    return  $typename."_".$$.".fn"; 
}

sub _diskspace {
    my ($self, $dir, $limit) =@_;
    my $block_size; #could be used where block size != 512 ?
    my $Gb = 1024 ** 3;
    
    open DF, "df $dir |" or $self->throw ("Can't open 'du' pipe");
    while (<DF>) 
    {
        if ($block_size) 
        {
            my @L = split;
            my $space_in_Gb = $L[3] * 512 / $Gb;
            return 0 if ($space_in_Gb < $limit);
            return 1;
        } 
        else 
        {
            ($block_size) = /(\d+).+blocks/i
                || $self->throw ("Can't determine block size from:\n$_");
        }
    }
    close DF || $self->throw("Error from 'df' : $!");
}

sub _deletefiles {
    my ($self, @files) = @_;
    
    my $unlinked = unlink(@files);

    if ($unlinked == @files) {
        return 1;
    } else {
        my @fails = grep -e, @files;
#        $self->throw("Failed to remove @fails : $!\n");
        $self->warn("Failed to remove @fails : $!\n");
    }
}

sub est_genome {
   my ($self,$arg) = @_;

  if (defined($arg)) {
     $self->{'_est_genome'} = $arg;
  }
  return $self->{'_est_genome'};
}

1;
