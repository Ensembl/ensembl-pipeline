#!/usr/local/bin/perl -w

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
# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
#compile time check for executable
use Bio::EnsEMBL::Analysis::Programs qw(est_genome); 
use Bio::PrimarySeq;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    $self->{'_clone'}  = undef;        #location of Bio::Seq object
    $self->{'_est_genome'} = undef;    #location of est2genome
    $self->{'_workdir'}   = undef;     #location of temp directory
    $self->{'_filename'}  =undef;      #file to store Bio::Seq object
    $self->{'_estfilename'} = undef    #file to store EST Bio::Seq object
    $self->{'_results'}   =undef;      #file to store results of analysis
    $self->{'_protected'} =[];         #a list of files protected from deletion
    $self->{'_arguments'} =undef;      #arguments for est2genome
    
    my( $genomic, $est, $est_genome, $arguments ) = 
        $self->_rearrange(['GENOMIC','EST', 'E2G', 'ARGS'], @args);
       
    $self->genomic_sequence($genomic) if $genomic; #create & fill key to Bio::Seq
    $self->est_sequence($est) if $est; #create & fill key to Bio::Seq
    if ($est_genome) 
    {   $self->est_genome($est_genome) ;}
    else
    {   
        eval 
        { $self->est_genome( $self->locate_executable ('est_genome')); };
        if ($@)
        { $self->est_genome('/usr/local/pubseq/bin/est_genome'); }
    }
    if ($arguments)  {
       $self->arguments($arguments) ;
    } else {
       $self->arguments(' -reverse ') ;      
    }    
    return $self; # success - we hope!
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
        $value->isa("Bio::PrimarySeq") || $self->throw("Input isn't a Bio::PrimarySeq");
        $self->{'_genomic_sequence'} = $value;
        $self->filename($value->id .".$$.seq");
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
        $value->isa("Bio::PrimarySeq") || $self->throw("Input isn't a Bio::PrimarySeq");
        $self->{'_est_sequence'} = $value;
        $self->estfilename($value->id.".$$.est.seq");
    }
    return $self->{'_est_sequence'};
}

sub estfilename {
    my ($self, $estfilename) = @_;
    $self->{_estfilename} = $estfilename if ($estfilename);
    return $self->{_estfilename};
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

=head2 est_genome

    Title   :   est_genome
    Usage   :   $obj->est_genome('~humpub/scripts/est_genomeHum');
    Function:   Get/set method for the location of est_genomeHum script
    Args    :   File path (optional)

=cut

sub est_genome {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("est_genome not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_est_genome} = $location ;
    }
    return $self->{_est_genome};
}

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
        $self->{_arguments} = $args ;
    }
    return $self->{_arguments};
}

###########
# Analysis methods
##########

## UPDATING: need to do some careful rearranging here to avoid breaking module 

=head2 run

  Title   : run
  Usage   : $self->run()
            or
            $self->run("genomic.seq", "est.seq", 'tmp')
  Function: Runs est2genome and stores results as FeaturePairs
  Returns : none
  Args    : Temporary filenames for genomic and est sequences and working directory

=cut

sub run {
    my ($self, @args) = @_;
    my ($genname, $estname, $dir) = $self->_rearrange(['genomic', 'est', 'dir'], @args);
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Genscan\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequences to file - pass method names required for access to Bio::Seq's
    $self->writefile('genomic_sequence', 'filename');
    $self->writefile('est_sequence', 'estfilename');
    #run genscan       
    $self->run_est_genome();
    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
}
 
sub run_est_genome {
    my ($self) = @_;        
    #run est_genome 
    #The -reverse switch ensures correct numbering on EST seq in either orientation

    print (STDERR "Running est_genome\n");
    system ($self->est_genome." ".$self->arguments
                      ." -genome ".$self->filename." -est "
                      .$self->estfilename. " > ".$self->results)
	or $self->throw("Failed running from est_genome : $!");
}

sub parse_results {
      my ($self) = @_;
      #Use the first line to get EST orientation
      open (ESTGENOME, '<'.$self->results) 
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
      my $firstline = <ESTGENOME>;
      
      my $estOrientation;
      if ($firstline =~ /reverse/i) {$estOrientation = -1;}
      else {$estOrientation = 1}
      
      #read output
      while (<ESTGENOME>) {

	  if ($_ =~ /^Segment/) {

	  print STDERR "$_";
	  #split on whitespace
	  my @elements = split;
      
      my (%feat1, %feat2);
	  #extract values from output line
	  $feat1 {score}  = $elements[2];
      $feat2 {score}  = $feat1 {score};
	  $feat1 {start}  = $elements[3];
	  $feat1 {end}    = $elements[4];
	  $feat1 {name}   = $elements[5];
	  $feat2 {name}   = $elements[8];
	  $feat1 {source} = 'est_genome';
	  $feat2 {source} = 'est_genome';
	  $feat1 {strand} = 1;
	  $feat2 {strand} = $estOrientation;
	  $feat1 {primary} = 'similarity';
	  $feat2 {primary} = 'similarity';
	  #ensure start is always less than end
	  if ($elements[6] < $elements[7])
	      {
	        $feat2 {start} =  $elements[6];
	        $feat2 {end} = $elements[7];
	      }
	  else
	      {
	        $feat2 {start} =  $elements[7];
	        $feat2 {end} = $elements[6];
	      }
      $feat2 {db} = undef;
      $feat2 {db_version} = undef;
      $feat2 {program} = 'est_genome';
      $feat2 {p_version} = 'unknown';
	  #create array of featurepairs              
	  $self->createfeaturepair (\%feat1, \%feat2);
        }    
      }
      close(ESTGENOME);
}

##############
# input/output methods
#############

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_fplist'}};
}
<<<<<<< Est2Genome.pm

sub _createfeatures {
    my ($self, $f1score, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id,
        $f1source, $f2source, $f1strand, $f2strand, $f1primary, $f2primary) = @_;
    
    #create analysis object
    my $analysis_obj    = new Bio::EnsEMBL::Analysis
                                (-db              => undef,
                                 -db_version      => undef,
                                 -program         => "est_genome",
                                 -program_version => "unknown",
                                 -gff_source      => $f1source,
                                 -gff_feature     => $f1primary,);
    

    #create features
    my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start =>   $f1start,
                                              -end =>      $f1end,
                                              -seqname =>  $f1id,
                                              -strand =>   $f1strand,
                                              -score =>    $f1score,
                                              -source_tag =>   $f1source,
                                              -primary_tag =>  $f1primary,
                                              -analysis => $analysis_obj );
     
   my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start =>    $f2start,
                                                -end =>      $f2end,
                                                -seqname =>  $f2id,
                                                -strand =>   $f2strand,
                                                -score =>    $f1score,
                                                -source_tag =>   $f2source,
                                                -primary_tag =>  $f2primary,
                                                -analysis => $analysis_obj );
    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
                                             -feature2 => $feat2) ;

    $self->_growfplist($fp); 
}

sub _growfplist {
    my ($self, $fp) =@_;
    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
}

sub _createfiles {
    my ($self, $genfile, $estfile, $dirname)= @_;
    
    #check for diskspace
    my $spacelimit = 0.01; # 0.01Gb or about 10 MB
    my $dir ="./";
    unless ($self->_diskspace($dir, $spacelimit)) 
    {
        $self->throw("Not enough disk space ($spacelimit Gb required)");
    }
            
    #if names not provided create unique names based on process ID    
    $genfile = $self->_getname("genfile") unless ($genfile);
    $estfile = $self->_getname("estfile") unless ($estfile);    
    
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
        $self->throw("Failed to remove @fails : $!\n");
    }
}



sub est_genome {
  my ($self,$arg) = @_;

  if (defined($arg)) {
      $self->{_est_genome} = $arg;
  }
  return $self->{_est_genome};
}
