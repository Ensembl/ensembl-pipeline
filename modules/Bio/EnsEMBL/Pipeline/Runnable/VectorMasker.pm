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

Bio::EnsEMBL::Pipeline::Runnable::VectorMasker

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq';
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::VectorMasker object
  my $vectmask = Bio::EnsEMBL::Pipeline::Runnable::VectorMasker->new (-CLONE => $seq);
  $vectmask->workdir($workdir);
  $vectmask->run();
  @featurepairs = $vectmask->output();

=head1 DESCRIPTION

VectorMasker takes a Bio::Seq (or Bio::PrimarySeq) object and runs blastn with vectors_etc, the
output is parsed by MSPcrunch and stored as Bio::EnsEMBL::FeaturePairs. 
Arguments can be passed to MSPcrunch through the arguments() method. 

=head2 Methods:

new($seq_obj)
blastn($path_to_blastn)
mspcrunch($path_to_mspcrunch)
vector($path_to_vector);
workdir($directory_name)
arguments($args)
run()
output()

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::VectorMasker;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object );

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::VectorMasker->new (-CLONE => $seq);
    Function:   Initialises VectorMasker object
    Returns :   a VectorMasker Object
    Args    :   A Bio::Seq object (-CLONE), any arguments for blastn (-ARGS) 

=cut

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{_fplist} = [];           #an array of feature pairs
    $self->{_clone}  = undef;        #location of Bio::Seq object
    $self->{_blastn} = undef;        #location of Blastn script
    $self->{_mspcrunch} = undef;     #location of MSPcrunch
    $self->{_vector} = undef;       #location of vector sequences
    $self->{_workdir}   = undef;     #location of temp directory
    $self->{_filename}  =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;      #file to store results of analysis
    $self->{_protected} =[];         #a list of files protected from deletion
    $self->{_arguments} =undef;      #arguments for MSPcrunch
    
    
    my( $clonefile, $arguments, $blastn, $msp, $vector) = 
            $self->_rearrange(['CLONE', 'ARGS', 'BLAST','MSP', 'VECT'], @args);
    $self->clone($clonefile) if ($clonefile);       
    if ($blastn) 
    {   $self->blastn($blastn) ;}
    else
    {   
        eval 
        { $self->blastn($self->locate_executable('blastn'));  }; 
        if ($@) 
        { $self->blastn('/usr/local/pubseq/bin/blastn'); }  
    }
    if ($msp) 
    {   $self->mspcrunch($msp) ;}
    else
    {   
        eval 
        { $self->mspcrunch($self->locate_executable('MSPcrunch')); }; 
        if ($@)
        { $self->mspcrunch('/usr/local/pubseq/bin/MSPcrunch'); }  
    }
    if ($vector) 
    {   $self->vector($vector) ;}
    else
    {   $self->vector('/tmp_mnt/nfs/disk100/humpub/blast/vectors_etc');   }
    if ($arguments) 
    {   $self->arguments($arguments) ;}
    else
    { $self->arguments(' -I 95 -d -') ;      }
    return $self; # success - we hope!
}

#################
# get/set methods 
#################

sub clone {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeq") || $seq->isa("Bio::Seq")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{_clone} = $seq ;
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".VectMask.out");
    }
    return $self->{_clone};
}

=head2 clonename

    Title   :   clonename
    Usage   :   $obj->clonename('AC00074');
    Function:   Get/set method for clone name. 
                This must be set manually when a file or pipe is parsed and the clonename is 
                not present in the executable output
    Args    :   File suffixes

=cut

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

=head2 blastn

    Title   :   blastn
    Usage   :   $obj->blastn('/usr/local/pubseq/bin/blastn');
    Function:   Get/set method for the location of blastn
    Args    :   File path (optional)

=cut

sub blastn {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("blastn not found at $location: $!\n") 
                                                    unless (-e $location && -x $location);
        $self->{_blastn} = $location ;
    }
    return $self->{_blastn};
}

=head2 mspcrunch

    Title   :   mspcrunch
    Usage   :   $obj->mspcrunch('/usr/local/pubseq/bin/MSPcrunch');
    Function:   Get/set method for the location of mspcrunch
    Args    :   File path (optional)

=cut

sub mspcrunch {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("MSPcrunch not found at $location: $!\n") 
                                                    unless (-e $location && -x $location);
        $self->{_mspcrunch} = $location ;
    }
    return $self->{_mspcrunch};
}

=head2 vector

    Title   :   vector
    Usage   :   $obj->vector('/tmp_mnt/nfs/disk100/humpub/blast/vectors_etc');
    Function:   Get/set method for the location of vector
    Args    :   File path (optional)

=cut

sub vector {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("vector_etc not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_vector} = $location ;
    }
    return $self->{_vector};
}

=head2 workdir

    Title   :   workdir
    Usage   :   $obj->wordir('~humpub/temp');
    Function:   Get/set method for the location of a directory to contain temp files
    Args    :   File path (optional)

=cut

=head2 arguments

    Title   :   arguments
    Usage   :   $obj->arguments(' -I 95 -x -');
    Function:   Get/set method for MSPcrunch arguments
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

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blastn and MSPcrunch and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir) = @_;
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Vectormasker\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run genscan       
    $self->run_analysis();
    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses blastn and MSPcrunch output to give a set of feature pairs
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub run_analysis {
    my ($self) = @_;
    print "Running Blastn and MSPcrunch on ".$self->filename."\n";
    open (OUTPUT, $self->blastn.' '.$self->vector.' '.$self->filename.' | '.
                    $self->mspcrunch.' '.$self->arguments.'|')
            or $self->throw("Couldn't open pipe to blastn & MSPcrunch: $!\n");  
    open (RESULTS, ">".$self->results)
            or $self->throw("Couldn't create file for Vectormasker results: $!\n");
    print RESULTS <OUTPUT>;
    close OUTPUT;
    close RESULTS;
}

sub parse_results {
    my ($self) = @_;
    my $filehandle;
    
    if (ref ($self->results) !~ /GLOB/)
    { 
        open (VECTOR, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*VECTOR;
    }
    else
    {
        $filehandle = $self->results;
    }
    
    unless (<$filehandle>)
    {
        print "No hit found with blastn and MSPcrunch\n";
        return;
    }
    while (<$filehandle>)
    {
        my @element = split;
        my (%feat1, %feat2);
        $feat1 {name} = $element[4];
        $feat1 {score} = $element[1];
        
        if ($element[2] < $element[3])
        {
            $feat1 {start} = $element[2];
            $feat1 {end} = $element[3];
            $feat1 {strand} = 1;
            $feat2 {strand} = 1;
        }
        else
        {
            $feat1 {start} = $element[3];
            $feat1 {end} = $element[2];
            $feat1 {strand} = 1;
            $feat2 {strand} = -1;
        }
        
        $feat2 {name} = $element[7];
        $feat2 {score} = $element[1];
        $feat2 {start} = $element[5];
        $feat2 {end} = $element[6];
        $feat1 {primary} = 'similarity';
        $feat1 {source} = 'VectorMasker';
        $feat2 {primary} = 'similarity';
        $feat2 {source} = 'vector';
        
        $feat2 {db} = 'vector_etc';
        $feat2 {db_version} = 'unknown';
        $feat2 {program} = 'VectorMasker';
        $feat2 {p_version} = 'unknown';
        
        $self->createfeaturepair(\%feat1, \%feat2); #may need to use references
    }
    close $filehandle;       
}



##############
# input/output methods
#############
=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of feature pairs
    Returns :   Returns an array of feature pairs
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_fplist'}};
}

