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

VectorMasker takes a Bio::Seq object and runs blastn with vectors_etc, the
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
    {   $self->blastn('/usr/local/pubseq/bin/blastn');   }
    if ($msp) 
    {   $self->mspcrunch($msp) ;}
    else
    {   $self->mspcrunch('/usr/local/pubseq/bin/MSPcrunch');   }
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
        $seq->isa("Bio::Seq") || $self->throw("Input isn't a Bio::Seq");
        $self->{_clone} = $seq ;
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".VectMask.out");
    }
    return $self->{_clone};
}

sub filename {
    my ($self, $filename) = @_;
    $self->{_filename} = $filename if ($filename);
    return $self->{_filename};
}

sub results {
    my ($self, $results) = @_;
    $self->{_results} = $results if ($results);
    return $self->{_results};
}

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

sub protect {
    my ($self, @filename) =@_;
    push (@{$self->{_protected}}, @filename) if (@filename);
    return @{$self->{_protected}};
}

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
                                                    unless (-e $location);
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
                                                    unless (-e $location);
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

sub workdir {
    my ($self, $directory) = @_;
    if ($directory)
    {
        mkdir ($directory, '777') unless (-d $directory);
        $self->throw ("$directory doesn't exist\n") unless (-d $directory);
        $self->{_workdir} = $directory;
    }
    return $self->{_workdir};
}

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
    $self->workdir('vectormasker') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run genscan       
    $self->run_analysis();
    #parse output and create features
    $self->parse_analysis();
    $self->deletefiles();
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses blastn and MSPcrunch output to give a set of feature pairs
    Returns :   none
    Args    :   optional filename

=cut

sub parsefile {
    my ($self, $filename) = @_;
    $self->results($filename) if ($filename);
    $self->parse_repmask();
}

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

sub parse_analysis {
    my ($self) = @_;
    open (VECTOR, "<".$self->results)
        or $self->throw ("Couldn't open file ".$self->results.": $!\n");
    unless (<VECTOR>)
    {
        print "No hit found with blastn and MSPcrunch\n";
        return;
    }
    while (<VECTOR>)
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
    close VECTOR;       
}

sub createfeaturepair {
    my ($self, $feat1, $feat2) = @_;
    #some contant strings
    
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => $feat2->{db},
                            -db_version      => $feat2->{db_version},
                            -program         => $feat2->{program},
                            -program_version => $feat2->{p_version},
                            -gff_source      => $feat2->{source},
                            -gff_feature     => $feat2->{primary});
    
    #create and fill Bio::EnsEMBL::Seqfeature objects
    my $seqfeature1 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat1->{name},
                            -start   => $feat1->{start},
                            -end     => $feat1->{end},
                            -strand  => $feat1->{strand},
                            -score   => $feat1->{score},
                            -source_tag  => $feat1->{source},
                            -primary_tag => $feat1->{primary},
                            -analysis => $analysis_obj);
    
    my $seqfeature2 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat2->{name},
                            -start   => $feat2->{start},
                            -end     => $feat2->{end},
                            -strand  => $feat2->{strand},
                            -score   => $feat2->{score},
                            -source_tag  => $feat2->{source},
                            -primary_tag => $feat2->{primary},
                            -analysis => $analysis_obj);
    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $seqfeature1,
                                             -feature2 => $seqfeature2) ;
    $self->growfplist($fp);                             
}

sub growfplist {
    my ($self, $fp) =@_;    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
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

sub writefile {
    my ($self) = @_;
    print "Writing sequence to ".$self->filename."\n";
    #create Bio::SeqIO object and save to file
    my $clone_out = Bio::SeqIO->new(-file => ">".$self->filename , '-format' => 'Fasta')
            or $self->throw("Can't create new Bio::SeqIO from ".$self->filename.":$!\n");
    $clone_out->write_seq($self->clone) 
            or $self->throw("Couldn't write to file ".$self->filename.":$!");
}

sub deletefiles {
    my ($self) = @_;
    #delete all analysis files 
    my @list = glob($self->filename."*");
    foreach my $result (@list)
    {
        my $protected = undef; #flag for match found in $protected
        foreach my $suffix ($self->protect)
        {        
            $protected = 'true' if ($result eq $self->filename.$suffix);
        }
        unless ($protected)
        {
            unlink ($result) or $self->throw ("Couldn't delete $result :$!");    
        }
    }
}

sub checkdir {
    my ($self) = @_;
    #check for disk space
    my $spacelimit = 0.01;
    $self->throw("Not enough disk space ($spacelimit required):$!\n") 
                        unless ($self->diskspace('./', $spacelimit));
    my $dir = $self->workdir();
    chdir ($dir) or $self->throw("Cannot change to directory $dir ($!)\n");
}

sub diskspace {
    my ($self, $dir, $limit) =@_;
    my $block_size; #could be used where block size != 512 ?
    my $Gb = 1024 ** 3;
    
    open DF, "df $dir |" or $self->throw ("Can't open 'du' pipe ($!)\n");
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
    close DF || $self->throw("Error from 'df' : $!\n");
}
