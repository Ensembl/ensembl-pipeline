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

Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker

=head1 SYNOPSIS
#create and fill Bio::Seq object
my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq'; 
my $seq = Bio::Seq->new();
my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
$seq = $seqstream->next_seq();
#create Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker object
my $repmask = Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->new (-CLONE => $seq);
$repmask->workdir($workdir);
$repmask->run();
my @results = $repmask->output();
    
=head1 DESCRIPTION
RepeatMasker takes a Bio::Seq object and runs RepeatMaskerHum on it. The
resulting .out file is parsed to produce a set of feature pairs.
Arguments can be passed to RepeatMaskerHum through the arguments() method. 

=head2 Methods:
new($seq_obj)
repeatmasker($path_to_RepeatMaskerHum)
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

package Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker;

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
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->new (-CLONE => $seq);
    Function:   Initialises RepeatMasker object
    Returns :   a RepeatMasker Object
    Args    :   A Bio::Seq object (-CLONE), any arguments for RepeatMaskerHum (-ARGS) 

=cut
sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{_fplist} = [];           #an array of feature pairs
    $self->{_clone}  = undef;        #location of Bio::Seq object
    $self->{_repeatmasker} = undef;  #location of RepeatMaskerHum script
    $self->{_workdir}   = undef;     #location of temp directory
    $self->{_filename}  =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;      #file to store results of RepeatMaskerHum
    $self->{_protected} =[];        #a list of files protected from deletion
    $self->{_arguments} =undef;      #arguments for RepeatMaskerHum
    
    my( $clonefile, $arguments, $repmask) = $self->_rearrange(['CLONE', 'ARGS', 'REPM'], @args);
    $self->clone($clonefile) if ($clonefile);       
    if ($repmask)
    {   $self->repeatmasker($repmask);  }
    else
    {  
        $self->repeatmasker('/tmp_mnt/nfs/disk100/humpub/scripts/RepeatMaskerHum');
    }
    
    if ($arguments) 
    {   $self->arguments($arguments) ;}
    else
    { $self->arguments('-low') ;      }
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
        $self->results($self->filename.".RepMask.out");
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

=head2 repeatmasker
    Title   :   repeatmasker
    Usage   :   $obj->repeatmasker('~humpub/scripts/RepeatMaskerHum');
    Function:   Get/set method for the location of RepeatMaskerHum script
    Args    :   File path (optional)
    
=cut
sub repeatmasker {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("RepeatMaskerHum not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_repeatmasker} = $location ;
    }
    return $self->{_repeatmasker};
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
    Usage   :   $obj->arguments('-init wing -pseudo -caceh -cut 25 -aln 200 -quiet');
    Function:   Get/set method for getz arguments arguments
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
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs RepeatMaskerHum script and creates array of featurepairs
    Returns :   none
    Args    :   optional $workdir and $args (e.g. '-ace' for ace file output)

=cut
sub run {
    my ($self, $dir, $args) = @_;
    #set arguments for repeatmasker
    $self->arguments($args) if ($args);
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for RepeatMasker\n");
    #set directory if provided
    $self->workdir('RepeatMaskerOutput') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile();        
    $self->run_repeatmasker();
    #parse output of repeat masker 
    $self->parse_repmask();
    $self->deletefiles();
}

=head2 parsefile
Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses RepeatMaskerHum output to give a set of feature pairs
    Returns :   none
    Args    :   optional filename

=cut
sub parsefile {
    my ($self, $filename) = @_;
    $self->results($filename) if ($filename);
    $self->parse_repmask();
}

sub run_repeatmasker {
    my ($self) = @_;
    #run RepeatMaskerHum
    print "Running RepeatMasker\n";
    system ($self->repeatmasker.' '.$self->arguments.' '.$self->filename); 
    #or $self->throw("Error running RepeatMasker: $!\n")
    #open repeat predictions
    open (REPOUT, "<".$self->results)
            or $self->throw($self->results." not created by RepeatMaskerHum\n");   
    close REPOUT;
}

sub parse_repmask {
    my ($self) = @_;
    print "Parsing output\n";
    open (REPOUT, "<".$self->results)
        or $self->throw("Error opening ".$self->results."\n");
    #check if no repeats found
    if (<REPOUT> =~ /no repetitive sequences detected/)
    {
        print "RepeatMaskerHum didn't find any repetitive sequences\n";
        close REPOUT;
        return;
    }
    #extract values
    my @output = <REPOUT>;
    for (my $index = 2; $index < scalar(@output); $index++) #loop from 3rd line
    {  
        my @element = split (/\s+/, $output[$index]);  
        my (%feat1, %feat2);
        $feat1 {name} = $element[5];
        $feat1 {score} = $element[1];
        $feat1 {start} = $element[6];
        $feat1 {end} = $element[7];
        #set strand ('+' = 1 and 'c' = -1)
        $feat2 {strand} = 1 if ($element[9] eq '+');
        $feat2 {strand} = -1 if ($element[9] eq 'C');
        $feat2 {name} = $element[10];
        $feat2 {score} = $element[1];
        ($element[13] =~ tr/()//d); #strip away parentheses from $element[13]
        ($element[14] =~ tr/()//d); #strip away parentheses from $element[14]
        $feat2 {start} = $element[13];     
        $feat2 {end} = $element[14];
        $feat1 {strand} = 1;
        $self->createfeaturepair(\%feat1, \%feat2); #may need to use references
    }
    close REPOUT;   
}

sub createfeaturepair {
    my ($self, $feat1, $feat2) = @_;
    #some contant strings
    my $source = 'RepeatMaskerHum';
    my $primary = 'similarity';
    
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => $source,
                            -program_version => "unknown",
                            -gff_source      => $source,
                            -gff_feature     => $primary,);
    
    #create and fill Bio::EnsEMBL::Seqfeature objects
    my $seqfeature1 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat1->{name},
                            -start   => $feat1->{start},
                            -end     => $feat1->{end},
                            -strand  => $feat1->{strand},
                            -score   => $feat1->{score},
                            -source_tag  => $source,
                            -primary_tag => $primary,
                            -analysis => $analysis_obj);
    
    my $seqfeature2 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat2->{name},
                            -start   => $feat2->{start},
                            -end     => $feat2->{end},
                            -strand  => $feat2->{strand},
                            -score   => $feat2->{score},
                            -source_tag  => $source,
                            -primary_tag => $primary,
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
