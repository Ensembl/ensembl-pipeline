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

=over4

=item new($seq_obj)

=item repeatmasker($path_to_RepeatMaskerHum)

=item workdir($directory_name)

=item arguments($args)

=item run()

=item output()

=back

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
    $self->{_protected} =[];         #a list of files protected from deletion
    $self->{_arguments} =undef;      #arguments for RepeatMaskerHum
    
    my( $clonefile, $arguments, $repmask) = $self->_rearrange(['CLONE', 'ARGS', 'REPM'], @args);
    $self->clone($clonefile) if ($clonefile);       
    if ($repmask)
    {   $self->repeatmasker($repmask);  }
    else
    {   
        eval 
        { $self->repeatmasker($self->locate_executable('RepeatMaskerHum')); };
        if ($@)
        { $self->repeatmasker('~humpub/scripts/RepeatMaskerHum'); }
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

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

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
    $self->workdir('tmp') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile();        
    $self->run_repeatmasker();
    #parse output of repeat masker 
    $self->parse_results();
    $self->deletefiles();
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses RepeatMaskerHum output to give a set of feature pairs
    Returns :   none
    Args    :   optional filename

=cut

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

sub parse_results {
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
        next if ($element[12-14] =~ /-/); # ignore features with negatives
        my (%feat1, %feat2);
        $feat1 {name} = $element[5];
        $feat1 {score} = $element[1];
        $feat1 {start} = $element[6];
        $feat1 {end} = $element[7];
        #The start and end values are in different columns depending on orientation!
        if ($element[9] eq '+')
        {
            $feat2 {strand} = 1;
            $feat2 {start} = $element[12];     
            $feat2 {end} = $element[13];
        }
        elsif ($element[9] eq 'C')
        {
            $feat2 {strand} = -1 ;
            $feat2 {start} = $element[14];     
            $feat2 {end} = $element[13];
        }
        $feat2 {name} = $element[10];
        $feat2 {score} = $element[1];
        $feat1 {strand} = 1;
        #misc
        $feat2 {db} = undef;
        $feat2 {db_version} = undef;
        $feat2 {program} = 'RepeatMaskerHum';
        $feat2 {p_version}='unknown';
        $feat2 {source}= 'ReapeatMaskerHum';
        $feat2 {primary}= 'similarity';
        $feat1 {source}= 'RepearMaskerHum';
        $feat1 {primary}= 'similarity';
        
        $self->createfeaturepair(\%feat1, \%feat2); #may need to use references
    }
    close REPOUT;   
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

