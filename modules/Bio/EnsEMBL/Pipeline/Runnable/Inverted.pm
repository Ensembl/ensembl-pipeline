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

Bio::EnsEMBL::Pipeline::Runnable::Inverted

=head1 SYNOPSIS

 #create and fill Bio::Seq object
 my $seq = Bio::Seq->new();
 my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
 $seq = $seqstream->next_seq();

 #create Bio::EnsEMBL::Pipeline::Runnable::Inverted object
 my $inverted = Bio::EnsEMBL::Pipeline::Runnable::Inverted->new (-CLONE => $seq);
 $inverted->workdir($workdir);
 $inverted->run();
 @featurepairs = $inverted->output();

=head1 DESCRIPTION

Inverted takes a Bio::Seq (or Bio::PrimarySeq) object and runs einverted. The
resulting inverted repeats are used to create a pair of Bio::FeaturePairs per
repeat. The output is returned as an array of feature pairs.  

=head2 Methods:

 new($seq_obj)
 workdir($directory_name)
 arguments($args)
 parsefile()
 run()
 output()

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Inverted;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Programs;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Inverted->new (-CLONE => $seq);
    Function:   Initialises Inverted object
    Returns :   an Inverted Object
    Args    :   A Bio::Seq object (-CLONE), any arguments for einverted (-ARGS) 

=cut

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{_fplist} = [];           #an array of feature pairs
    $self->{_clone}  = undef;        #location of Bio::Seq object
    $self->{_einverted} = undef;       #location of einverted
    $self->{_workdir}   = undef;     #location of temp directory
    $self->{_filename}  =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;      #file to store results of analysis
    $self->{_protected} =[];         #a list of files protected from deletion
    $self->{_arguments} =undef;      #arguments for einverted
    
    
    my( $clonefile, $arguments, $einverted, $equickinverted) = 
            $self->_rearrange(['CLONE', 'ARGS', 'ETAND', 'EQTAND'], @args);
    $self->clone($clonefile) if ($clonefile);       
    if ($einverted) 
    {   $self->einverted($einverted) ;}
    else
    {   
        eval 
        { $self->einverted($self->locate_executable('einverted')); };
        if ($@)
        { $self->einverted('/nfs/disk100/pubseq/emboss/osf/EMBOSS-0.0.4/emboss/einverted'); }
    }
    if ($arguments) 
    {   $self->arguments($arguments) ;}
    else
    { $self->arguments(' -mismatch -4 -threshold 50 -gap 12 -match 3 ') ;      }
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
        $self->clonename($self->clone->id);
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".inverted.out");
    }
    return $self->{_clone};
}

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

=head2 einverted

    Title   :   einverted
    Usage   :   $obj->einverted('/nfs/disk100/pubseq/emboss/osf/EMBOSS-0.0.4/emboss/einverted')
    Function:   Get/set method for the location of einverted
    Args    :   File path (optional)

=cut

sub einverted {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("einverted not found at $location: $!\n") 
                                                    unless (-e $location && -x $location);
        $self->{_einverted} = $location ;
    }
    return $self->{_einverted};
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
    Function:   Get/set method for einverted arguments
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

=head2 clonename

    Title   :   clonename
    Usage   :   $obj->clonename('AC00074');
    Function:   Get/set method for clone name. 
                This must be set manually when a file or pipe is parsed and the clonename is 
                not present in the executable output
    Args    :   File suffixes

=cut

###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs einverted to create array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir) = @_;
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Inverted\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run genscan       
    $self->run_inverted();
    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses einverted output to give a set of feature pairs
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub run_inverted {
    my ($self) = @_;
    system ($self->einverted.' '.$self->arguments.' -sequence '.$self->filename.' -outfile '.$self->results);
    $self->throw ("Error running einverted\n") unless (-e $self->results());
}

sub parse_results {
    my ($self) = @_;
    my (%feat1, %feat2);
    my $filehandle;
    if (ref ($self->results) !~ /GLOB/)
    {
        open (INVERTED, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*INVERTED;
    }
    else
    {
        $filehandle = $self->results;
    }
    
    unless (<$filehandle>)
    {
        print STDERR "No repeats found with einverted\n";
        return;
    }
    while (<$filehandle>)
    {
        if (m|^Score\s+(\d+)\D+(\d+)/(\d+)\D+(\d+)|i)
        {
            if ($self->clonename)
            {
                $feat1 {name} = $self->clonename;
            }
            else
            {
                $self->results =~ m!/.+/(.+)|(.+)!; #extract filename
                #($1) ? $feat1{name} = $1 : $feat1{name} = $2;
                if ($1) { $feat1{name} = $1; }
                elsif ($2) { $feat1{name} = $2; }
            }
            $feat2 {name} = $feat1{name}."_inv_repeat";
            $feat1 {score} = $1; #score as an inverted repeat
            $feat2 {score} = $1; 
            $feat1 {percent} = $4; #percentage identity
            $feat2 {percent} = $4;
        }
        elsif (/(\d+)\D+(\d+)/ && ($feat1{start}))
        { 
            if ($1 < $2)
            {
                $feat2 {start} = $1;
                $feat2 {end} = $2;
                $feat2 {strand} = 1;
            }
            else
            {
                $feat2 {start} = $2;
                $feat2 {end} = $1;
                $feat2 {strand} = -1;
            }
            $feat1 {primary} = 'repeat';
            $feat1 {source} = 'einverted';
            $feat2 {primary} = 'repeat';
            $feat2 {source} = 'einverted';
            $feat2 {db} = undef;
            $feat2 {db_version} = undef;
            $feat2 {program} = 'einverted';
            $feat2 {p_version} = 'unknown';
            $self->create_repeat(\%feat1, \%feat2);
            #reverse scores and names
            ($feat1 {name}, $feat2 {name}) = ($feat2 {name}, $feat1 {name});
            ($feat1 {score}, $feat2 {score}) = ($feat2 {score}, $feat1 {score});
            $self->create_repeat(\%feat2, \%feat1);
            #reset flag variables to allow correct reading of next repeat
            $feat1{name} = undef;
            $feat1{start} = undef;
        }
        elsif (/(\d+)\D+(\d+)/ && ($feat1{name}))
        {
            if ($1 < $2)
            {
                $feat1 {start} = $1;
                $feat1 {end} = $2;
                $feat1 {strand} = 1;
            }
            else
            {
                $feat1 {start} = $2;
                $feat1 {end} = $1;
                $feat1 {strand} = -1;
            }
        }
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

