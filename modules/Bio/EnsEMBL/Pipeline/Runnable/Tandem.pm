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

Bio::EnsEMBL::Pipeline::Runnable::Tandem

=head1 SYNOPSIS

 #create and fill Bio::Seq object
 my $seq = Bio::Seq->new();
 my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
 $seq = $seqstream->next_seq();

 #create Bio::EnsEMBL::Pipeline::Runnable::Tandem object
 my $tandem = Bio::EnsEMBL::Pipeline::Runnable::Tandem->new (-CLONE => $seq);
 $tandem->workdir($workdir);
 $tandem->run();
 @featurepairs = $tandem->output();

=head1 DESCRIPTION

Tandem takes a Bio::Seq (or Bio::PrimarySeq) object and runs equicktandem and etandem. The run strategy is taken from
tan_Search in hp.pl. equicktandem is run to get a range of sizes from 2-3. Theses sizes are used as the
only repeat length searched by etandem. The output is added to the output of etandem run with a repeat
length of 2 and the uniform switch set to 'yes'. The results are parsed to produce a set
ofBio::ensEMBL::FeaturePairs. 

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

package Bio::EnsEMBL::Pipeline::Runnable::Tandem;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Programs;
use Bio::Root::RootI;
use Bio::Seq;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Tandem->new (-CLONE => $seq);
    Function:   Initialises Tandem object
    Returns :   a Tandem Object
    Args    :   A Bio::Seq object (-CLONE), any arguments for eTandem (-ARGS) 

=cut

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{_fplist} = [];           #an array of feature pairs
    $self->{_clone}  = undef;        #location of Bio::Seq object
    $self->{_etandem} = undef;       #location of etandem
    $self->{_equicktandem} = undef;  #location of equicktandem
    $self->{_workdir}   = undef;     #location of temp directory
    $self->{_filename}  =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;      #file to store results of analysis
    $self->{_protected} =[];         #a list of files protected from deletion
    $self->{_arguments} =undef;      #arguments for etandem
    
    
    my( $clonefile, $arguments, $etandem, $equicktandem) = 
            $self->_rearrange(['CLONE', 'ARGS', 'ETAND', 'EQTAND'], @args);
    $self->clone($clonefile) if ($clonefile);       
    if ($etandem) 
    {   $self->etandem($etandem) ;}
    else
    {   
        eval 
        { $self->etandem($self->locate_executable('etandem')); };
        if ($@)
        { $self->etandem('/nfs/disk100/pubseq/emboss/osf/EMBOSS-0.0.4/emboss/etandem'); }
    }
    if ($equicktandem) 
    {   $self->etandem($equicktandem) ;}
    else
    {   
        eval 
        { $self->equicktandem($self->locate_executable('equicktandem')); };
        if ($@)
        { $self->equicktandem('/nfs/disk100/pubseq/emboss/osf/EMBOSS-0.0.4/emboss/equicktandem'); }
    }
    if ($arguments) 
    {   $self->arguments($arguments) ;}
    else
    { $self->arguments(' -mismatch yes -threshold 20 ') ;      }
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
        $self->results($self->filename.".tandem.out");
    }
    return $self->{_clone};
}

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

=head2 etandem

    Title   :   etandem
    Usage   :   $obj->etandem('/nfs/disk100/pubseq/emboss/osf/EMBOSS-0.0.4/emboss/etandem')
    Function:   Get/set method for the location of etandem
    Args    :   File path (optional)

=cut

sub etandem {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("etandem not found at $location: $!\n") 
                                                    unless (-e $location && -x $location);
        $self->{_etandem} = $location ;
    }
    return $self->{_etandem};
}

=head2 equicktandem

    Title   :   equicktandem
    Usage   :   $obj->equicktandem('/nfs/disk100/pubseq/emboss/osf/EMBOSS-0.0.4/emboss/equicktandem')
    Function:   Get/set method for the location of equicktandem
    Args    :   File path (optional)

=cut

sub equicktandem {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("equicktandem not found at $location: $!\n") 
                                                    unless (-e $location && -x $location);
        $self->{_equicktandem} = $location ;
    }
    return $self->{_equicktandem};
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
    Function:   Get/set method for etandem arguments
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
    Function:   Runs equicktandem and etandem to create array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir) = @_;
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Tandem\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run genscan       
    $self->run_tandem();
    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
}

=head2 clonename

    Title   :   clonename
    Usage   :   $obj->clonename('AC00074');
    Function:   Get/set method for clone name. 
                This must be set manually when a file or pipe is parsed and the clonename is 
                not present in the executable output
    Args    :   File suffixes

=cut

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses etandem output to give a set of feature pairs
    Returns :   none
    Args    :   optional filename

=cut

sub run_tandem {
    my ($self) = @_;
    #emulating tan_search used by run_tandem from hp.pl. This long because many files are created and read
    #tan_Search runs quicktandem to get repeat lengths and runs tandem with these lengths, it then runs tandem with a length of 2 with the -U switch
    my $tandem_repeats = '2';#keeps track of repeat lengths produced from equicktandem - default run uses 2
    my @results;             #stores all hits from all runs of etandem
    print STDERR "Running equicktandem on ".$self->filename."\n";
    system ($self->equicktandem().' -threshold 20 -maxrepeat 3 -sequence '.$self->filename()." -outfile ".$self->results().'.qt');
    open (QUICK, "<".$self->results.'.qt');
    while (<QUICK>)
    {
        my @element = split;
        chomp (@element);
        if ($element[3] && $element[3] > 1 && $element[3] < 4)
        {
            $tandem_repeats = join (' ', $tandem_repeats, $element[3]) 
                            unless ($tandem_repeats =~ /\b$element[3]\b/);        
        }
        
    }
    close QUICK;
    unlink ($self->results.'.qt'); #delete temporary quicktandem file
    #now run etandem on every repeat length picked up by equicktandem
    print STDERR "Running etandem on ".$self->filename."\n";
    my @repeatLength = split (/\s+/, $tandem_repeats);
    foreach my $repeat (@repeatLength)
    {
        system ($self->etandem." ".$self->arguments." -minrepeat $repeat -maxrepeat $repeat"
                ." -sequence ".$self->filename ." -outfile ".$self->results.$repeat);
        open (OUT, "<".$self->results.$repeat) 
                or $self->throw ("Error opening putput of etandem $repeat: $!\n");
        my @output = <OUT>;
        close OUT;
        push (@results, @output);
        unlink ($self->results.$repeat);    
    }
    system ($self->etandem." ".$self->arguments." -uniform yes -minrepeat 2 -maxrepeat 2"
                ." -sequence ".$self->filename ." -outfile ".$self->results);
    open (RESULTS, ">>".$self->results);
    print RESULTS @results;
    close RESULTS;
}

sub parse_results {
    my ($self) = @_;
    my $filehandle;
    if (ref ($self->results) !~ /GLOB/)
    { 
        open (TANDEM, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*TANDEM;
    }
    else
    {
        $filehandle = $self->results;
    } 
    
    unless (<$filehandle>)
    {
        print "No repeats found with etandem\n";
        return;
    }
    while (<$filehandle>)
    {
        my @element = split;
        my (%feat1, %feat2);
        if ($self->clonename)
            {
                $feat1{name}     = $self->clonename;
            }
            else
            {
                $self->results =~ m!/.+/(.+)|(.+)!; #extract filename
                #($1) ? $feat1{name} = $1 : $feat1{name} = $2;
                if ($1) { $feat1{name} = $1; }
                elsif ($2) { $feat1{name} = $2; }
            }
        $feat1 {score} = $element[0];
        $feat1 {start} = $element[1];
        $feat1 {end} = $element[2];
        $feat1 {strand} = 1;
        $feat2 {strand} = 1;
        $feat2 {name} = $element[6];
        $feat2 {score} = $element[5];
        $feat2 {start} = 1;
        $feat2 {end} = $feat1{end} - $feat1{start};
        $feat1 {primary} = 'similarity';
        $feat1 {source} = 'etandem';
        $feat2 {primary} = 'similarity';
        $feat2 {source} = 'etandem';
        $feat2 {db} = undef;
        $feat2 {db_version} = undef;
        $feat2 {program} = 'etandem';
        $feat2 {p_version} = 'unknown';
        
        $self->createfeaturepair(\%feat1, \%feat2);
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

