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

Bio::EnsEMBL::Pipeline::Runnable::Blast

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq';
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::Blast object
  my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (-CLONE => $seq);
  $blast->workdir($workdir);
  $blast->run();
  @featurepairs = $blast->output();

=head1 DESCRIPTION

Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs blast with, the
output is parsed by BPLite and stored as Bio::EnsEMBL::FeaturePairs. 
Arguments can be passed to BPLite through the arguments() method. 

=head2 Methods:

    new($seq_obj)
    blast($path_to_blast)
    BPLite($path_to_bplite)
    database($path_to_database);
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

package Bio::EnsEMBL::Pipeline::Runnable::Blast;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use BPlite;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object );

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new (-CLONE => $seq);
    Function:   Initialises Blast object
    Returns :   a Blast Object
    Args    :   A Bio::Seq object (-CLONE), any arguments for blast (-ARGS)
                The blast executable (-BLAST) and database (-DB).

=cut

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{_fplist} = [];           #an array of feature pairs
    $self->{_clone}  = undef;        #location of Bio::Seq object
    $self->{_blast}  = undef;        #location of Blast
    $self->{_database} = undef;      #name and location of database
    $self->{_workdir}   = undef;     #location of temp directory
    $self->{_filename}  =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;      #file to store results of analysis
    $self->{_protected} =[];         #a list of files protected from deletion
    $self->{_arguments} =undef;      #arguments for blast
      
    my( $clone, $blast, $database, $arguments, $threshold) = 
            $self->_rearrange(['CLONE', 'BLAST', 'DB', 'ARGS', 'THRESHOLD'], @args);
    
    $self->clone($clone) if ($clone);       
    
    if ($blast =~ m!/!) #path to blast is provided 
    {   $self->blast($blast) ;}
    elsif ($blast =~ /blastn|blastx|blastp|tblastn|tblastx/)
    {   
        eval 
        { $self->blast($self->locate_executable($blast));  }; 
        if ($@) 
        { $self->blast('/usr/local/pubseq/bin/'.$blast); }  
    }
    
    if ($database) 
    {   $self->database($database) ;}
    else
    {   $self->database('');   }
    
    if ($arguments) 
    {   $self->arguments($arguments) ;}
    else
    { $self->arguments(' ');   }
    
    if ($threshold) 
    {   $self->threshold($threshold) ;}
    else
    { $self->threshold(0);     }
    
    return $self; # success - we hope!
}

#################
# get/set methods 
#################

sub clone {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{_clone} = $seq ;
        $self->clonename($self->clone->id);
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".blast.out");
    }
    return $self->{_clone};
}


=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

=head2 threshold

    Title   :   threshold
    Usage   :   $obj->threshold($value);
    Function:   Get/set method for threshold p_value required for outputting
                Feature/FeaturePair
    Args    :   Optional value (blast p_value)

=cut

=head2 blast

    Title   :   blast
    Usage   :   $obj->blast('/usr/local/pubseq/bin/blastn');
    Function:   Get/set method for the location of blast
    Args    :   none

=cut

sub blast {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("executable not found at $location: $!\n") 
                                                    unless (-e $location && -x $location);
        $self->{_blast} = $location ;
    }
    return $self->{_blast};
}

=head2 database

    Title   :   database
    Usage   :   $obj->database('dbEST');
    Function:   Get/set method for the location of database
    Args    :   none

=cut

sub database {
    my ($self, $db) = @_;
    if ($db)
    {
        $self->{_database} = $db ;
    }
    return $self->{_database};
}

=head2 workdir

    Title   :   workdir
    Usage   :   $obj->wordir('~humpub/temp');
    Function:   Get/set method for the location of a directory to contain temp files
    Args    :   File path (optional)

=cut

=head2 arguments

    Title   :   arguments
    Usage   :   $obj->arguments(' -I ');
    Function:   Get/set method for blast arguments
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
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir) = @_;
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Blast\n");
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

    Title   :   parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses blast and BPLite output to give a set of feature pairs
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub run_analysis {
    my ($self) = @_;
    
    print STDERR ("Running blast and BPlite:\n ".$self->blast.' '
                    .$self->database.' '.$self->filename.' '
                    .$self->arguments.' > ' .$self->results."\n");
    
    system ($self->blast.' '.$self->database.' '.$self->filename
            .' '.$self->arguments.' > '.$self->results);
    
    $self->throw("Couldn't create file for Blast results: $!\n") 
                unless (-e $self->results);
}

#New and improved! takes filenames and handles, therefore pipe compliant!
sub parse_results {
    my ($self) = @_;
    my $filehandle;
    if (ref ($self->results) !~ /GLOB/)
    { 
        open (BLAST, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*BLAST;
    }
    else
    {
        $filehandle = $self->results;
    }    
    
    unless (<$filehandle>)
    {
        print "No hit found with blast \n";
        return;
    }
    my $parser = new BPlite ($filehandle);
    $parser->query();
    $parser->database();
    my $threshold = $self->threshold();
    while(my $sbjct = $parser->nextSbjct) 
    {
       while (my $hsp = $sbjct->nextHSP)
       {
            my (%feat1, %feat2);
            if ($self->clonename)
            {
                $feat1{name}     = $self->clonename;
            }
            else
            {
                $self->results =~ m!/.+/(.+)|(.+)!; #extract filename
                $feat1{name} = ($1) ?  $1 :  $2;
            }
            
            $feat1 {score}       = $hsp->score;
            $feat2 {score}       = $hsp->score;
            $feat1 {percent}     = $hsp->percent;
            $feat2 {percent}     = $hsp->percent;
            $feat1 {p}           = $hsp->P;
            $feat2 {p}           = $hsp->P;
            
            unless ($hsp->queryBegin > $hsp->queryEnd)
            {
                $feat1 {start}   = $hsp->queryBegin;
                $feat1 {end}     = $hsp->queryEnd;
                $feat1 {strand}  = 1;
            }
            else
            {
                $feat1 {end}     = $hsp->queryBegin;
                $feat1 {start}   = $hsp->queryEnd;
                $feat1 {strand}  = -1;
            }
            
            $sbjct->name =~ /[\||\s|:](\w+)[\||\s|:]/; #extract subjectname
            $feat2 {name}    = $1;
            
            unless ($hsp->sbjctBegin > $hsp->sbjctEnd)
            {
                $feat2 {start}   = $hsp->sbjctBegin;
                $feat2 {end}     = $hsp->sbjctEnd;
                $feat2 {strand}  = 1;
            }
            else
            {
                $feat2 {end}     = $hsp->sbjctBegin;
                $feat2 {start}   = $hsp->sbjctEnd;
                $feat2 {strand}  = -1;
            }
            
            if ($self->database)
            {
                $feat2 {db} = $self->database;
            }
            else
            {
                $feat2 {db} = 'unknown';
            }
                        
            if ($self->blast)
            {
                $self->blast =~ m!/.+/(.+)|(.+)!; #extract executable name
                if ($1) 
                    { 
                    $feat2 {program} = $1; }
                elsif ($2) 
                    { $feat2 {program} = $2; }
            }
            else
            {
                $feat2 {program} = 'unknown';
            }
            $feat2 {p_version} = '1';
            $feat2 {db_version} = '1';
            $feat1 {primary} = 'similarity';
            $feat1 {source}  = $feat2{program};
            $feat2 {primary} = 'similarity';
            $feat2 {source}  = $feat2{program};
            
            if ( defined $threshold && ($feat1{p} >= $threshold))
            { 
                $self->createfeaturepair(\%feat1, \%feat2); 
            }
            else
            {
                print STDERR "Discarding ".$feat2{name}." p ".$feat1{p}." threshold ($threshold)\t";
            }
        }

    }
           
}

# FeaturePair methods:
#sub createfeaturepair {
#    my ($self, $feat1, $feat2) = @_;
#    #call Bio::EnsEMBL::Pipeline::RunnableI version
#    $self->SUPER::createfeaturepair($feat1, $feat2);
#    my $fp = $self->shrinkfplist();
#    $fp->percent_id     ($feat1->{percent});
#    $fp->hpercent_id    ($feat2->{percent});
#    $fp->p_value        ($feat1->{p});
#    $fp->hp_value       ($feat2->{p});
#    $self->growfplist($fp);                             
#}

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
