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

RepeatMasker takes a Bio::Seq (or Bio::PrimarySeq) object and runs RepeatMasker on it. The
resulting .out file is parsed to produce a set of feature pairs.
Arguments can be passed to RepeatMasker through the arguments() method. 

=head2 Methods:

=over4

=item new($seq_obj)

=item repeatmasker($path_to_RepeatMasker)

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

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->new (-CLONE => $seq);
    Function:   Initialises RepeatMasker object
    Returns :   a RepeatMasker Object
    Args    :   A Bio::Seq object (-CLONE), any arguments for RepeatMasker (-ARGS) 

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    $self->{'_fplist'} = [];           #an array of feature pairs
    $self->{'_clone'}  = undef;        #location of Bio::Seq object
    $self->{'_repeatmasker'} = undef;  #location of RepeatMasker executable
    $self->{'_workdir'}   = undef;     #location of temp directory
    $self->{'_filename'}  =undef;      #file to store Bio::Seq object
    $self->{'_results'}   =undef;      #file to store results of RepeatMasker
    $self->{'_protected'} =[];         #a list of files protected from deletion
    $self->{'_arguments'} =undef;      #arguments for RepeatMasker
    
    my( $clone, $arguments, $repmask) = $self->_rearrange([qw(CLONE
							      ARGS
							      REPM)], 
							  @args);
    
    $self->clone($clone) if ($clone);       

    my $bindir = $::pipeConf{'bindir'} || undef;

    if (-x $repmask) {
        # passed from RunnableDB (full path assumed)
        $self->repeatmasker($repmask);
    }
    elsif ($::pipeConf{'bin_RepeatMasker'} && -x ($repmask = "$::pipeConf{'bin_RepeatMasker'}")) {
        $self->repeatmasker($repmask);
    }
    elsif ($bindir && -x ($repmask = "$bindir/RepeatMasker")) {
        $self->repeatmasker($repmask);
    }
    else {
        # search shell $PATH
        eval {
            $self->repmask($self->locate_executable('RepeatMasker'));
        };
        if ($@) {
            $self->throw("Can't find executable RepeatMasker");
        }
    }

    if ($arguments) {
	$self->arguments($arguments);
    }
    else {
	$self->arguments('');
    }
    return $self;
}

#################
# get/set methods 
#################

sub clone {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
        }
        $self->{'_clone'} = $seq ;
        
        $self->clonename($self->clone->id);
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".out");
    }
    return $self->{'_clone'};
}

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

=head2 repeatmasker

    Title   :   repeatmasker
    Usage   :   $obj->repeatmasker('~humpub/scripts/RepeatMasker');
    Function:   Get/set method for the location of RepeatMasker script
    Args    :   File path (optional)

=cut

sub repeatmasker {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("RepeatMasker not found at $location: $!\n") 
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

=head2 clonename

    Title   :   clonename
    Usage   :   $obj->clonename('AC00074');
    Function:   Get/set method for clone name. 
                This must be set manually when a file or pipe is parsed and the clonename is 
                not present in the executable output
    Args    :   File suffixes

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
        $self->{'_arguments'} = $args ;
    }
    return $self->{'_arguments'};
}
###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs RepeatMasker script and creates array of featurepairs
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
    $self->workdir('/tmp') unless ($self->workdir($dir));
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
    Function:   Parses RepeatMasker output to give a set of feature pairs
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub run_repeatmasker {
    my ($self) = @_;
    #run RepeatMasker
    print "Running RepeatMasker; ";
    print "cmd: ", $self->repeatmasker.' '.$self->arguments.' '.$self->filename, "\n";
    $self->throw("Error running RepeatMasker on ".$self->filename."\n") 
        if (system ($self->repeatmasker.' '.$self->arguments.' '.$self->filename)); 
}

#New and improved! takes filenames and handles, therefore pipe compliant!
sub parse_results {
    my ($self) = @_;
    my $filehandle;
    #The results filename is different in RepeatMasker (8/14/2000) and RepeatMaskerHum
    my $results = (-e $self->results) ? $self->results : $self->filename.".out";
    if (ref ($self->results) !~ /GLOB/)
    {
        open (REPOUT, "<$results") or $self->throw("Error opening $results \n");
        $filehandle = \*REPOUT;
    }
    else
    {
        $filehandle = $self->results;
    } 
    
    
    #check if no repeats found
    if (<$filehandle> =~ /no repetitive sequences detected/)
    {
        print STDERR "RepeatMasker didn't find any repetitive sequences\n";
        close $filehandle;
        return;
    }
    #extract values
    
    while (<$filehandle>)
    {  
        if (/\d+/) #ignore introductory lines
        {
            my @element = split;
            # ignore features with negatives 
            next if ($element[11-13] =~ /-/); 
            my (%feat1, %feat2);
            $feat1 {name} = $element[4];
            $feat1 {score} = $element[0];
            $feat1 {start} = $element[5];
            $feat1 {end} = $element[6];
            #The start and end values are in different columns depending 
            #on orientation!
            if ($element[8] eq '+')
            {
                $feat2 {strand} = 1;
                $feat2 {start} = $element[11];     
                $feat2 {end} = $element[12];
            }
            elsif ($element[8] eq 'C')
            {
                $feat2 {strand} = -1 ;
                $feat2 {start} = $element[13];     
                $feat2 {end} = $element[12];
            }
            $feat2 {name} = $element[9];
            $feat2 {score} = $feat1 {score};
            $feat1 {strand} = $feat2 {strand};
            #misc
            $feat2 {db} = undef;
            $feat2 {db_version} = undef;
            $feat2 {program} = 'RepeatMasker';
            $feat2 {p_version}='1';
            $feat2 {source}= 'RepeatMasker';
            $feat2 {primary}= 'repeat';
            $feat1 {source}= 'RepeatMasker';
            $feat1 {primary}= 'repeat';
            $self->create_repeat(\%feat1, \%feat2);
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

1;
