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

Bio::EnsEMBL::Pipeline::Runnable::Scan

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq';
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::Scan object
  my $scan = Bio::EnsEMBL::Pipeline::Runnable::Scan->new (-CLONE => $seq);
  $scan->workdir($workdir);
  $scan->run();
  my @results = $scan->output();

=head1 DESCRIPTION

Scan takes a Bio::Seq (or Bio::PrimarySeq) object and runs Scan on it. The
resulting .out file is parsed to produce a set of feature pairs.
Arguments can be passed to Scan through the arguments() method. 

=head2 Methods:

=over4

=item new($seq_obj)

=item Scan($path_to_Scan)

=item workdir($directory_name)

=item arguments($args)

=item run()

=item output()

=back

=head1 CONTACT

Kerstin Jekosch <kj2@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Scan;

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
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Scan->new (-CLONE => $seq);
    Function:   Initialises Scan object
    Returns :   a Scan Object
    Args    :   A Bio::Seq object (-CLONE), any arguments for Scan (-ARGS) 

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    $self->{'_fplist'}    = [];        #an array of feature pairs
    $self->{'_clone'}     = undef;     #location of Bio::Seq object
    $self->{'_program'}   = undef;     #location of Scan executable
    $self->{'_workdir'}   = undef;     #location of temp directory
    $self->{'_filename'}  = undef;     #file to store Bio::Seq object
    $self->{'_results'}   = undef;     #file to store results of Scan
    $self->{'_protected'} = [];        #a list of files protected from deletion
    $self->{'_arguments'} = undef;     #arguments for Scan
    
    my( $clone, $arguments, $scan) = $self->_rearrange([qw(CLONE
							   ARGS
							   SCAN)], 
							  @args);

    $scan = '/nfs/farm/Worms/bin/scan_lib.pl' unless defined($scan);

    $self->clone($clone) if ($clone);       

#    $self->program ($self->find_executable ($self->analysis->program_file));
    my $scan = $self->find_executable($scan);

    $self->Scan($scan);

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

=head2 Scan

    Title   :   Scan
    Usage   :   $obj->Scan('~humpub/scripts/Scan');
    Function:   Get/set method for the location of Scan script
    Args    :   File path (optional)

=cut

sub Scan {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Scan not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_Scan} = $location ;
    }
    return $self->{_Scan};
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
    Function:   Runs Scan script and creates array of featurepairs
    Returns :   none
    Args    :   optional $workdir and $args (e.g. '-ace' for ace file output)

=cut

sub run {
    my ($self, $dir, $args) = @_;

    #set arguments for Scan
    $self->arguments($args) if ($args);

    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Scan\n");

    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();

    # reset filename and results as necessary (adding the directory path)
    my $tmp = $self->workdir;
    my $input = $tmp."/".$self->filename;
    $self->filename ($input);
    $tmp .= "/".$self->results;
    $self->results ($tmp);

    #write sequence to file
    $self->writefile();        
    $self->run_Scan();

    #parse output of repeat masker 
    $self->parse_results();
    $self->deletefiles();
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses Scan output to give a set of feature pairs
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub run_Scan {
    my ($self) = @_;
    #run Scan
    print "Running Scan; ";
#    $self->throw('line is: '.$self->Scan.' '.$self->arguments.' '.$self->filename.' > '.$self->results);
    print "cmd: ", $self->Scan.' '.$self->arguments.' '.$self->filename.' > '.$self->results, "\n";
    $self->throw("Error running Scan on ".$self->filename."\n") 
        if (system ($self->Scan.' '.$self->arguments.' '.$self->filename.' > '.$self->results)); 
}

#New and improved! takes filenames and handles, therefore pipe compliant!
sub parse_results {
    my ($self) = @_;
    my $filehandle;
    my $results = (-e $self->results) ? $self->results : $self->filename.".out";
    if (ref ($self->results) !~ /GLOB/) {
        open (REPOUT, "<$results") or $self->throw("Error opening $results \n");
        $filehandle = \*REPOUT;
    }
    else {
        $filehandle = $self->results;
    } 
    
    
    #check if no repeats found
#    if (<$filehandle> =~ /^$/) {
#	$self->throw("haha:$_");
#        print STDERR "Scan didn't find any repetitive sequences\n";
#        close $filehandle;
#        return;
#    }
#    else {
#    	$self->throw("hehe:$_");
#    }
    #extract values
    
    while (<$filehandle>)
    {  
        if (/\S+/) #ignore introductory lines
        {
            my @element = split;
            # ignore features with negatives 
#            next if ($element[11-13] =~ /-/); 
            my (%feat1, %feat2);
            $feat1 {name}  = $element[2];
            $feat1 {score} = $element[1];
            if ($element[3] < $element[4]) {
		    $feat1 {start}  = $element[3];
        	    $feat1 {end}    = $element[4];
	            $feat2 {strand} = 1;
            }
	    elsif ($element[4] < $element[3]) {
		    $feat1 {start}  = $element[4];
        	    $feat1 {end}    = $element[3];
   	            $feat2 {strand} = -1;
            }
	    $feat2 {start}  = $element[5];     
	    $feat2 {end}    = $element[6];
            $feat2 {name}   = $element[0];
            $feat2 {score}  = $feat1 {score};
            $feat1 {strand} = $feat2 {strand};
            #misc
            $feat2 {db}         = undef;
            $feat2 {db_version} = undef;
            $feat2 {program}    = 'Scan';
            $feat2 {p_version}  = '1';
            $feat2 {source}     = 'Scan';
            $feat2 {primary}    = 'repeat';
            $feat1 {source}     = 'Scan';
            $feat1 {primary}    = 'repeat';
#	    $self->throw("$_\n	clone $element[2],score $element[1],start $element[3],stop $element[4],score $feat2{score},strand $feat1{strand}\n");
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
