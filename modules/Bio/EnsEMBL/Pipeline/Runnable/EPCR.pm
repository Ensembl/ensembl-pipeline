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

Bio::EnsEMBL::Pipeline::Runnable::EPCR

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq';
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::EPCR object
  my $pcr = Bio::EnsEMBL::Pipeline::Runnable::EPCR->new (-CLONE => $seq);
  $pcr->workdir($workdir);
  $pcr->run();
  my @results = $pcr->output();

=head1 DESCRIPTION

EPCR takes a Bio::Seq (or Bio::PrimarySeq) object and runs e-PCR on it. The
resulting .out file is parsed to produce a set of feature pairs.
Arguments can be passed to e-PCR through the arguments() method. 

=head2 Methods:

=over4

=item new($seq_obj)

=item repeatmasker($path_to_EPCRHum)

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

package Bio::EnsEMBL::Pipeline::Runnable::EPCR;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::Object;
#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object );

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::EPCR->new (-CLONE => $seq);
    Function:   Initialises EPCR object
    Returns :   a EPCR Object
    Args    :   A Bio::Seq object (-CLONE), a database (-DB).

=cut

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{_fplist} = [];           #an array of feature pairs
    $self->{_clone}  = undef;        #location of Bio::Seq object
    $self->{_epcr} = undef;  #location of EPCRHum script
    $self->{_workdir}   = undef;     #location of temp directory
    $self->{_filename}  =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;      #file to store results of EPCRHum
    $self->{_protected} =[];         #a list of files protected from deletion
    $self->{_db}        =undef;
    
    my( $clone, $epcr, $db) 
        = $self->_rearrange(['CLONE', 'PCR', 'DB'], @args);
    
    $self->clone($clone) if ($clone);       
    if ($epcr)
    {   $self->epcr($epcr);  }
    else
    {   
        eval 
        { $self->epcr($self->locate_executable('e-PCR')); };
        if ($@)
        { $self->epcr('/usr/local/badger/bin/e-PCR'); }
    }
    
    $self->throw("sts or primer Database required ($db)\n") unless ($db);
    $self->db($db);
    
    return $self; # success - we hope!
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
        $self->{_clone} = $seq ;
        
        $self->clonename($self->clone->id);
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".PCR.out");
    }
    return $self->{_clone};
}

=head2 epcr

    Title   :   epcr
    Usage   :   $obj->epcr('~humpub/scripts/EPCR');
    Function:   Get/set method for the location of EPCR
    Args    :   File path (optional)

=cut

sub epcr {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("e-PCR not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_epcr} = $location ;
    }
    return $self->{_epcr};
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

sub db {
    my ($self, $db) = @_;
    if ($db)
    {
        $self->throw("DB not found! Check path and name ($db)\n") 
                                                    unless (-e $db);
        $self->{'_db'} = $db;
    }
    return $self->{'_db'};
}

###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs EPCRHum script and creates array of featurepairs
    Returns :   none
    Args    :   optional $workdir and $args (e.g. '-ace' for ace file output)

=cut

sub run {
    my ($self, $dir, $args) = @_;
    #set arguments for repeatmasker
    $self->arguments($args) if ($args);
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for EPCR\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile();        
    $self->run_epcr();
    #parse output of repeat masker 
    $self->parse_results();
    $self->deletefiles();
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses EPCRHum output to give a set of feature pairs
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub run_epcr {
    my ($self) = @_;
    #run EPCRHum
    my $command = $self->epcr.' '.$self->db.' '.$self->filename.' > '.$self->results;
    print STDERR "Running EPCR ($command)\n";
    system ($command); 
    #or $self->throw("Error running EPCR: $!\n")
    #check results exist
    $self->throw($self->results." not created by EPCR\n") unless (-e $self->results);   
}

#New and improved! takes filenames and handles, therefore pipe compliant!
sub parse_results {
    my ($self) = @_;
    my $filehandle;
    if (ref ($self->results) !~ /GLOB/)
    {
        open (REPOUT, "<".$self->results)
            or $self->throw("Error opening ".$self->results."\n");
        $filehandle = \*REPOUT;
    }
    else
    {
        $filehandle = $self->results;
    } 
    #print STDERR <$filehandle>;    
    my @output = <$filehandle>;
    
    if (scalar (@output) == 0)
    {
        print STDERR "e-PCR didn't find any hits on this sequence\n";
        close $filehandle;
        return;
    }
    
    foreach my $output (@output) #loop from 3rd line
    {  
        my @element = split (/\s+/, $output);  
        my (%feat1, %feat2);
        
        $feat1 {'name'}         = $element[0];
        $feat2 {'name'}         = ($element[3]) ? $element[3] : $element[2];
        my ($start, $end)       = split (/\.\./, $element[1]);
        $feat1 {'start'}        = $start;
        $feat1 {'end'}          = $end;
        $feat2 {'start'}        = 1;     
        $feat2 {'end'}          = $end - $start;
        $feat1 {'strand'}       = 0;
        $feat2 {'strand'}       = 0;
        #misc
        $feat1 {'score'}        = -1000;
        $feat2 {'score'}        = -1000;
        $feat2 {'db'}           = $self->db || undef;
        $feat2 {'db_version'}   = 1;
        $feat2 {'program'}      = 'e-PCR';
        $feat2 {'p_version'}    ='1';
        $feat2 {'source'}       = 'e-PCR';
        $feat2 {'primary'}      = 'sts';
        $feat1 {'source'}       = 'e-PCR';
        $feat1 {'primary'}      = 'sts';
        
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

