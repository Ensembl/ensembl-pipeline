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
  my $pcr = Bio::EnsEMBL::Pipeline::Runnable::EPCR->new (-QUERY => $seq);
  $pcr->workdir($workdir);
  $pcr->run();
  my @results = $pcr->output();

=head1 DESCRIPTION

EPCR takes a Bio::Seq (or Bio::PrimarySeq) object and runs e-PCR on it. The
resulting .out file is parsed to produce a set of feature pairs.
Arguments can be passed to e-PCR through the arguments() method. 
Options can be passed to e-PCR through the options() method. 

=head2 Methods:

=over4

=item new($seq_obj)

=item epcr($path_to_EPCR)

=item workdir($directory_name)

=item arguments($args)

=item options($args)

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
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI );

=head2 new

    Title   :   new
    Usage   :   my $obj = Bio::EnsEMBL::Pipeline::Runnable::EPCR->new(-QUERY => $seq);
    Function:   Initialises EPCR object
    Returns :   a EPCR Object
    Args    :   A Bio::Seq object (-QUERY), a database (-DB).

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    $self->{'_fplist'}    = [];    # an array of feature pairs
    $self->{'_query'}     = undef; # location of Bio::Seq object
    $self->{'_epcr'}      = undef; # location of EPCR binary
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of EPCR
    $self->{'_options'}   = undef; # file to store results of EPCR
    $self->{'_protected'} = [];    # a list of files protected from deletion
    $self->{'_db'}        = undef;
    
    my( $query, $epcr, $db, $options) = $self->_rearrange([qw(
	QUERY PCR DB OPTIONS
    )], @args);

    $epcr = 'e-PCR' unless ($epcr);

    $self->query  ($query)   if ($query);       
    $self->options($options) if ($options);       

    $self->epcr   ($self->find_executable($epcr));

    if (defined($db)) {
      $self->db     ($self->find_file($db));
    } else {
      $self->throw("No primer database entered\n");
    }
    
    return $self; 
}

#################
# get/set methods 
#################

sub query {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
        }
        $self->{'_query'} = $seq ;
        
        $self->queryname($self->query->id);
        $self->filename($self->query->id.".$$.seq");
        $self->results($self->filename.".PCR.out");
    }
    return $self->{'_query'};
}

=head2 epcr

    Title   :   epcr
    Usage   :   $obj->epcr('~humpub/bin/EPCR');
    Function:   Get/set method for the location of EPCR
    Args    :   File path (optional)

=cut

sub epcr {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("e-PCR not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_epcr'} = $location ;
    }
    return $self->{'_epcr'};
}

=head2 options

    Title   :   options
    Usage   :   $obj->options('M=500');
    Function:   Get/set method for e-PCR options
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
        $self->{'_arguments'} = $args ;
    }
    return $self->{'_arguments'};
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

    Title   :   run
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs EPCR and creates array of featurepairs
    Returns :   none
    Args    :   optional $workdir and $args (e.g. '-ace' for ace file output)

=cut

sub run {
    my ($self, $dir, $args) = @_;

    #set arguments for epcr
    $self->arguments($args) if ($args);

    #check clone
    my $seq = $self->query() || $self->throw("Clone required for EPCR\n");

    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();

    #write sequence to file
    $self->writefile();        
    $self->run_epcr();

    #parse output of epcr
    $self->parse_results();
    $self->deletefiles();
}

sub run_epcr {
    my ($self) = @_;
    #run EPCR

    my $command = $self->epcr.' '.$self->db.' '.$self->filename.' '.
     $self->options.' > '.$self->results;

    print STDERR "Running EPCR ($command)\n";
    $self->throw("Error running EPCR on ".$self->filename."\n")
     if system($command); 
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
        open (PCROUT, "<".$self->results)
            or $self->throw("Error opening ".$self->results."\n");
        $filehandle = \*PCROUT;
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
	# OK? - this split will break if line doesn't start with whitespace
	# needs checking ...
        my @element = split (/\s+/, $output);  
        my (%feat1, %feat2);
        
        $feat1 {'name'}         = $element[0];
        $feat2 {'name'}         = ($element[3]) ? $element[3] : $element[2];
	# nasty hack - escape the "'" - sql barfs with things
	# like "3'UTR" - should be fixed in FeatureAdaptor...
	$feat2{name} =~ s{\'}{\\\'};
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
        
        $self->create_FeaturePair(\%feat1, \%feat2);
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
