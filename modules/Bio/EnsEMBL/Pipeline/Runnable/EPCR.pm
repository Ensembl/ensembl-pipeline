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
  my $seqfile = '/nfs/disk65/mq2/temp/bA151E14.seq';
  my $seqstream = Bio::SeqIO->new(-file => $seqfile, -fmt => 'Fasta');
  my $seq = $seqstream->next_seq();

  #create and run Bio::EnsEMBL::Pipeline::Runnable::EPCR object
  my $pcr = Bio::EnsEMBL::Pipeline::Runnable::EPCR->new (-QUERY => $seq);
  $pcr->workdir($workdir);
  $pcr->run();

  # collect results
  my @results = $pcr->output();

=head1 DESCRIPTION

EPCR takes a Bio::Seq (or Bio::PrimarySeq) object and runs e-PCR on it. The
resulting .out file is parsed to produce a set of feature pairs.
Options can be passed to e-PCR through the options() method. 

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EPCR;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

Initializes EPCR object. Takes a Bio::Seq object and an STS database.

  my $epcr = Bio::EnsEMBL::Pipeline::Runnable::EPCR->new(
      -query => $seq
      -sts   => 'sts_file'
  );

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
    $self->{'_sts'}       = undef;
    
    my( $query, $epcr, $sts, $options) = $self->_rearrange([qw(
	QUERY PCR STS OPTIONS
    )], @args);

    $epcr = 'e-PCR' unless ($epcr);

    $self->query  ($query)   if ($query);       
    $self->options($options) if ($options);       

    $self->epcr   ($self->find_executable($epcr));

    if (defined($sts)) {
      $self->sts     ($self->find_file($sts));
    } else {
      $self->throw("No primer database entered\n");
    }
    
    return $self; 
}


=head2 Get/set methods

=over 4

=item epcr

name of binary (e.g. 'e-PCR'),
or explicit location if a full path is given

=item options

options to pass to e-PCR, such as margin size 'M=200'

=item query

a Bio::Seq like object on which the analysis will be run

=item sts

Location of STS file containing primer pairs and names of markers
(e.g. 'mapprimer_maps_110')

=back

=cut


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


sub epcr {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("e-PCR binary not found at $location: $!\n") 
         unless (-e $location);
        $self->{'_epcr'} = $location ;
    }
    return $self->{'_epcr'};
}


sub sts {
    my ($self, $sts) = @_;
    if ($sts)
    {
        $self->throw("STS file not found! Check path and name ($sts)\n") 
         unless (-e $sts);
        $self->{'_sts'} = $sts;
    }
    return $self->{'_sts'};
}


=head2 run

Runs EPCR. Writes seqeuence to a temporary file, runs EPCR, parses
output file and deletes temporary files.

=cut

sub run {
    my ($self) = @_;

    #check seq
    my $seq = $self->query() || $self->throw("Seq required for EPCR\n");

    #set directory if provided
    $self->workdir('/tmp') unless $self->workdir();
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

    my $command = $self->epcr.' '.$self->sts.' '.$self->filename.' '.
     $self->options.' > '.$self->results;

    print STDERR "Running EPCR ($command)\n";
    $self->throw("Error running EPCR on ".$self->filename."\n")
     if system($command); 
    #or $self->throw("Error running EPCR: $!\n")
    #check results exist
    $self->throw($self->results." not created by EPCR\n") unless (-e $self->results);   
}

=head2 parse_results

Takes either an ouput file or a file handle and parses output.
Creates FeaturePair's;

=cut

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
    
    foreach (@output) #loop from 3rd line
    {  
        my @element = split;
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
        $feat2 {'end'}          = $end - $start + 1;
        $feat1 {'strand'}       = 0;
        $feat2 {'strand'}       = 0;
        #misc
        $feat1 {'score'}        = -1000;
        $feat2 {'score'}        = -1000;
        $feat2 {'db'}           = $self->sts || undef;
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


=head2 output

Returns an array of FeaturePair's (must be called after parse_results()).

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_fplist'}};
}

1;
