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
Options are passed to e-PCR through M, W, NMIN, NMAX.

=head2 Methods:

=over4

=item new($seq_obj)

=item epcr($path_to_EPCR)

=item workdir($directory_name)

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
    Usage   :   my $obj = Bio::EnsEMBL::Pipeline::Runnable::EPCR->new(-CLONE => $seq);
    Function:   Initialises EPCR object
    Returns :   a EPCR Object
    Args    :   A Bio::Seq object (-CLONE), a database (-DB).

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    $self->{'_fplist'}    = [];    # an array of feature pairs
    $self->{'_clone'}     = undef; # location of Bio::Seq object
    $self->{'_epcr'}      = undef; # location of EPCR binary
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of EPCR
    $self->{'_word_size'} = undef;     # e-PCR 'W' parameter
    $self->{'_margin'}    = undef;     # e-PCR 'M' parameter
    $self->{'_min_mismatch'}  = undef;     # e-PCR 'N' parameter
    $self->{'_max_mismatch'}  = undef;     # e-PCR 'N' parameter
    $self->{'mismatch'}   = undef;     # e-PCR 'N' parameter
    $self->{'_protected'} = [];    # a list of files protected from deletion
    $self->{'_db'}        = undef;
    $self->{'_hitlist'}   = {};
    
    my( $clone, $epcr, $db, $margin, $word_size, $min_mismatch,
    $max_mismatch) = $self->_rearrange([qw(
	CLONE PCR DB M W NMIN NMAX
    )], @args);

    $epcr = 'e-PCR' unless ($epcr);

    $self->clone  ($clone)       if ($clone);       
    $self->min_mismatch($min_mismatch)   if defined $min_mismatch;       
    $self->max_mismatch($max_mismatch)   if defined $max_mismatch;       
    $self->word_size($word_size) if defined $word_size;       
    $self->margin($margin)       if defined $margin;       

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
        $self->results($self->filename.".PCR.out");
    }
    return $self->{'_clone'};
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


sub word_size {
    my ($self, $args) = @_;
    if (defined $args)
    {
        $self->{'_word_size'} = $args ;
    }
    return $self->{'_word_size'};
}

sub min_mismatch {
    my ($self, $args) = @_;
    if (defined $args)
    {
        $self->{'_min_mismatch'} = $args ;
    }
    return $self->{'_min_mismatch'};
}

sub max_mismatch {
    my ($self, $args) = @_;
    if (defined $args)
    {
        $self->{'_max_mismatch'} = $args ;
    }
    return $self->{'_max_mismatch'};
}

sub mismatch {
    my ($self, $args) = @_;
    if (defined $args)
    {
        $self->{'mismatch'} = $args ;
    }
    return $self->{'mismatch'};
}

sub margin {
    my ($self, $args) = @_;
    if (defined $args)
    {
        $self->{'_margin'} = $args ;
    }
    return $self->{'_margin'};
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

# copes with running with different values of N (primer bp mismatch)
# for different values of W still need a separate analysis

sub run {
    my ($self, $dir) = @_;

    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for EPCR\n");

    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();

    #write sequence to file
    $self->writefile();        

    # first run - M = M_MIN
    $self->mismatch($self->min_mismatch);
    $self->run_epcr();
    $self->parse_results();

    # successive runs if needed
    for (my $mm = $self->min_mismatch + 1; $mm <= $self->max_mismatch; $mm++) {
	my $sts_tmp = $self->workdir . '/' . $self->clonename . "_$mm";
	$self->_cp_sts_file($self->db, $sts_tmp);
	$self->mismatch($mm);
	my $db = $self->db;
	$self->db($sts_tmp);
	$self->run_epcr;
	$self->db($db);  # reset to original STS file
			 # bit inefficient, better to work from
			 # current STS file, but this is simpler
			 # and STS files are normally huge
	$self->parse_results();
	$self->_rm_sts_file($sts_tmp);
    }

    # clean up
    $self->deletefiles();
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses EPCR output to give a set of feature pairs
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

sub run_epcr {
    my ($self) = @_;
    #run EPCR

    my $mismatch  = $self->mismatch;
    my $margin    = $self->margin;
    my $word_size = $self->word_size;

    my $options = "";
    $options .= " M=$margin " if defined $margin;
    $options .= " W=$word_size " if defined $word_size;
    $options .= " N=$mismatch " if defined $mismatch;

    my $command = $self->epcr.' '.$self->db.' '.$self->filename.' '.
     $options.' > '.$self->results;

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
	$self->_add_hit($feat2{name});
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

# make a temporary copy of the STS file to run against
# with those markers already hit removed

sub _cp_sts_file {
    my ($self, $source, $dest) = @_;

    my %hits = $self->_get_hits;
    eval {
        open SOURCE, "< $source";
        open DEST, "> $dest";
        while (<SOURCE>) {
	    my $id = (split)[0];
	    print DEST unless defined $hits{$id};
        }
        close DEST;
        close SOURCE;
    };
    if ($@) {
	$self->throw("Unable to copy STS file");
    }
}

sub _add_hit {
    my ($self, $hit) = @_;

    $self->{'_hitlist'}->{$hit} = 1;
}


# return a hash of hit IDs

sub _get_hits {
    my ($self) = @_;

    return %{$self->{'_hitlist'}};
}

sub _rm_sts_file {
    my ($self, $file) = @_;

    unlink $file;
}

1;
