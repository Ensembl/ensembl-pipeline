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

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EPCR;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Map::MarkerFeature;
use Bio::EnsEMBL::Map::Marker;

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
           
    $self->{'_mflist'}    = [];    # an array of feature pairs
    $self->{'_query'}     = undef; # location of Bio::Seq object
    $self->{'_epcr'}      = undef; # location of EPCR binary
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of EPCR
    $self->{'_word_size'} = undef;     # e-PCR 'W' parameter
    $self->{'_margin'}    = undef;     # e-PCR 'M' parameter
    $self->{'_min_mismatch'}  = 0;     # e-PCR 'N' parameter
    $self->{'_max_mismatch'}  = 0;     # e-PCR 'N' parameter
    $self->{'_mismatch'}  = undef;     # e-PCR 'N' parameter
    $self->{'_protected'} = [];    # a list of files protected from deletion
    $self->{'_sts'}       = undef;
    $self->{'_hitlist'}   = {};
    
    my( $query, $epcr, $sts, $margin, $word_size, $min_mismatch,
    $max_mismatch) = $self->_rearrange([qw(
	QUERY EPCR STS M W NMIN NMAX
    )], @args);

    $epcr = 'e-PCR' unless ($epcr);

    $self->query  ($query)   if ($query);       
    $self->min_mismatch($min_mismatch)   if defined $min_mismatch;       
    $self->max_mismatch($max_mismatch)   if defined $max_mismatch;       
    $self->word_size($word_size) if defined $word_size;       
    $self->margin($margin)       if defined $margin;       

    $self->epcr   ($self->find_executable($epcr));

    $self->throw("No primer database or list of STSs entered\n") unless $sts;

    if (ref $sts) {
      $self->sts($sts);
    }
    else {
      my $file = $self->find_file($sts);
      $self->throw("STS file not found! Check path and name ($sts)\n") 
       unless -e $file;
      $self->sts($file);
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
        
        $self->filename($self->query->id.".$$.seq");
        $self->results($self->filename.".PCR.out");
        $self->file($self->results);
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
        $self->{'_mismatch'} = $args ;
    }
    return $self->{'_mismatch'};
}


sub sts {
    my ($self, $sts) = @_;
    if ($sts)
    {
        $self->{'_sts'} = $sts;
    }
    return $self->{'_sts'};
}


sub margin {
    my ($self, $args) = @_;
    if (defined $args)
    {
        $self->{'_margin'} = $args ;
    }
    return $self->{'_margin'};
}


=head2 run

Runs EPCR. Writes seqeuence to a temporary file, runs EPCR, parses
output file and deletes temporary files.

=cut

# copes with running with different values of N (primer bp mismatch)
# for different values of W still need a separate analysis

sub run {
    my ($self) = @_;

    #check seq
    my $seq = $self->query() || $self->throw("Seq required for EPCR\n");

    #set directory if provided
    $self->workdir('/tmp') unless $self->workdir();
    $self->checkdir();

    #write sequence to file
    $self->writefile();        

    # first run - M = M_MIN
    $self->mismatch($self->min_mismatch);

    if (ref $self->sts) {
        my $sts_tmp = $self->workdir . '/' . $self->query->id . "_0";
        $self->_dump_sts_file($sts_tmp);
        $self->run_epcr($sts_tmp);
        $self->_rm_sts_file($sts_tmp);
    }
    else {
        $self->run_epcr($self->sts);
    }
    $self->parse_results();

    # successive runs if needed
    for (my $mm = $self->min_mismatch + 1; $mm <= $self->max_mismatch; $mm++) {
	my $sts_tmp = $self->workdir . '/' . $self->query->id. "_$mm";
	if (ref $self->sts) {
	  $self->_dump_sts_file($sts_tmp);
	}
	else {
	  $self->_cp_sts_file($sts_tmp);
	}
	$self->mismatch($mm);
	$self->run_epcr($sts_tmp);
	$self->parse_results();
	$self->_rm_sts_file($sts_tmp);
    }

    # clean up
    $self->deletefiles();

    return 1;
}


sub run_epcr {
    my ($self, $file) = @_;
    #run EPCR

    my $mismatch  = $self->mismatch;
    my $margin    = $self->margin;
    my $word_size = $self->word_size;

    my $options = "";
    $options .= " M=$margin " if defined $margin;
    $options .= " W=$word_size " if defined $word_size;
    $options .= " N=$mismatch " if defined $mismatch;

    my $command = $self->epcr.' '.$file.' '.$self->filename.' '.
     $options.' > '.$self->results;

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
        my $feat;
	my ($name, $start, $end, $dbID) = $_ =~ m!(\S+)\s+(\d+)\.\.(\d+)\s+(\w+)!;

	$self->_add_hit($dbID);

        $feat->{'name'} = $name;
        $feat->{'start'} = $start;
        $feat->{'end'} = $end;
        $feat->{'strand'} = 0;
        $feat->{'score'} = -1000;
        $feat->{'dbID'} = $dbID;
        $feat->{'source'} = 'e-PCR';
        $feat->{'primary'} = 'sts';

	print "$name : $start : $end : $dbID\n";

	$self->create_MarkerFeature($feat);
    }
    close $filehandle;   
}


=head2 output

Returns an array of FeaturePair's (must be called after parse_results()).

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_mflist'}};
}

# make a temporary copy of the STS file to run against
# with those markers already hit removed

sub _cp_sts_file {
    my ($self, $dest) = @_;

    my $source = $self->sts;
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


sub _dump_sts_file {
    my ($self, $dest) = @_;

    my %hits = $self->_get_hits;
    eval {
        open DEST, "> $dest";
	foreach my $m (@{$self->sts}) {
	    next if $hits{$m->dbID};
	    next if $m->max_primer_dist == 0;
	    next unless length($m->left_primer) > 0;
	    next unless length($m->right_primer) > 0;
	    print DEST join("\t",
		$m->dbID,
		$m->left_primer,
		$m->right_primer,
		join("-", $m->min_primer_dist, $m->max_primer_dist),
	    ), "\n";
        }
        close DEST;
    };
    if ($@) {
	$self->throw("Unable to dump STS file");
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


sub create_MarkerFeature {
    my ($self, $feat) = @_;

    my $analysis = Bio::EnsEMBL::Analysis->new(
        -db              => undef,
        -db_version      => undef,
        -program         => $feat->{'program'},
        -program_version => $feat->{'program_version'},
        -gff_source      => $feat->{'source'},
        -gff_feature     => $feat->{'primary'}
    );

    my $m = Bio::EnsEMBL::Map::Marker->new;
    $m->dbID($feat->{'dbID'});

    my $mf = Bio::EnsEMBL::Map::MarkerFeature->new;
    $mf->analysis($analysis);
    $mf->score   ($feat->{'score'});
    $mf->seqname ($feat->{'name'});
    $mf->start   ($feat->{'start'});
    $mf->end     ($feat->{'end'});
    $mf->strand  ($feat->{'strand'});
    $mf->dbID    ($feat->{'dbID'});
    $mf->marker  ($m);

    # display_label must be a null string, and not undef
    # can't be set above as it is not known to SeqFeature
    # (SimpleFeature->new uses SeqFeature->new)
    # $mf->display_label($feat->{'hit'});

    if ($mf) {
	$mf->validate();

	# add to _mflist
	push(@{$self->{'_mflist'}}, $mf);
    }
}

1;
