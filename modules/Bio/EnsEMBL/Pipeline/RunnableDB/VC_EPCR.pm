#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# Modified 11.2001 by SCP to run on Virtual Contigs

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::VC_EPCR

=head1 SYNOPSIS

my $db   = Bio::EnsEMBL::DBLoader->new($locator);
my $epcr = Bio::EnsEMBL::Pipeline::RunnableDB::VC_EPCR->new( 
    -dbobj      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis
);
$epcr->fetch_input();
$epcr->run();
$epcr->output();
$epcr->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::EPCR to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::VC_EPCR;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::EPCR;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::VC_EPCR object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::VC_EPCR object
    Args    :   -dbobj:     A Bio::EnsEMBL::DBSQL::DBAdaptor, 
                -input_id:  A virtual contig (e.g. chr_name.start.end)
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    $self->{'_parameters'}  = undef;
    
    $self->throw("Analysis object required") unless ($self->analysis);
    
    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::EPCR');
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for epcr from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $vc_str = $self->input_id;
    my ($chr, $start, $end, $sgp) = $vc_str =~ m!(\S+)\.(\d+)\.(\d+):?([^:]*)!;

    $self->db->assembly_type($sgp) if $sgp;

    my $slice = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr, $start, $end);
    $self->slice($slice);

    my $genseq = $slice->primary_seq() or $self->throw("Unable to fetch contig");
    $self->genseq($genseq);
}

#get/set for runnable and args
sub runnable {
    my ($self, $runnable) = @_;
    my $arguments = "";

    # $self->analysis->parameters is a comma-delimited list.
    # Anything of the form a => b is passed to the runnable's new method.
    # Other text is given to the runnable as "options"
    
    if ($runnable)
    {
        #extract parameters into a hash
        my ($parameter_string) = $self->analysis->parameters() ;
        my %parameters;
        if ($parameter_string)
        {
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs)
            {
		my ($key, $value) = split (/=>/, $pair);
		if ($key && $value) {
                    $key   =~ s/^\s+//g;
                    $key   =~ s/\s+$//g;
                    $value =~ s/^\s+//g;
                    $value =~ s/\s+$//g;
		    $parameters{$key} = $value;
		}
		else
		{
		    $arguments .= " $key ";
		}
            }
        }
        $parameters {'-sts'}     = $self->analysis->db_file();  
        $parameters {'-options'} = $arguments;
        $parameters {'-pcr'}     = $self->analysis->program_file();  
        #creates empty Bio::EnsEMBL::Runnable::EPCR object
        $self->{'_runnable'} = $runnable->new(%parameters);
    }
    return $self->{'_runnable'};
}


sub write_output {
    my($self) = @_;
    my $start;
    my $contig;
    my ($raw_start, $raw_end);
    my ($hit_start, $hit_end);
    my $max_walk = 15;   # max bp to walk until start of raw contig

    my $db=$self->db();
    my $simple_f_a = $self->db->get_SimpleFeatureAdaptor();
    
    my @features = $self->output();
    my %feat_by_contig;
  
    my $vc = $self->slice;

    foreach my $f (@features) {
	$f->analysis($self->analysis);
	$hit_start = $f->start;
	$hit_end   = $f->end;
	my ($raw1, $raw1_pos);  # contig with hit start
	my ($raw2, $raw2_pos);  # contig with hit end
	my $raw;

	# if a primer has N's at the 5' end, the start of the STS as
	# reported may be in a gap region. need to walk along until
	# a raw contig is found.

	# get raw ctg/pos correspondong to start of STS
	# walk along until we can map the STS start to a raw contig
	my $walkies = 0;
	while ($walkies < $max_walk and ! ref $raw1) {
	    ($raw1, $raw1_pos) = $vc->raw_contig_position($hit_start);
	    $hit_start++;
	}
	$hit_start--;

	$walkies = 0;
	# get raw ctg/pos correspondong to end of STS
	# walk along until we can map the STS end to a raw contig
	while ($walkies < $max_walk and ! ref $raw2) {
	    ($raw2, $raw2_pos) = $vc->raw_contig_position($hit_end);
	    $hit_end--;
	}
	$hit_end++;

	# exclude 'sticky' markers - spanning two raw contigs
	# we should deal with this somehow ...
	unless ($raw1->id eq $raw2->id) {
	    # yikes - sticky marker
	    print STDERR "Ignoring marker spanning two raw contigs: ";
	    print STDERR $raw1->id, " ", $raw2->id, "\n";
	    next;
	}
	$raw       = $raw1;
	$raw_start = $raw1_pos;
	$raw_end   = $raw2_pos;

	if ($raw->static_golden_ori != 1) {
	    ($raw_start, $raw_end) = ($raw_end, $raw_start);
	}

	$feat_by_contig{$raw->id} = [] unless defined $feat_by_contig{$raw->id};
	$f->start($raw_start);
	$f->end($raw_end);
	push @{$feat_by_contig{$raw->id}}, $f;
    }

    foreach my $contig_id (keys %feat_by_contig) {
	eval {
	    $contig = $db->get_RawContigAdaptor->fetch_by_name($contig_id);
	};

	if ($@) {
	    print STDERR "Contig not found, skipping writing output to db: $@\n";
	}
	elsif (@features = @{$feat_by_contig{$contig_id}}) {
	  foreach my $f(@features){
	    $f->analysis($self->analysis);
	    $simple_f_a->store($contig->dbID, $f);
	  }
	}

    }
    return 1;
}


=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetches output data from a frozen perl object
                stored in file $file_name
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut


1;
