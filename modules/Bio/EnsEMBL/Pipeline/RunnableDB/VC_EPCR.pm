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
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
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

    $self->dbobj->static_golden_path_type($sgp) if $sgp;

    my $vc = $self->dbobj->get_StaticGoldenPathAdaptor->
     fetch_VirtualContig_by_chr_start_end($chr, $start, $end);

    my $genseq = $vc->primary_seq() or $self->throw("Unable to fetch contig");
    $self->genseq($genseq);
    $self->vc($vc);
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
            $parameter_string =~ s/\s+//g;
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs)
            {
		my ($key, $value) = split (/=>/, $pair);
		if ($key && $value) {
		    $parameters{$key} = $value;
		}
		else
		{
		    $arguments .= " $key ";
		}
            }
        }
        $parameters {'-db'}      = $self->analysis->db_file();  
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

    my $db=$self->dbobj();
    my @features = $self->output();
    my %feat_by_contig;
  
    my $vc = $self->vc;

    foreach my $f (@features) {
	$f->analysis($self->analysis);
	$start  = $f->start;
	my ($raw, $raw_pos) = $vc->raw_contig_position($start);
	if ($raw && $raw_pos) {
	    if ($raw->static_golden_ori == 1) {
		$raw_start = $raw_pos;
		$raw_end = $f->end + $raw_start - $start;
	    }
	    else {
		$raw_end = $raw_pos;
		$raw_start = $start - $f->end + $raw_end;
	    }
	    my $raw_end = $f->end + $raw_start - $start;
	    $feat_by_contig{$raw->id} = [] unless defined $feat_by_contig{$raw->id};
	    $f->start($raw_start);
	    $f->end($raw_end);
	    push @{$feat_by_contig{$raw->id}}, $f;
	}
    }

    foreach my $contig_id (keys %feat_by_contig) {
	eval {
	    $contig = $db->get_Contig($contig_id);
	};
	print $contig, "\n";
	print length($@), "\n";

	if ($@) {
	    print STDERR "Contig not found, skipping writing output to db: $@\n";
	}
	elsif (@features = @{$feat_by_contig{$contig_id}}) {
	    print STDERR "Writing features to database\n";
	    my $feat_adp=Bio::EnsEMBL::DBSQL::FeatureAdaptor->new($db);
	    $feat_adp->store($contig, @features);
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
