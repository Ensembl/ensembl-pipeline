#
#
# Cared for by Ensembl <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::VC_MarkovRegions

=head1 SYNOPSIS

my $db     = Bio::EnsEMBL::DBLoader->new($locator);

my $markov = Bio::EnsEMBL::Pipeline::RunnableDB::VC_MarkovRegions->new ( 
                                    -dbobj      => $db,
			            -input_id   => $input_id
                                    -analysis   => $analysis 
                                    );

$markov->fetch_input();

$markov->run();

$markov->write_output(); # writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::MarkovRegions.pm
to add functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Pipeline::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::VC_MarkovRegions;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::MarkovRegions;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a
                Bio::EnsEMBL::Pipeline::RunnableDB::VC_MarkovRegions
		object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::VC_MarkovRegions
                object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:  A virtual contig (e.g. 'chr_name.start.end')
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}      = []; # ???   
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    $self->throw("Analysis object required") unless ($self->analysis);
    
    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::MarkovRegions');
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for MarkovRegions from the
                database
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
    if ($runnable)
    {
        #extract parameters into a hash
        my ($parameter_string) = $self->parameters() ;
        my %parameters;
        if ($parameter_string)
        {
            $parameter_string =~ s/\s+//g;
            my @pairs = split (/,/, $parameter_string);
            
            foreach my $pair (@pairs)
            {
                my ($key, $value) = split (/=>/, $pair);
                $parameters{$key} = $value;
            }
        }
        #creates empty Bio::EnsEMBL::Runnable::MarkovRegions object
        $self->{'_runnable'} = $runnable->new
	    ( '-threshold'  => $parameters{'-threshold'},
	      '-iterations' => $parameters{'-iterations'},
	      '-motif_file' => $parameters{'-motif_file'},
	      );
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


1;
