#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Blast

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $blast   = Bio::EnsEMBL::Pipeline::RunnableDB::Blast->new ( 
                                                    -db         => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
$blast->fetch_input();
$blast->run();
$blast->output();
$blast->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Blast to add
functionality for reading and writing to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Blast;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Config::Blast;

use vars qw(@ISA);

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
   
    $self->throw("No input id") unless defined($self->input_id);

    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    my $genseq    = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");
    
    $self->query($genseq);
  
    my $seq = $self->query->seq;

    if ($seq =~ /[CATG]{3}/) {
        $self->input_is_void(0);
    } else {
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }

    my $ungapped;

    if($::pipeConf{'ungapped'}){
      $ungapped = 1;
    } else {
      $ungapped = undef;
    }

    my $run = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(-query          => $self->query,
							   -database       => $self->analysis->db_file,
							   -program        => $self->analysis->program,
							   -options        => $self->analysis->parameters,
							   -threshold_type => 'PVALUE',
							   -threshold      => 1,
							   -ungapped       => $ungapped,
							  );

    $self->runnable($run);

    return 1;
}

1;
