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
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);



@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

my %UNGAPPED;
my %UNMASKED;

foreach my $db (@$DB_CONFIG) {
  my ($name, $ungapped, $unmasked) = ($db->{'name'}, $db->{'ungapped'}, $db->{min_unmasked});
  
  if($db && $name){
    $UNGAPPED{$name} = $ungapped;
    $UNMASKED{$name} = $unmasked;
  }else{
    my($p, $f, $l) = caller;
    warn("either db ".$db." or name ".$name." isn't defined so can't work $f:$l\n");
  }
}

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

    $self->fetch_sequence($PIPELINE_REPEAT_MASKING);
  
    my $seq = $self->query->seq;
    my $unmasked;
    if($UNMASKED{$self->analysis->db_file}){
      $unmasked = 1;
    } else {
      $unmasked = 3;
    }
    if ($seq =~ /[CATG]{$unmasked}/) {
        $self->input_is_void(0);
    } else {
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }

    my $ungapped;

    if($UNGAPPED{$self->analysis->db_file}){
      $ungapped = 1;
    } else {
      $ungapped = undef;
    }
    # note the parameters column of the analysis table for a blast
    # run must follow this format
    # -runnable_arg => value, -runnable_arg => value
    #for example -filter => 0, -options => -cpus 1 -spoutmax 1 
    #-hitdist 40 W=4 T=16 V=700000 B=700000 Y=320000000 Z=500000000
    my $run = Bio::EnsEMBL::Pipeline::Runnable::Blast->new
      (-query          => $self->query,
       -database       => $self->analysis->db_file,
       -program        => $self->analysis->program,
       $self->parameter_hash 
       -threshold_type => 'PVALUE',
       -threshold      => 1,
       -ungapped       => $ungapped,
      );
    
    $self->runnable($run);
    
    return 1;
}

1;
