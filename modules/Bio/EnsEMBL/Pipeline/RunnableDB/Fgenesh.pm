
#
#
# Cared for by Laura Clarke  <lec@sanger.ac.uk>
#
# Copyright Laura Clarke

#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB:Fgenesh
=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $fgenesh     = Bio::EnsEMBL::Pipeline::RunnableDB::Fgenesh->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$fgenesh->fetch_input();
$fgenesh->run();
$fgenesh->output();
$fgenesh->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Fgenesh to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Fgenesh;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Fgenesh;
use Bio::EnsEMBL::Pipeline::Config::General;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub fetch_input {
    my( $self) = @_;
    

    $self->throw("No input id") unless defined($self->input_id);
    
    $self->fetch_sequence;
    
    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Fgenesh
      (-query   => $self->query,
       -fgenesh => $self->analysis->program_file,
       -matrix  => $self->analysis->db_file,
       -param   => $self->arguments
      );
    $self->runnable($runnable);

    return 1;
}


1;






