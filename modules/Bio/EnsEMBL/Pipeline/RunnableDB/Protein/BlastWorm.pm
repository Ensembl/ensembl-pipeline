# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Protein::BlastWorm

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::Protein::BlastWorm->new ( -dbobj      => $db,
 	    	                                                          -input_id   => $input_id,
                                                                          -analysis   => $analysis,
                                                                        );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->output;
  $seg->write_output; # writes features to to DB

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Pipeline::Runnable::Protein::BlastWorm
  to add functionality to read and write to databases.
  A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor is required for database access (dbobj).
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Pipeline::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

  Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

  The rest of the documentation details each of the object methods. 
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::BlastWorm;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastWorm;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

 Title    : new
 Usage    : $self->new ( -DBOBJ       => $db
                         -INPUT_ID    => $id
                         -ANALYSIS    => $analysis,
                       );
 Function : creates a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::BlastWorm object
 Example  : 
 Returns  : a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::BlastWorm object
 Args     : -dbobj    :  a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
            -input_id :  input id
            -analysis :  a Bio::EnsEMBL::Pipeline::Analysis
 Throws   :

=cut

sub new {
    my ($class, @args) = @_;

    # this new method also parses the @args arguments,
    # and verifies that -dbobj and -input_id have been assigned
    my $self = $class->SUPER::new(@args);
    $self->throw ("Analysis object required") unless ($self->analysis);
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    # get the path to the binaries from the Analysis object (analysisprocess table)
    my $params .= "-program=>".$self->analysis->program_file.",";
    # get the analysisId from the Analysis object (analysisprocess table)
    $params .= "-analysisid=>".$self->analysis->dbID.",";
    # define the database
    $params .= "-database=>".$self->analysis->db_file.",";
    # set the filter
    $params .= "-filter=>1,";
    # define some threshold
    $params .= "-threshold=>0.1,";
    $params .= "-threshold_type=>PVALUE,";
    # define the options
    $params .= "-options=>".$self->analysis->parameters;
    
    $self->parameters($params);

    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::BlastWorm');
    return $self;
}

# IO methods

=head2 fetch_input

 Title    : fetch_input
 Usage    : $self->fetch_input
 Function : fetches the query sequence from the database
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub fetch_input {
    my ($self) = @_;
    my $protein_id = $self->input_id;
    my $proteinAdaptor = $self->dbobj->get_ProteinAdaptor;
    my $protein = $proteinAdaptor->fetch_Protein_by_dbid ($protein_id)
	|| $self->throw ("couldn't get the protein sequence from the database");
    $self->genseq($protein);
}


=head2 write_output

 Title    : write_output
 Usage    : $self->write_output
 Function : writes the features to the database
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub write_output {
    my ($self) = @_;
    my $proteinFeatureAdaptor = $self->dbobj->get_Protfeat_Adaptor;
    my @featurepairs = $self->output;
    foreach my $featurepair (@featurepairs) {
        $proteinFeatureAdaptor->write_Protein_Blast_feature ($featurepair);
    }
}

# runnable method

=head2 runnable

 Title    :  runnable
 Usage    :  $self->runnable($arg)
 Function :  sets a runnable for this RunnableDB
 Example  :
 Returns  :  Bio::EnsEMBL::Pipeline::RunnableI
 Args     :  Bio::EnsEMBL::Pipeline::RunnableI
 Throws   :

=cut

sub runnable {
    my ($self, $runnable) = @_;
    
    if ($runnable) {
        # extract parameters into a hash
        my ($parameter_string) = $self->parameters;
        my %parameters;
        if ($parameter_string) {
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs) {
                my ($key, $value) = split (/=>/, $pair);
		$key =~ s/\s+//g;
                # no, we need the spaces for the blast options
#                $value =~ s/\s+//g;
                $parameters{$key} = $value;
            }
        }
        $self->{'_runnable'} = $runnable->new (%parameters);
    }
    return $self->{'_runnable'};
}

1;
