# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::BlastWorm

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::BlastWorm->new ( -dbobj      => $db,
 	    	                                                 -input_id   => $input_id,
                                                                 -analysis   => $analysis,
                                                               );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->output;
  $seg->write_output; # writes features to to DB

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Pipeline::Runnable::BlastWorm
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

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastWorm;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastWorm;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

 Title    : new
 Usage    : $self->new ( -dbobj       => $db
                         -inpout_id    => $id
                         -analysis    => $analysis,
                       );
 Function : creates a Bio::EnsEMBL::Pipeline::RunnableDB::Blastp object
 Example  : 
 Returns  : a Bio::EnsEMBL::Pipeline::RunnableDB::Blastp object
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
    
    # set up seg specific parameters,
    # my $params = $self->parameters;
    # if ($params ne "") { $params .= ","; }
    # get the path to the binaries from the Analysis object (analysisprocess table)
    my $params .= "-program=>".$self->analysis->program_file.",";
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
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);
    my $genseq    = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");

    print STDERR "Setting genseq to " . $genseq. "\n";

    $self->genseq($genseq);
    print STDERR "Set genseq to " . $self->genseq. "\n";
    # input sequence needs to contain at least 3 consecutive nucleotides
    my $seq = $self->genseq->seq;
    $self->throw("Need at least 3 nucleotides") unless ($seq =~ /[CATG]{3}/);
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

    my @featurepairs = $self->output;

    unless (@featurepairs >= 1) {
        return;
    }

    # get the internal id of the contig
    my $sth = $self->dbobj->prepare ( q{ SELECT internal_id
                                           FROM contig
                                          WHERE id = ?
                                       } );
    $sth->execute ($featurepairs[0]->seqname);
    my $internalId = ($sth->fetchrow_array)[0];


    $sth = $self->dbobj->prepare ( q{ INSERT INTO feature
                                                  (id, contig, seq_start, seq_end,
                                                   score, strand, analysis, name,
                                                   hstart, hend, hid, evalue, perc_id, cigar)
                                           VALUES ('NULL', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                                    } );


    # get AnalysisAdaptor
    my $analysisAdaptor = $self->dbobj->get_AnalysisAdaptor;

    # write analysis to the database
    my $analysis = $featurepairs[0]->analysis;
    my $analysisId;

    unless ($analysis) {
        $self->throw ("Feature ".$featurepairs[0]->id ." doesn't have analysis. Cannot write to database");
    }

    unless ($analysisId = $analysisAdaptor->exists ($analysis)) {
        $analysisId = $analysisAdaptor->store ($analysis);
    }

    # loop over all featurepairs
    foreach my $featurepair (@featurepairs) {

        $featurepair->feature1->validate_prot_feature;
        $featurepair->feature2->validate_prot_feature;

        my $cigar_string;
        if ($featurepair->feature1->has_tag ('cigar')) {
            my @cigar_tags = $featurepair->feature1->each_tag_value ('cigar');
            $cigar_string = $cigar_tags[0];

        }

        ################
        my $target_seqname;

        # a little temporary hack to modify gadfly id's (Fly database)
        if ($featurepair->hseqname =~ /^\S+\|FB\S+\|(CT\d+)\|FB\S+/) {
            $target_seqname = "$1";
        }
        else {
            $target_seqname = $featurepair->hseqname;
        }

        ################

        $sth->execute ($internalId, $featurepair->start, $featurepair->end,
                       $featurepair->score, $featurepair->strand,
                       $analysisId, $self->analysis->program, 
#                       $featurepair->hstart, $featurepair->hend, $featurepair->hseqname,
                      $featurepair->hstart, $featurepair->hend, $target_seqname, 
                      $featurepair->p_value, $featurepair->percent_id, $cigar_string);
        
    }
    $sth->finish;
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
        $self->{'_runnable'} = $runnable->new (-analysis => $self->analysis, %parameters);
    }
    return $self->{'_runnable'};
}

1;

