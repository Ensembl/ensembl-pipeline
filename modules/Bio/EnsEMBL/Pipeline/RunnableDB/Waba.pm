#
#
# Cared for by Marc Sohrmann  <ms2@sanger.ac.uk>
#
# Copyright Marc Sohrmann
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

 Bio::EnsEMBL::Pipeline::RunnableDB::Waba

=head1 SYNOPSIS

 my $trna = Bio::EnsEMBL::Pipeline::RunnableDB::Waba->new ( -dbobj      => $db,
			                                    -input_id   => $input_id
                                                            -analysis   => $analysis );
 $trna->fetch_input();
 $trna->run();
 $trna->output();
 $trna->write_output(); #writes to DB

=head1 DESCRIPTION

 This object wraps Bio::EnsEMBL::Pipeline::Runnable::Waba to add
 functionality to read and write to databases.
 The appropriate Bio::EnsEMBL::Analysis object must be passed for
 extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
 required for database access.

=head1 CONTACT

 Describe contact details here

=head1 APPENDIX

 The rest of the documentation details each of the object methods. 
 Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Waba;

use strict;
$|=1;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Waba;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Waba object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Waba object
    Args    :   -dbobj:     A Bio::EnsEMBL::DBSQL::DBAdaptor, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}      = [];  
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    

    $self->throw("Analysis object required") unless ($self->analysis);
    
    # set up waba specific parameters (none for the moment)
    my $params = $self->parameters();
    if ($params ne "") { $params .= " , "; }

    $self->parameters($params);

    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::Waba');
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Waba from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);
    
    my $ctg_name  = $self->input_id;
    my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($ctg_name);
    #put the contig in to the runnable
    $self->runnable->query($contig);

    $self->query($contig);
}

=head2 runnable

    Title   :   runnable
    Usage   :   $self->runnable($arg)
    Function:   Sets a runnable for this RunnableDB
    Returns :   Bio::EnsEMBL::Pipeline::RunnableI
    Args    :   Bio::EnsEMBL::Pipeline::RunnableI

=cut


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

    my @refs  = $self->output;

    unless (@refs >= 1) {
        return;
    }

    # prepare sql
    my $sth_c = $self->db->prepare ( q{ SELECT contig_id
                                             FROM contig
                                            WHERE name = ?
				         } );

    my $sth = $self->db->prepare ( q{ INSERT INTO waba_feature
                                                  (id, contig_id, seq_start, seq_end,
                                                   score, strand, analysis_id, name,
                                                   hit_start, hit_end, hit_name, perc_ident, state, cigar)
					   VALUES ('NULL', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
				       } );

    my $sth_fset = $self->db->prepare ( q{ INSERT INTO waba_fset
                                                          (id, contig_id, seq_start, seq_end,
                                                           score, strand, analysis_id, name,
                                                           hit_start, hit_end, hit_name, perc_ident)
                                                   VALUES ('NULL', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
					    } );

    my $sth_fset_feature = $self->db->prepare ( q{ INSERT INTO waba_fset_feature
                                                                  (feature, fset, rank)
                                                           VALUES (?, ?, ?)
					            } );

    my $sth_last = $self->db->prepare ( q{ SELECT last_insert_id() } );

    # loop over all arrays
    foreach my $aref (@refs) {
        my @fps = @$aref;
        my $fset_fp = shift @fps;

      #my $fset_fp = $aref;

        # get the internal id of the contig
        $sth_c->execute ($fset_fp->seqname);
        my $internalId = ($sth_c->fetchrow_array)[0];

        # get AnalysisAdaptor
        my $analysisAdaptor = $self->db->get_AnalysisAdaptor;

        # write analysis to the database
        my $analysis = $fset_fp->analysis;
        my $analysisId;
        unless ($analysis) {
            $self->throw ("Feature ".$fset_fp->id ." doesn't have analysis. Cannot write to database");
        }
        unless ($analysisId = $analysisAdaptor->exists ($analysis)) {
            $analysisId = $analysisAdaptor->store ($analysis);
        }

        # write waba_fset (total match) to the db
        $sth_fset->execute ($internalId, $fset_fp->start, $fset_fp->end,
                            $fset_fp->score, $fset_fp->strand,
                            $analysisId, $self->analysis->program, 
                            $fset_fp->hstart, $fset_fp->hend, $fset_fp->hseqname,
                            $fset_fp->percent_id);

        $sth_last->execute;
        my $fset_id = ($sth_last->fetchrow_array)[0];

        # write the features to the db
        my $rank = 0;
        foreach my $featurepair (@fps) {
            $featurepair->feature1->validate_prot_feature;
            $featurepair->feature2->validate_prot_feature;
  
            my $state;
            my $cigar;
            if ($featurepair->feature1->has_tag ('state')) {
                my @state_tags = $featurepair->feature1->each_tag_value ('state');
                $state = $state_tags[0];
	    }
            else {
                $self->throw ("Waba feature needs the state info");
	    }
            if ($featurepair->feature1->has_tag ('cigar')) {
                my @cigar_tags = $featurepair->feature1->each_tag_value ('cigar');
                $cigar = $cigar_tags[0];
	    }
            $sth->execute ($internalId, $featurepair->start, $featurepair->end,
                           $featurepair->score, $featurepair->strand,
                           $analysisId, $self->analysis->program, 
                           $featurepair->hstart, $featurepair->hend, $featurepair->hseqname,
                           $featurepair->percent_id, $state, $cigar);

            $sth_last->execute;
            my $feature_id = ($sth_last->fetchrow_array)[0];

            $sth_fset_feature->execute ($feature_id, $fset_id, ++$rank); 
        }
    }
    $sth->finish;
    $sth_c->finish;
    $sth_fset->finish;
    $sth_fset_feature->finish;
    $sth_last->finish;
}


#get/set for runnable and args
sub runnable {
    my ($self, $runnable) = @_;
    
    if ($runnable) {
        #extract parameters into a hash
        my ($parameter_string) = $self->parameters() ;
        my %parameters;
        if ($parameter_string) {
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs) {
                my ($key, $value) = split (/=>/, $pair);
		$key =~ s/\s+//g;
                $value =~ s/\s+//g;
                $parameters{$key} = $value;
            }
        }
        $self->{'_runnable'} = $runnable->new(-analysis => $self->analysis, %parameters);
    }
    return $self->{'_runnable'};
}

1;
