# Mar 6, 2006 5:20:22 PM
#
# Created by Mustapha Larbaoui <ml6@sanger.ac.uk>

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::Finished::StateinfoContainer

=head1 SYNOPSIS

my $dbobj = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(
	-host   => $dbhost,
	-dbname => $dbname,
	-user   => $dbuser,
	-pass   => $dbpass,
	-port   => $dbport,
);

$sic = $dbobj->get_StateInfoContainer;

=head1 DESCRIPTION

Module which inherits from StateInfoContainer and overrides specific methods to implement the incremental
updating of the Blast analysis used in the Finished Pipeline.

=head1 FEEDBACK

=head1 AUTHOR - Mustapha Larbaoui

Mustapha Larbaoui E<lt>ml6@sanger.ac.ukE<gt>

=head1 CONTACT

Post general queries to B<anacode@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::Finished::StateInfoContainer;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;

@ISA = qw( Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer );

=head2 fetch_db_version

Fetches the db version of a specific input ID and analysis
which have been completed.
Takes one string - input ID
and one analysis object - Bio::EnsEMBL::Analysis
and returns the db version as a string.
Throws an exception if the analysis object does not have a type.

=cut

sub fetch_db_version {
	my ( $self, $inputId, $analysis) = @_;
	my $db_version = '';

	throw("[$analysis] is not a Bio::EnsEMBL::Pipeline::Analysis object")
	  unless $analysis->isa("Bio::EnsEMBL::Pipeline::Analysis");

	my $sth = $self->prepare(
	q {
	    SELECT db_version
    	FROM input_id_analysis
    	WHERE input_id = ?
    	AND analysis_id = ? }
	);
	$sth->execute($inputId,$analysis->dbID);

	while ( my $row = $sth->fetchrow_arrayref ) {
		$db_version = $row->[0];
	}

	return $db_version;
}

=head2 fetch_analysis_by_input_id

Fetches all analyses on a specific input ID
which have been completed and don't need to be updated
Takes one string - input ID
and returns a ref to a list of Bio::EnsEMBL::Analysis.

=cut

sub fetch_analysis_by_input_id {
	my ( $self, $inputId ) = @_;

	my @result;
	my @row;

	my $anaAd = $self->db->get_AnalysisAdaptor();

	my $sth = $self->prepare(
		q{
    SELECT analysis_id , db_version
    FROM input_id_analysis
    WHERE input_id = ? }
	);
	$sth->execute($inputId);

	while ( my $row = $sth->fetchrow_arrayref ) {
		my $analysis = $anaAd->fetch_by_dbID( $row->[0] );
		my $version  = $row->[1];
		next unless ($analysis);
		if (   !$version
			|| $analysis->logic_name =~ /Halfwise/
			|| $analysis->db_version eq $version )
		{
			#print "Analysis " . $analysis->logic_name . " V. $version\n";
			push( @result, $analysis );
		}
		else {
			#print "Analysis ". $analysis->logic_name." V. $version need to be updated to V. ". $analysis->db_version."\n";
		}
	}

	return \@result;
}

=head2 store_input_id_analysis

Stores an input ID, type, analysis and db version searched [optionally runtime info].
Takes an input ID (as string), Bio::EnsEMBL::Analysis object, is_dbversion_saved Boolean
and [optionally runtime_info as string].
Throws an exception if any of the inputs are invalid or if
the analysis object does not have a type.
 called by: Job->run_module

=cut

sub store_input_id_analysis {
	my ( $self, $inputId, $analysis, $host, $is_dbversion_saved, $run_time ) =
	  @_;
	throw("[$analysis] is not a Bio::EnsEMBL::Pipeline::Analysis object")
	  unless $analysis->isa("Bio::EnsEMBL::Pipeline::Analysis");

	throw("Invalid inputId [$inputId]")
	  unless $inputId;

	throw("No type defined in analysis obj")
	  unless defined( $analysis->input_id_type );

	if ( defined($run_time) ) {
		my $db_version = $is_dbversion_saved ? $analysis->db_version : 0;
		my $sth = $self->prepare(
			qq{
      		REPLACE INTO input_id_analysis
      		(input_id, input_id_type, analysis_id, created, runhost, db_version , result)
      		values (?, ?, ?, now(), ?, ? , ?)
      	}
		);
		$sth->execute( $inputId, $analysis->input_id_type, $analysis->dbID,
			$host, $db_version, $run_time );
	}
	else {
		my $db_version = $is_dbversion_saved ? $analysis->db_version : 0;
		my $sth = $self->prepare(
			qq{
      	REPLACE INTO input_id_analysis
      	(input_id, input_id_type, analysis_id, created, runhost, db_version)
      	values (?, ?, ?, now(), ?, ?)
      	}
		);
		$sth->execute( $inputId, $analysis->input_id_type, $analysis->dbID,
			$host, $db_version );
	}
}

1;
