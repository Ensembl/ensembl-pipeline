
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
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

Module which inherits from StateInfoContainer and overrides specific
methods to implement the incremental updating of the Blast analysis used
in the Finished Pipeline.

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

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

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
	my $ana_hash;
	my $anaAd = $self->db->get_AnalysisAdaptor();

	my $sth = $self->prepare(
		q{
    SELECT analysis_id , db_version
    FROM input_id_analysis
    WHERE input_id = ? }
	);
	$sth->execute($inputId);
    # Get a list of completed analysis for this input_id
	while ( my $row = $sth->fetchrow_arrayref ) {
		my $analysis = $anaAd->fetch_by_dbID( $row->[0] );
		my $version  = $row->[1];
		next unless ($analysis);
		$ana_hash->{$analysis->logic_name} = [$analysis,$version];
	}
	# Process the list to return completed analysis that are:
	# 1) non-updatable (i.e. Trf, RepeatMasker, ...)
	# 2) up to date (i.e. same current and saved db versions)
	# 3) out of date but depending analysis not run yet
	# For example: Est2genome_other_raw has 
	# saved_version = 04-Oct-09 (101) and current_version = 02-Jun-10 (103)
	# but depending analysis Est2genome_other not run so include it in the list
	# of completed analysis so that the depending analysis get run. 
	foreach(keys %$ana_hash){
		my $analysis        = $ana_hash->{$_}->[0];
		my $saved_version   = $ana_hash->{$_}->[1];
		my $current_version = $analysis->db_version;
		if ( !$current_version ||                  # test points 1
		      $saved_version eq $current_version ) # test points 2
        {
            push( @result, $analysis );
        } elsif(my ($filter_ana) = /(.*)_raw/) { # test point 3
            push( @result, $analysis ) unless (
                $ana_hash->{$filter_ana}       ||
                $ana_hash->{$filter_ana.'_SW'} ||
                $ana_hash->{$filter_ana.'_TR'}
            );
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

	throw("Must provide inputId")
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
