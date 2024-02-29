
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::Pipeline::Finished::PipeQueue;

=head1 NAME

Bio::EnsEMBL::Pipeline::Finished::PipeQueue - access to pipe_queue.queue

=head1 DESCRIPTION

This module wraps up access to the otterlive.pipe_queue database.

It is still not in the Ensembl canonical "DBAdaptor" style, but it is
better than having multiple copies of the C<< DBI->connect >> code.

=head1 METHODS

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Pipeline::Config::BatchQueue qw( QUEUE_HOST QUEUE_NAME );

use DBI;


=head2 configure($queue_host, $queue_name)

Replace the default configuration, which is taken from
L<Bio::EnsEMBL::Pipeline::Config::BatchQueue>

=cut

{
    # State used by the class.  Overridden by dequeuer.pl
    my ($queue_host, $queue_name) = ($QUEUE_HOST, $QUEUE_NAME);

    sub configure {
        my $called = shift;
        ($queue_host, $queue_name) = @_;
    }

    sub _get_config {
        return ($queue_host, $queue_name);
    }
}


=head2 qdbh()

Class method.  Returns a new L<DBI> pointing at the currently
configured C<pipe_queue>.

Password comes from L<Net::Netrc> per team practice.

=cut

sub qdbh {
    my ($called) = @_;
    my ( $dbhost, $dbname ) = $called->_get_config;

    my ($dbuser, $dbpass, $dbport) = &get_db_param( $dbhost );
    my $dsn = "DBI:mysql:host=$dbhost;dbname=$dbname;port=$dbport";

    return DBI->connect($dsn, $dbuser, $dbpass,
                        { PrintError => 1, AutoCommit => 1 }); # current defaults
}


=head2 ensure_sth( \$sth, $sql )

Class method.  If the DBI statement handle $sth is not connected,
prepare it from L</qdbh>.

Deprecated: this is to assist cleanup existing code, not intended for
writing new code.

=cut

sub ensure_sth {
    my ($called, $sth_ref, $sql) = @_;

    if (!$$sth_ref || not $$sth_ref->{Database}->ping) {
        $$sth_ref = $called->qdbh->prepare($sql);
    }

    return 1;
}


# XXX:DRY there are already many. fix later.
use Net::Netrc;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
sub get_db_param {
	my ( $dbhost ) = @_;
	my ( $dbuser, $dbpass, $dbport );

	my $ref = Net::Netrc->lookup($dbhost);
	throw("$dbhost entry is missing from ~/.netrc") unless ($ref);
	$dbuser = $ref->login;
	$dbpass = $ref->password;
	$dbport = $ref->account;
	throw(
		"Missing parameter in the ~/.netrc file:\n
			machine " .  ( $dbhost || 'missing' ) . "\n
			login " .    ( $dbuser || 'missing' ) . "\n
			password " . ( $dbpass || 'missing' ) . "\n
			account "
		  . ( $dbport || 'missing' )
		  . " (should be used to set the port number)"
	  )
	  unless ( $dbuser && $dbpass && $dbport );

	return ( $dbuser, $dbpass, $dbport );
}


=head2 push_job($qdbh, $job, $priority)

Add a job to the queue.  If $priority is true it overrides
C<<$job->priority>>.

(This class method has an awkward interface, it is an intermediate
step in extracting pipequeue logic.)  The caller must supply a
returnvalue from L</qdbh>, presumably cached.

=cut

sub push_job {
	my ($called, $dbq, $job, $priority) = @_;

	my $insert = $dbq->prepare(
		qq {
			INSERT INTO queue (created, priority, job_id, host, pipeline, analysis, is_update)
			VALUES ( NOW() , ? , ? , ? , ? , ? , ?)
		}
	);
	my $job_id = $job->dbID;
	my $dbc    = $job->adaptor->db->dbc();
	my $dbname = $dbc->dbname;
	my $host   = $dbc->host;
	my $analysis = $job->analysis->logic_name;
	my $job_priority = $job->priority;
	$job_priority = $priority if $priority;

	my $update = $job->update;

	return $insert->execute( $job_priority, $job_id, $host, $dbname, $analysis, $update );
}


1;

=head1 AUTHOR

Matthew Astley mca@sanger.ac.uk

=cut
