#!/usr/local/bin/perl5.6.1 -w
package Bio::EnsEMBL::Pipeline::AnaSubmission;

use strict;
use DBI;
use Carp;
use Time::Local qw( timelocal );
use Exporter;
use vars qw( @ISA @EXPORT_OK );
use Hum::Submission qw{ get_user };
use Hum::SubmissionConf;

@ISA = ('Exporter');
@EXPORT_OK = qw( 
                set_db_args
                sub_db
                prepare_statement
                get_db
                );

=pod

=head2 ref_from_query( SQL )

Returns a reference to an array of anonymous
arrays containing the results from running the
B<SQL> query on the database.  

=cut
##
##sub ref_from_query {
##    my( $query ) = @_;
##
    ##my $dbh = sub_db();

##    my $sth = $dbh->prepare( $query );
##    $sth->execute;
##    return $sth->fetchall_arrayref;
##}
##
##sub_db returns the database handle
{
    my ($user , $dbname , $host, $port);
    my( $db );

    sub set_db_args{
        ## take in user name password port etc, and use them to connect to the db
        my %args_hash = @_;
        
        $user   = $args_hash{user}  ;        
        $dbname = $args_hash{db_name};
        $host   = $args_hash{host};
        $port   = $args_hash{port};
        
    }

   
    
    sub sub_db{
        unless ($db) {
            
            # go to default values unless given already
            unless ($user) { $user= get_user();};
            unless ($dbname) { $dbname = 'submissions'} ;
            unless (defined $host && defined $port){
                ($host, $port)=Hum::SubmissionConf::localisation();
            }
                  
            
            
            if (my $test = $ENV{'SUBMISSION_TEST_DB'}) {
                # Legal formats:
                #
                # submissions_test
                # submissions_test:3306
                # submissions_test@ecs2b
                # submissions_test@ecs2b:3306
            
                if ($test =~ /^([\w_\$]+)(?:\@([\w\-\.]+))?(?::(\d+))?$/) {
                    $dbname = $1;
                    $host   = $2 if $2;
                    $port   = $3 if $3;
                } else {
                    confess "Can't parse SUBMISSION_TEST_DB environment variable '$test'";
                }
            }
#            warn "DBI:mysql:host=$host;port=$port;database=$dbname:", $user;
            # Make the database connection
            $db = DBI->connect("DBI:mysql:host=$host;port=$port;database=$dbname",
                               $user, undef, {RaiseError => 1, PrintError => 0})
                or die "Can't connect to submissions database as '$user' ",
                DBI::errstr();
        }
        return $db;
    }

    my( @active_statement_handles );
    
    sub prepare_statement {
        my( $text ) = @_;
        
        my $sth = sub_db()->prepare($text);
        push(@active_statement_handles, $sth);
        return $sth;
    }


     sub get_db{
        if ($db){
            return $db;
        }else{
            warn "Not returning a db;";
            return;
        }
        
            
     }

}



1;

__END__

=head1 NAME - Submission

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk
