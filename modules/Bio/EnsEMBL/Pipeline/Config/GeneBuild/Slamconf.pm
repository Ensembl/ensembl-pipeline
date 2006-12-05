# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::GeneBuild::Slamconf - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Slamconf;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Slamconf qw(  );

=head1 DESCRIPTION

Slamconf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Slamconf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%Slamconf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::Config::GeneBuild::Slamconf;

use strict;
use vars qw( %Slamconf );



# Hash containg information for Slam-run. The base working-database is passed by SlamDB.pm

%Slamconf = (  

              SLAM_ORG1_NAME => 'H.sapiens'  ,                     # (valid species : R.Norvegicus, H.sapiens or M.musulus)
              SLAM_ORG2_NAME => 'M.musculus' ,


              # Slam-options
              
              SLAM_BIN => '/acari/work6a/jhv/project_slam/slam_prog/slam',
              SLAM_PARS_DIR => '/acari/work6a/jhv/project_slam/slam_prog/Pars',
              SLAM_MAX_MEMORY_SIZE => '2572864',                   # default is 1572864 (1.5 gig)
              SLAM_MINLENGTH => '250',
              SLAM_MAXLENGTH => '300000',                          # Slam.pm- default is 100500 bp. max length of regions to compare

              # DNA-database for first species
              SLAM_ORG1_DNA_DB_USER => 'ensro' ,
              SLAM_ORG1_DNA_DB_PASS => '' ,
              SLAM_ORG1_DNA_DB_NAME => 'homo_sapiens_core_30_35c',
              SLAM_ORG1_DNA_DB_HOST => 'ecs2' ,
              SLAM_ORG1_DNA_DB_PORT => '3364' ,


              # DNA-database for second species 

              SLAM_ORG2_DNA_DB_USER => 'ensro' ,
              SLAM_ORG2_DNA_DB_PASS => '' ,
              SLAM_ORG2_DNA_DB_NAME => 'mus_musculus_core_30_33f' ,
              SLAM_ORG2_DNA_DB_HOST => 'ecs2' ,
              SLAM_ORG2_DNA_DB_PORT => '3364' ,

              # database where the slam-predictions for the second species are written to
              SLAM_ORG2_RESULT_DB_USER  => '' ,
              SLAM_ORG2_RESULT_DB_PASS  => '' ,
              SLAM_ORG2_RESULT_DB_NAME  => '' ,
              SLAM_ORG2_RESULT_DB_HOST  => '',
              SLAM_ORG2_RESULT_DB_PORT  => ''
           );




    ################################################################################


sub import {
my ($callpack) = caller(0); # Name of the calling package
my $pack = shift; # Need to move package off @_

# Get list of variables supplied, or else

my @vars = @_ ? @_ : keys( %Slamconf );
return unless @vars;

# Predeclare global variables in calling package
eval "package $callpack; use vars qw("
. join(' ', map { '$'.$_ } @vars) . ")";
die $@ if $@;


foreach (@vars) {
if ( defined $Slamconf{ $_ } ) {
no strict 'refs';
# Exporter does a similar job to the following
# statement, but for function names, not
# scalar variables:
*{"${callpack}::$_"} = \$Slamconf{ $_ };
} else {
die "Error: File Slamconf.pm : $_ not known\n";
}
}
}

1;
