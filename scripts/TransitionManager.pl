# Copyright EMBL-EBI 2000
# Author: Arne Stabenau
# Creation: 11.07.2000

# walk through object_state and
# state_description table
# Fire off jobs.
# if states are involved maybe good to have statenames in key.

use Getopt::Long;
use strict;
use Bio::EnsEMBL::Pipeline::ControlDB;


use vars qw/ $cdb $portnumber $sleepseconds $configurationFile /;

$ENV{PATH} = "/nfs/disk21/stabenau/src/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Transitions:".$ENV{PATH};


&GetOptions( 
	     'port=i'    => \$portnumber,
	     'sleep=i' => \$sleepseconds,
	     'c=s' => \$configurationFile,
	    );

# if portnumber, open port for external control
$cdb = Bio::EnsEMBL::Pipeline::ControlDB->new();

$SIG{CHLD} = sub { wait };

while( 1 ) {
  # check the port for any commands
  # postponed
  
  # read in all inTransition lines from states which need monitoring
  my @toMonitor = $cdb->get_inTransit_toMonitor;
  # for each state_nickname call the TransitionMonitor
  my %toMonitorHash;
  # print STDERR ("passed get_inTransit_toMonitor\n" );
  # print STDERR (@toMonitor);
  
  for (@toMonitor) {
    push( @{$toMonitorHash{$_->{object_class}."_".$_->{state_nickname}}}, $_ );
  }

  for my $module ( keys %toMonitorHash ) {
    my @idList = map { $_->{object_state_id} } @{$toMonitorHash{$module}};
      
    eval {
      print STDERR ("Monitor: $module\n" );
      require "Bio/EnsEMBL/Pipeline/Monitor/${module}.pm";
      "Bio::EnsEMBL::Pipeline::Monitor::${module}"->check( @idList );
    }; 
    if( $@ ) {
      print STDERR ( "Problem with $module!\n" );
      print STDERR ( $@, "\n" );
    }
  }
  

  # read in all objects not in Transition and not in endState
  my @fireTransit = $cdb->get_nonFinal_nonTransit;
  # print STDERR ( "passsed get_nonFinal_nonTransit.\n" );
  print STDERR @fireTransit,"\n";
  my %fireTransitHash;
  for (@fireTransit) {
    push( @{$fireTransitHash{$_->{object_class}."_".$_->{state_nickname}}}, $_ );
  }
  
  for my $stateKey ( keys %fireTransitHash ) {
    my @idList;
    my @objIdList;
    
    @objIdList = map { $_->{object_state_id}} @{$fireTransitHash{$stateKey}};

    if( $fireTransitHash{$stateKey}->[0]->{needsMonitoring} eq 'true' ) {
      @idList = map { $_->{object_id}} @{$fireTransitHash{$stateKey}};
    } else {
      @idList = @objIdList;
    }
    
    # start TransitionModules from the result
    my $pid = fork;
    if ( $pid == 0 ) {
      $cdb->object_toTransit( \@objIdList );
      print STDERR ( "Transit: ". $fireTransitHash{$stateKey}->[0]->{transition_module}, "\n" );
      
      exec( $fireTransitHash{$stateKey}->[0]->{transition_module}, @idList );
      # uups, we are here, something failed on exec.
      $cdb->object_reset_transit( \@objIdList );
      die( "Couldnt exec script ".$fireTransitHash{$stateKey}->[0]->{transition_module} );
    } else {
      if( ! defined $pid ) {
        print STDERR ( "Failed to fork!\n" );
      }
    }
  }
  
  # no wait
  if( defined $sleepseconds ) {
    print STDERR ( "Sleeping $sleepseconds\n" );
    my $sleeptime = $sleepseconds;
    while( $sleeptime > 1 ) {
      $sleeptime -= sleep $sleeptime;
    }
  } 
}
  
