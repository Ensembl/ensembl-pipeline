use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::Load;

=head2 new

  Arg        : None
  Example    : $load = Bio::EnsEMBL::Pipeline::Load->new();
  Description: Creates an object that can be used to query the load on the current machine.
  Returntype : Bio::EnsEMBL::Pipeline::Load
  Exceptions :
  Caller     :

=cut

sub new {

  my $class = shift;

  my $self = bless {}, $class;

  return $self;

}


=head2 get

  Arg [1]    : None
  Example    : $
  Description: Returns the load of the machine over the last 1, 5 and 15 minutes.
  Returntype : List of integers
  Exceptions : 
  Caller     : 

=cut

sub get {

  my $self = shift;

  my $uptime = `uptime`;

  #print "uptime=" . $uptime . "\n";

  $uptime =~ /^.*load average:\s*(\d+\.\d\d),\s(\d+\.\d\d),\s(\d+\.\d\d)$/;  # why doesn't backreferencing work?
  #print "bits: " . $1 . " " . $2 . " ". $3 . "\n";

  return ($1, $2, $3);

}


1;
