use strict;
use warnings;

#require 'syscall.ph';

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
  Returntype : List of integers first 3 are the load averages.
  Exceptions : 
  Caller     : 

=cut

sub get {

  my $self = shift;

  # use /proc/loadavg to avoid shelling out to uptime every time
  open (LOAD, "/proc/loadavg");
  my @result = split / /, <LOAD>;
  close (LOAD);

  #print "$result[0] $result[1] $result[2]\n";

  return @result;

}


1;
