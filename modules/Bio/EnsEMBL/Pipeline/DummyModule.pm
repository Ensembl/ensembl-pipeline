
use strict;


package Bio::EnsEMBL::Pipeline::DummyModule;

sub new {
  my $class = shift;

  return bless {}, $class;
}



sub run {
  print "RUNNING DUMMY\n";
  #sleep(int(rand(50)));
  return;
}


sub write_output {
  print "WRITING DUMMY\n";
  #sleep(int(rand(10)));
  return;
}

sub fetch_input {
  print "READING DUMMY\n";
  #sleep(int(rand(10)));
  return;
}


1;
