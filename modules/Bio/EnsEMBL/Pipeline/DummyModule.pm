
use strict;


package Bio::EnsEMBL::Pipeline::DummyModule;

sub new {
  my $class = shift;

  return bless {}, $class;
}



sub run {
  print "RUNNING DUMMY\n";
  sleep(10);
  return;
}


sub write_output {
  print "WRITING DUMMY\n";
  sleep(20);
  return;
}

sub fetch_input {
  print "READING DUMMY\n";
  sleep(30);
  return;
}


1;
