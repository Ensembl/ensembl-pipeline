# A simple wrapper script for running jobs through GridEngine.
#
#!/usr/local/bin/perl

print "Running on machine " . `hostname` . "\n";
my $ok = 1;
if (exists($ENV{PREEXEC})) {
  print "Pre exec = " . $ENV{PREEXEC} . "\n";
  if (system($ENV{PREEXEC})) {
    print STDERR "Failed pre exec " . $ENV{PREEXEC} . "\n";
    $ok = 0;
  }
}
if ($ok) {
  foreach my $arg (@ARGV) {
    print "Arg = $arg\n";
  }
  if (system(@ARGV)) {
    print STDERR "system @ARGV failed: $?";
  }
} else {
  print "Didn't run job\n";
}

print "Output file: " . $ENV{FINAL_STDOUT} . "\n";
print "Error  file: " . $ENV{FINAL_STDERR} . "\n";
print "Output path is " . $ENV{SGE_STDOUT_PATH} . "\n";
print "Error path is " . $ENV{SGE_STDERR_PATH} . "\n";
my $cmd = "rcp " . $ENV{SGE_STDOUT_PATH} ." " . $ENV{SGE_O_HOST} . ":" . $ENV{FINAL_STDOUT};
system ($cmd);
my $cmd = "rcp " . $ENV{SGE_STDERR_PATH} ." " . $ENV{SGE_O_HOST} . ":" .$ENV{FINAL_STDERR};
system ($cmd);

unlink $ENV{SGE_STDERR_PATH};
unlink $ENV{SGE_STDOUT_PATH};
