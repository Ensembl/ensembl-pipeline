#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( throw );
use Getopt::Long;

use IO::File;

sub usage {
  print <<USAGE_END
list_tests.pl is a script for displaying what tests can be run for
a given species.  This is what the command line should look like:

  ./list_tests.pl -species homo_sapiens -conf_file TestDB.conf

USAGE_END
}

my $species   = 'homo_sapiens';
my $conf_file = 'TestDB.conf';
my $help      = 0;

if ( !GetOptions( 'species|s=s'   => \$species,
                  'conf_file|c=s' => \$conf_file,
                  'help|h!'       => \$help, ) )
{
  usage();
  exit(1);
}

if ($help) { usage(); exit(0) }

if ( !-f $conf_file ) {
  throw("Can not find configuration file '$conf_file': $!");
}

my $conf = do $conf_file;

my $data_dir = $conf->{'data_dir'};

if ( !-d $data_dir ) {
  throw("Can not find data directory '$data_dir': $!");
}

my $curr_dir = $ENV{'PWD'};
my $zip_file = sprintf( "%s.zip", $species );
my $zip_path = sprintf( "%s/%s", $data_dir, $zip_file );
my $dest_dir = sprintf( "%s/%s", $curr_dir, $species );

if ( !-f $zip_path ) {
  throw("Unable to unpack '$zip_path': $!");
}

if ( !-d $dest_dir ) {
  mkdir($dest_dir);
  my @cmd = ( 'unzip', '-q', $zip_path, '-d', $dest_dir );
  print("Unpacking '$zip_file'...\n");
  system(@cmd) == 0 or
    throw("Error unpacking '$zip_file' using '@cmd': $!");
}
else {
  print("Will use the existing unpacked directory '$dest_dir'\n");
}

my $analysis_file = sprintf( "%s/analysis", $dest_dir );

my $fh = IO::File->new($analysis_file);
if ( !defined($fh) ) {
  throw("Can not open analysis file '$analysis_file': $!");
}

printf( "%-15s %-15s\n", 'logic_name', 'module' );
printf( "%-15s %-15s\n", '----------', '------' );
while ( my $line = $fh->getline() ) {
  chomp($line);

  my @values = split( /\t/, $line );

  my $logic_name = $values[2];
  if ( $logic_name =~ /^submit/ ) { next }

  my $module = $values[10];

  printf( "%-15s %-15s\n", $logic_name, $module );
}

$fh->close();
