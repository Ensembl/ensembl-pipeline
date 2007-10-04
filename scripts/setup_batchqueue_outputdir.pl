#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Getopt::Long;

my $h;
my $help;

GetOptions(
           'help!' => \$help,
           'h!' => \$h,
          ) or perldoc();

if($help || $h){
  perldoc();
}

if(! -d $DEFAULT_OUTPUT_DIR){
  throw("DEFAULT_OUTPUT_DIR is not defined") if(!$DEFAULT_OUTPUT_DIR);
  create_directory($DEFAULT_OUTPUT_DIR);
}

foreach my $hash(@$QUEUE_CONFIG){
  if($hash->{output_dir}){
    if(! -d $hash->{output_dir}){
      create_directory($hash->{output_dir});
    }
  }
}


if(! -d $DEFAULT_OUTPUT_DIR){
  throw($DEFAULT_OUTPUT_DIR." doesn't exist");
}

foreach my $hash(@$QUEUE_CONFIG){
  if($hash->{output_dir}){
    if(! -d $hash->{output_dir}){
     throw($hash->{output_dir}." doesn't exist");
    }
  }
}



sub create_directory{
  my ($dir_path) = @_;
  my $cmd = "mkdir -p $dir_path";
  eval{
    print "Running ".$cmd."\n";
    system($cmd);
  };
  if($@){
    throw("Failed to create $dir_path $@");
  }
}

sub perldoc{
	exec('perldoc', $0);
	exit(0);
}


=pod

=head1 NAME

setup_batchqueue_outputdir.pl

This script will setup the output directories from
Bio::EnsEMBL::Pipeline::Config::BatchQueue

or mail ensembl-dev@ebi.ac.uk

=head1 SYNOPSIS

This script goes through the output directory settings in BatchQueue and
make sures those output directories exist and if they dont creating them
plush any directories that lead up to them

As this script using mkdir -p to construct the directorys though you need
to ensure you have permission to create directories in the locations you
are asking the script to otherwise it wont work

=head1 OPTIONS

-h/-help print out the perdocs

=head1 EXAMPLES

perl setup_batchqueue_outputdir.pl

=head1 SEE ALSO

Bio::EnsEMBL::Pipeline::Config::BatchQueue and the_ensembl_pipeline_infrastructure.txt from ensembl-doc

=cut
