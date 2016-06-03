#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 NAME

  check_node.pl

=head1 SYNOPSIS
 
  check_node.pl
  Checks to make sure the estfile and input_id file (if it does not look like chr_name.start-end))
  can be found.
  Intended to be run as the -E command in a bsub

=head1 DESCRIPTION


=head1 OPTIONS

=cut

use warnings ;
use strict;
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Exonerate qw (
							      EST_CHUNKDIR
							      EST_GENOMIC
							      
							     );

my $chunkname = $ARGV[0];

if(!defined $chunkname){
  print STDERR "Usage: check_node.pl chunkname\n";
  exit(1);
}

my $estfiledir = $EST_CHUNKDIR;
my $estfile = $estfiledir . "/" . $chunkname;

my $input_id   = $EST_GENOMIC;

# check to see if can find $estfile
if (! -e $estfile){
    print STDERR "can't find $estfile\n";
    exit(1);
}

# check to see if $input_id is likely to be a file, and if it is, can we find it?
if(!($input_id =~/^\S+\.\d+-\d+$/)){
  if (! -e $input_id){
    print STDERR "can't find $input_id\n";
    exit(1);
  }
}
exit (0);


