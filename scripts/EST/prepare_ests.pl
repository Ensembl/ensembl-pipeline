#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::SeqIO;
use Getopt::Long;

my @infiles;
my $outfile;
my $kill_list_file;

&GetOptions(
            'kill_list:s'    => \$kill_list_file,
            'infiles:s'      => \@infiles,
            'outfile:s'      => \$outfile,
           );


@infiles = split(/,/,join(',',@infiles));


my %kill_list;

if ($kill_list_file) {
  %kill_list = %{get_kill_list($kill_list_file)};
}

my %ids;

my $out = Bio::SeqIO->new(
                          -format => 'fasta',
                          -file   => ">".$outfile,
                         );


foreach my $file (@infiles){
  my $in    = Bio::SeqIO->new(
			      -format => 'fasta',
			      -file   => $file,
                             );
 print $file . "\n";
 SEQ:while(my $seq = $in->next_seq){
    my $id = $seq->id;
    
    if($kill_list{$id}){
      print $id." on kill list\n";
      next SEQ;
    }

    if($ids{$id}){
      print $id." has already appeared\n";
      next SEQ;
    }

    $ids{$id} = 1;
    $out->write_seq($seq);

  }
  $in->close;
}
$out->close;

sub get_kill_list {
  my ($kill_list) = @_;
  my %kill_list = ();
  die "kill list is not defined" if(!$kill_list);
  open KILL_LIST_FH, $kill_list or die "can't open $kill_list";
  while (<KILL_LIST_FH>) {
    my @element = split;
    if (scalar(@element) == 0) {	# blank or empty line
      next;
    }
    $kill_list{$element[0]} = 1;
  }
  close KILL_LIST_FH or die "file error for $kill_list";
  return \%kill_list;
}

