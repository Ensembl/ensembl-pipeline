use strict;

=head1 NAME

  prepare_proteome.pl

=head1 SYNOPSIS

  prepare_proteome.pl

=head1 DESCRIPTION

  prepare_proteome.pl prepares a fasta file of protein sequences
  input files (also in fasta format) specified in the Scripts.pm
  config file. This file is needed for pmatch comparisons and 
  its creation is the first part of the GeneBuild.

  The file has a description line consisting solely of the 
  accession number after the leading >
  All U are replaced by X to prevent pmatch complaining.

  The final part of the script does a tiny pmatch test run to 
  reveal any problems persisting in the file that would prevent 
  later full scale pmatches from running smoothly.

=head1 OPTIONS

  Options are to be set in GeneBuild config files

     GB_KILL_LIST   location of text file listing Swissprot IDs to ignore
     GB_PROTEOME_FILES the locations of the files to parse and the regex to parse them with
     GB_PMATCH location of pmatch executable
     GB_PFASTA where to write the output file too
=cut


use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);


my @file_info = @$GB_PROTEOME_FILES;
my $protfile  = $GB_PFASTA;
my $pmatch    = $GB_PMATCH;

my %kill_list = %{&get_kill_list($GB_KILL_LIST)} if($GB_KILL_LIST);

my %ids;

foreach my $file_info(@file_info){
  my $file  = $file_info->{file_path};
  my $regex = $file_info->{header_regex};
  my $in    = Bio::SeqIO->new(
			      -format => 'fasta',
			      -file   => $file,
                          );
  my $out = Bio::SeqIO->new(
                            -format => 'fasta',
                            -file   => ">>".$protfile,
                           );
 SEQ:while(my $seq = $in->next_seq){
    my $parseable_string = $seq->id." ".$seq->desc;
    my ($id) = $parseable_string =~ /$regex/;
    if(!$id){
      throw($regex." failed to parse an id out of ".
            $parseable_string)
    }
    if($kill_list{$id}){
      print $id." on kill list\n";
      next SEQ;
    }
    if($ids{$id}){
      print $id." has already appeared\n";
      next SEQ;
    }
    $seq->desc('');
    $seq->id($id);
    $ids{$id} = 1;
    my $seq_string = $seq->seq;
    $seq_string =~	s/U/X/g;
    $seq->seq($seq_string);
    $out->write_seq($seq);
  }
}

&test_protfile;

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



sub test_protfile {

  # set up a temporary file
  my $time = time;
  chomp ($time);
  my $tmpfile = "cup.$$.$time.fa";
  open (SEQ, ">$tmpfile") or die "can't open $tmpfile\n";
  print SEQ ">test_seq\n";
  print SEQ 'cctgggctgcctggggaagcacccagggccagggagtgtgaccctgcaggctccacacaggactgccagaggcacac';
  close SEQ;

  # do a pmatch test run
  print "starting pmatch test ... \n";
  open(PM, "$GB_PMATCH -D $protfile $tmpfile | ") or die "Can't run $GB_PMATCH\n";
  while(<PM>) {
    print $_;
  }
  close PM;

  # tidy up
  unlink $tmpfile;

}
