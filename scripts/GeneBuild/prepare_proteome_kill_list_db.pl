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

use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::KillList::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => 'PROTEIN');
my %kill_list = %{$kill_list_object->get_kill_list()};

my @file_info = @$GB_PROTEOME_FILES;
my $protfile  = $GB_PFASTA;
my $pmatch    = $GB_PMATCH;

my %ids;

foreach my $file_info(@file_info){
  my $file  = $file_info->{file_path};
  my $regex = $file_info->{header_regex};
  my $refseq = $file_info->{refseq};
  my $in    = Bio::SeqIO->new(
			      -format => 'fasta',
			      -file   => $file,
                          );
  if(-e $protfile){
    print "Protfile ".$protfile." already exists, these ".
      "entries will be appended to the end of the file\n";
    print "Do you want this? answer y/n\n";
    my $reply = <>;
    chomp;
    if($reply =~ /^n/i){
      print "You need to delete or change the name of ".$protfile." before rerunning\n";
      exit(0);
    }
  }
  my $out = Bio::SeqIO->new(
                            -format => 'fasta',
                            -file   => ">>".$protfile,
                           );
 SEQ:while(my $seq = $in->next_seq){
    my $parseable_string = $seq->id." ".$seq->desc;


    my ($id) = $parseable_string =~ /$regex/;
    if(!$id){
#      throw($regex." failed to parse an id out of ".
#            $parseable_string)
      warn($regex." failed to parse an id out of ".
            $parseable_string);
        next SEQ;
    }
    if($seq->seq =~ /XXXXX/){
      print "protein ".$id." has 5 or more consecutive Xs ".
        "in the sequence skipping\n";
      next SEQ;
    }
    if($refseq && !($id =~ /NP/)){
      print "Refseq protein which is not an NM skipping $id\n";
      next SEQ;
    }
    
    my $no_version_id = $id;
    $no_version_id =~ s/\.\d+//;
    if(exists $kill_list{$no_version_id}){
      print $id." on kill list as $no_version_id\n";
      next SEQ;
    }
    if($ids{$id}){
      print $id." has already appeared\n";
      next SEQ;
    }

    $seq->desc('');
    $seq->id($id);
    print "Adding ".$id." to hash\n";
    $ids{$id} = 1;
    my $seq_string = $seq->seq;
    $seq_string =~	s/U/X/g;
    $seq->seq($seq_string);
    $out->write_seq($seq);
  }
}

&test_protfile;

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


