package Bio::EnsEMBL::Pipeline::RunnableDB::Fathom;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::Fathom;
use Bio::EnsEMBL::Pipeline::Runnable::Fathom;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


sub fetch_input{
  my ($self) = @_;

  $self->throw("No input id") unless defined($self->input_id);
  my ($chr, $start, $end) = $self->input_id =~ m!$SLICE_INPUT_ID_REGEX!;
  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my $slice =  $slice_adaptor->fetch_by_chr_start_end
    ($chr, $start, $end) or $self->throw("Unable ot fetch Slice");
  $self->query($slice);
  my @genes = @{$slice->get_all_Genes};
  my $fathom = Bio::EnsEMBL::Pipeline::Runnable::Fathom->new
    (
     -query => $self->query,
     -genes => \@genes,
     -fathom => $self->analysis->program_file,
     -hmmfile => $self->analysis->db_file,
    );

  $self->runnable($fathom);

  return 1;
}



sub write_output{
  my ($self) = @_;

  $self->throw("Can't write results to a dir which doesn't exist ".
               $FATHOM_RESULTS_DIR) unless(-d $FATHOM_RESULTS_DIR);

  my $file = $FATHOM_RESULTS_DIR."/".
    $self->input_id.".$$.".$self->analysis->logic_name.".results";
  print STDERR "have file ".$file."\n";
  open(OUT, ">".$file) or $self->throw("Couldn't open ".$file);
  foreach my $runnable ($self->runnable){
    my $output = $runnable->output;
    print STDERR "have output ".$output."\n";
    my %hash = %$output;
    print STDERR "have hash ".%hash."\n";
    my @names = keys(%hash);
    print STDERR "have names ".@names."\n";
    foreach my $name (@names){
      my $score = $hash{$name};
      print OUT $name." ".$score."\n";
      print STDERR $name." ".$score."\n";
    }
  }
}


sub output{
  my ($self) = @_;

  return $self->runnable->output;
}

1;
