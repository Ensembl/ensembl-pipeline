#
#
# Cared for by Laura Clarke <lec@sanger.ac.uk>
#
# Copyright Laura Clarke
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::C

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq'; 
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::GC object
  my $gc = Bio::EnsEMBL::Pipeline::Runnable::GC->new (-CLONE => $seq);
  $gc->workdir($workdir);
  $gc->run();
  my @genes = $gc->output();
  my @exons = $gc->output_exons();
  my $seqfeature = $gc->output_singlefeature();

=head1 DESCRIPTION



=head2 Methods:

=over 4

=item new($seq_obj)

=item    

=item    workdir($directory_name)

=item    run()

=item    output()

=item   

=item    

=back

=head1 CONTACT

ensembl-dev.ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 


=cut

package Bio::EnsEMBL::Pipeline::Runnable::GC;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::GC->new (-CLONE => $seq);
    Function:   Initialises GC object
    Returns :   a GC object
    Args    :   A Bio::Seq object 
                (sizw of window to be analysed)

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    
  
  $self->{'_flist'} = []; #an array of Bio::SeqFeatures
  $self->{'_sequence'} = undef; #location of bio::seq object
  $self->{'_window'} = 0; #size of window to be analysed
  
  print STDERR "args: ", @args, "\n";  

  my($sequence, $window) = $self->_rearrange([qw(CLONE
						 WINDOW)],
					     @args);
  #print "sequence : ".$sequence."\n";
  $self->clone($sequence) if ($sequence);

  if (defined $window && $window>=0){
      $self->window($window);
  } else {
      $self->window(250); 
  }
  
  return $self; 
}

#################
#GET/SET METHODS#
#################

=head2 clone

    Title   :   clone
    Usage   :    $GC->clone($seq);
    Function:   sets the sequence the fgenesh object will run on
  and checks it is a Bio::Seq
    Returns :   a seq
    Args    :   A Bio::Seq object 
                

=cut


sub clone{
    my ($self, $seq) = @_;
    #print "sequence: ".$seq."\n";
    if(defined $seq){
     unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
      {
	  $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
      }
      $self->{'_sequence'} = $seq ;
 } 
      

  return $self->{'_sequence'};    

}

=head2 window

  Title    : Window
  Usage    : $GC->window(250);
  Function : sets the size of dna sequence to be analysised at once
  Returns  : The window size
  Args     : The window size

=cut


sub window {
    my ($self, $window) = @_;
    if (defined $window) {
      $self->{'_window'} = $window ;
    }
    return $self->{'_window'};
}



##################                
#analysis methods#
##################

=head2 run

    Title   :  run
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs gc script and creates array of features
    Returns :   none
    Args    : none

=cut

sub run {
    my ($self) = @_;
    #set arguments for cpg
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for cpg\n");
   
    
    #running the gc analysis
    $self->{'_flist'} = $self->run_gc();
    
}

=head2 run

    Title   :  run_gc
    Usage   :   $obj->run_gc()
    Function:   calculates the gc contents of the specified windows of dna sequence
    Returns :  array of features
    Args    : none

=cut

sub run_gc {

    my ($self) = @_;
    
    #print "running GC content analysis\n";
   

    my $length = $self->clone->length;
    my $seq = $self->clone->seq;
    my $window = $self->window;
    #my $gccount = $seq =~ tr/GC/GC/;
    #my $Ncount = $seq =~ tr/N/N/;
    #my $total_gcfraction = $gccount/($length-$Ncount);

    

    my $start_point = 0;
    
    my @features;
    while($start_point<$length)
    {
	my $end_point=$start_point+($self->window)-1;
	if($end_point>$length)
	{
	    $end_point=$length;
	}
	
	my $chunk = substr($seq, $start_point, $window);
	my $gccount = $chunk =~ tr/gcGC/gcGC/;
	my $Ncount = $chunk =~ tr/nN/nN/;
	my $windowsize = $end_point-$start_point;
	my $gcfraction = $gccount/($windowsize-$Ncount)*100;
	#print "gcfraction : ".$gcfraction."\n gcount = ".$gccount."\n N count = ".$Ncount."\n window size = ".$windowsize."\n length = ".length($chunk)."\n";
	my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => 'gc',
                            -program_version => '1',
                            -gff_source      => 'gc',
                            -gff_feature     => 'gc_content');
        my $gc = Bio::EnsEMBL::SeqFeature->new
	    (   -seqname => $self->clone->id,
		-start   => $start_point,
		-end     => $end_point,
		-strand => 0,
		-score   => $gcfraction,
		-source_tag  => 'gc',
		-primary_tag => 'gc_content',
		-analysis => $analysis_obj);  

	$start_point = $start_point + $window + 1;
	#$end_point = $end_point + $window + 1;

	push (@features, $gc);
    }

    return \@features;

}



################
#output methods#
################

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of SeqFeatures representing predicted genes 
                with exons stored as SubSeqFeatures.
    Returns :   An array of SeqFeatures the gc content of each widnow
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_flist'}};
}













