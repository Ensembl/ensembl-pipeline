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

Bio::EnsEMBL::Pipeline::Runnable::GC

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $seqfile = '/nfs/disk65/mq2/temp/bA151E14.seq'; 
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $seqfile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::GC object
  my $gc = Bio::EnsEMBL::Pipeline::Runnable::GC->new (-QUERY => $seq);
  $gc->workdir($workdir);
  $gc->run();
  my @genes = $gc->output();

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

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 


=cut

package Bio::EnsEMBL::Pipeline::Runnable::GC;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::GC->new (-QUERY => $seq);
    Function:   Initialises GC object
    Returns :   a GC object
    Args    :   A Bio::Seq object 
                (size of window to be analysed)

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);    
  
  $self->{'_flist'} = []; #an array of Bio::SimpleFeatures
  $self->{'_sequence'} = undef; #location of bio::seq object
  $self->{'_window'} = 0; #size of window to be analysed
  
  #print STDERR "args: ", @args, "\n";  

  my($sequence, $window) = $self->_rearrange([qw(QUERY
						 WINDOW)],
					     @args);
  #print "sequence : ".$sequence."\n";
  $self->query($sequence) if ($sequence);

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

=head2 query

    Title   :   query
    Usage   :    $GC->query($seq);
    Function:   sets the sequence the fgenesh object will run on
  and checks it is a Bio::Seq
    Returns :   a seq
    Args    :   A Bio::Seq object 
                

=cut


sub query{
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
    #check seq
    my $seq = $self->query() || $self->throw("Seq required for cpg\n");
   
    
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
   

    my $length = $self->query->length;
    my $seq = $self->query->seq;
    my $window = $self->window;
    #my $gccount = $seq =~ tr/GC/GC/;
    #my $Ncount = $seq =~ tr/N/N/;
    #my $total_gcfraction = $gccount/($length-$Ncount);

    

    my $start_point = 0;
    
    my @features;
    while($start_point<$length)
    {
	my $end_point = $start_point + ($self->window) - 1;
	if($end_point>$length)
	{
	    $end_point=$length;
	}
	my $gcfraction;
	my $chunk = substr($seq, $start_point, $window);
	my $GC_count = $chunk =~ tr/gcGC//;
	my $AT_count = $chunk =~ tr/atAT//;
	#print "gcount = ".$gccount."\n N count = ".$Ncount."\n window size = ".$windowsize."\n length = ".length($chunk)."\n";
	my $tot_bases = $GC_count + $AT_count;
	if($tot_bases == 0){
	  $gcfraction = 0;
	}else{
	  #print "division = ".$division."\n";
	  $gcfraction = $GC_count/$tot_bases;
	}
	#print "gcfraction : ".$gcfraction."\n";# gcount = ".$gccount."\n N count = ".$Ncount."\n window size = ".$windowsize."\n length = ".length($chunk)."\n";
	my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => 'gc',
                            -program_version => '1',
                            -gff_source      => 'gc',
                            -gff_feature     => 'gc_content');
        
        
        $start_point++; ## substr takes first character as '0' - necessary to +1 to write correct co-ordinates to DB
        $end_point ++; ## same as above
        my $gc = Bio::EnsEMBL::SimpleFeature->new
	    (   -seqname => $self->query->id,
		-start   => $start_point,
		-end     => $end_point,
		-strand => 0,
		-score   => $gcfraction,
		-analysis => $analysis_obj);  

	$gc->display_label('');

	$start_point = $end_point ; # +1 already added to $end point above - not necesssary here
	#$end_point = $end_point + $window ;

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
    Ahem    :   does it really?!?

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_flist'}};
}















