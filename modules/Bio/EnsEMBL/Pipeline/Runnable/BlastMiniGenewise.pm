#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new
    ('-genomic'    => $genseq,
     '-features'   => $features,
     '-seqfetcher' => $seqfetcher);
    
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::PrimarySeqI;
use Bio::SeqIO;
use Bio::DB::RandomAccessI;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_INPUTID_REGEX
					);

require "Bio/EnsEMBL/Pipeline/pipeConf.pl";

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my( $genomic, $ids, $seqfetcher, $endbias) = $self->_rearrange([qw(GENOMIC
																																		 IDS
																																		 SEQFETCHER
																																		 ENDBIAS)],
                                           @args);
  
  $self->throw("No genomic sequence input")            unless defined($genomic);
  $self->throw("No seqfetcher provided")               unless defined($seqfetcher);
	$self->throw("No ids arrary ref provided")           unless defined($ids);

  $self->throw("[$genomic] is not a Bio::PrimarySeqI") unless $genomic->isa("Bio::PrimarySeqI");
	
	$self->ids($ids)                                     if defined($ids);
  $self->genomic_sequence($genomic)                    if defined($genomic);
  $self->endbias($endbias)                             if defined($endbias);
  $self->seqfetcher($seqfetcher)                       if defined($seqfetcher);

  return $self;
}

sub ids {
  my ($self,$ids) = @_;

	if (!defined($self->{_idlist})) {
		$self->{_idlist} = [];
	}
	if (defined($ids)) {
    if (ref($ids) eq "ARRAY") {
      push(@{$self->{'_idlist'}},@$ids);
    } else {
      $self->throw("[$ids] is not an array ref.");
    }
  }
	return @{$self->{_idlist}};
}

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 endbias

    Title   :   endbias
    Usage   :   $self->endbias($endbias)
    Function:   Get/set method for genewise endbias
    Returns :   
    Args    :   

=cut

sub endbias {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_endbias'} = $arg;
    }

    if (!defined($self->{'_endbias'})) {
      $self->{'_endbias'} = 0;
    }    

    return $self->{'_endbias'};
}

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::EnsEMBL::Pipeline::SeqFetcher object
    Args    :   Bio::EnsEMBL::Pipeline::SeqFetcher object

=cut

sub seqfetcher {
  my( $self, $value ) = @_;    
  if ($value) {
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
}


=head2 run

  Title   : run
  Usage   : $self->run()
  Function: 
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;

    my @features = $self->run_blast;
    print STDERR "there are ".@features." features from reblast\n";
    unless (@features) {
        print STDERR "Contig has no associated features\n";
        return;
    }

		my @newf;
		foreach my $f( @features) {
			my $newf = new Bio::EnsEMBL::FeaturePair(-feature1 => $f->feature2,
																							-feature2 => $f->feature1);
			#print $newf->gffstring . "\n";
			push(@newf,$newf);
		}

    my $mmg = new Bio::EnsEMBL::Pipeline::Runnable::MultiMiniGenewise('-genomic'    => $self->genomic_sequence,
																																			'-features'   => \@newf,
																																			'-seqfetcher' => $self->seqfetcher,
																																			'-endbias'    => $self->endbias);

    $mmg->run;
    
    my @f = $mmg->output;
    #print STDERR "Multi Mini Genewise has ".@f." output\n";
    push(@{$self->{'_output'}},@f);

}

sub run_blast {
    my ($self) = @_;

    my @seq         = $self->get_Sequences;
    my @valid_seq   = $self->validate_sequence(@seq);
    #print STDERR "there are ".@valid_seq." valid sequences\n";
    my $blastdb     = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(-sequences => [$self->genomic_sequence],
								    -type      => 'DNA');
    print STDERR "\n";
    $blastdb->run;
    print STDERR "\n";
    my @features;
    my @aligns;
    my $dbname = $blastdb->dbname;
    my @sorted_seqs = sort {$a->id cmp $b->id} @valid_seq;
    foreach my $seq (@sorted_seqs) {
      # First sort out the header parsing. Blergh! cb25.NA_057.31208-61441 Slice, no descrtipion 
     
      #print STDERR "ID ".$self->genomic_sequence->id."\n";
      if($GB_INPUTID_REGEX && $self->genomic_sequence->id =~ /$GB_INPUTID_REGEX/){
	$::fasta_header_re{$dbname} = $GB_INPUTID_REGEX;
      }elsif ($self->genomic_sequence->id =~ /^(.*)\|(.*)\|(.*)/) {
	$::fasta_header_re{$dbname} = '^.*\|(.*)\|.*';
      } elsif ($self->genomic_sequence->id =~ /^..\:(.*)/) {
	$::fasta_header_re{$dbname} = '^..\:(.*)';
      }else {
	$::fasta_header_re{$dbname} = '^(\w+)\s+';
      }
      
      my $run = new Bio::EnsEMBL::Pipeline::Runnable::Blast(-query    => $seq,
							    -program  => 'wutblastn',
							    -database => $blastdb->dbfile,
							    -filter => 0,
							   );
      $run->run;
      print STDERR $run->output." gapped alignments come out of blast\n";
      push(@aligns,$run->output);
    }
    print STDERR "there are ".@aligns." gapped alignments\n";
    foreach my $align(@aligns){
     # my @fps = $align->ungapped_features;
      push(@features, $align);
    }
    #print STDERR @features." ungapped alignments come out of blast\n";
    $blastdb->remove_index_files;
    unlink $blastdb->dbfile;
    #print STDERR "my first feature ".$features[0]."\n";
    return @features;
}
    
sub get_Sequences {
    my ($self) = @_;

    my @seq;

    foreach my $id ($self->ids) {
        my $seq = $self->get_Sequence($id);

        if (defined($seq) && $seq->length > 0) {
            push(@seq,$seq);
        } else {
            print STDERR "Invalid sequence for $id - skipping\n";
        }
    }

    return @seq;

}

sub validate_sequence {
    my ($self,@seq) = @_;
    my @validated;

    foreach my $seq (@seq) {

        my $sequence = $seq->seq;

        if ($sequence !~ /[^acgtn]/i) {
            push (@validated, $seq);
        } else {
            $_ = $sequence;
            my $len = length ($_);
            my $invalidCharCount = tr/bB/xX/;

            if ($invalidCharCount / $len > 0.05) {
                $self->warn("Ignoring ".$seq->display_id()
                    ." contains more than 5% ($invalidCharCount) "
                    ."odd nucleotide codes ($sequence)\n Type returns "
                    .$seq->moltype().")\n");
            } else {
                $seq->seq($_);
                push (@validated, $seq);
            }
        }
    } 
    return @validated;  
}

=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id)
  Function: Fetches sequences with id $id
  Returns : Bio::PrimarySeq
  Args    : none

=cut
    
sub get_Sequence {
    my ($self,$id) = @_;
    my $seqfetcher = $self->seqfetcher;
    my $seq;

    if (!defined($id)) {
      $self->warn("No id input to get_Sequence");
    }  
    
    eval {
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };

    if($@) {
      $self->warn("Problem fetching sequence for id [$id] $@\n");
      return undef;
    }
    
    if(!defined($seq)){
      $self->warn("Could not find sequence for [$id]");
    }

    return $seq;
	}

1;
