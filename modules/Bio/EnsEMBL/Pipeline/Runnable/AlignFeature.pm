#!/usr/local/bin/perl

#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::AlignFeature

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::AlignFeature->new(
                                             -genomic  => $genseq,
					     -features => $features,			  
                                             );
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::AlignFeature->new(-genomic => $seq);


    foreach my $f (@features) {
	$obj->addFeature($f);
    }

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

package Bio::EnsEMBL::Pipeline::Runnable::AlignFeature;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
use Bio::EnsEMBL::Pipeline::MiniSeq;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;

#compile time check for executable
use Bio::EnsEMBL::Analysis::Programs qw(est_genome pfetch); 
use Bio::Seq;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $genomic, $features ) = $self->_rearrange(['GENOMIC',
						   'FEATURES'], @args);
       
    $self->throw("No genomic sequence input")           unless defined($genomic);
    $self->throw("[$genomic] is not a Bio::PrimarySeq") unless $genomic->isa("Bio::PrimarySeq");

    $self->genomic_sequence($genomic) if $genomic; 

    $self->{_features} = [];

    if (defined($features)) {
	if (ref($features) eq "ARRAY") {
	    my @f = @$features;
	    
	    foreach my $f (@f) {
		
		if ($f->isa("Bio::EnsEMBL::FeaturePair")) {
		    $self->addFeature($f);
		} else {
		    $self->warn("Can't add feature [$f]. Not a Bio::EnsEMBL::FeaturePair");
		}
	    }
	} else {
	    $self->throw("[$features] is not an array ref.");
	}
    }
    
    return $self; # success - we hope!
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
        $value->isa("Bio::PrimarySeq") || $self->throw("Input isn't a Bio::PrimarySeq");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 addFeature 

    Title   :   addFeature
    Usage   :   $self->addFeature($f)
    Function:   Adds a feature to the object for realigning
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   Bio::EnsEMBL::FeaturePair

=cut

sub addFeature {
    my( $self, $value ) = @_;
    
    if ($value) {
        $value->isa("Bio::EnsEMBL::FeaturePair") || $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair");
	push(@{$self->{_features}},$value);
    }
}


=head2 get_all_FeaturesbyId

    Title   :   get_all_FeaturesById
    Usage   :   $hash = $self->get_all_FeaturesById;
    Function:   Returns a ref to a hash of features.
                The keys to the hash are distinct feature ids
    Returns :   ref to hash of Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut

sub get_all_FeaturesById {
    my( $self) = @_;
    
    my  %idhash;

    FEAT: foreach my $f ($self->get_all_Features) {
	print("Feature is $f " . $f->seqname . "\t" . $f->hseqname . "\n");
	if (!(defined($f->hseqname))) {
	    $self->warn("No hit name for " . $f->seqname . "\n");
	    next FEAT;
	} 
	if (defined($idhash{$f->hseqname})) {
	    push(@{$idhash{$f->hseqname}},$f);
	} else {
	    $idhash{$f->id} = [];
	    push(@{$idhash{$f->hseqname}},$f);
	}
    }

    return (\%idhash);
}


=head2 get_all_Features

    Title   :   get_all_Features
    Usage   :   @f = $self->get_all_Features;
    Function:   Returns the array of features
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut


sub get_all_Features {
    my( $self, $value ) = @_;
    
    return (@{$self->{_features}});
}


=head2 get_all_FeatureIds

  Title   : get_all_FeatureIds
  Usage   : my @ids = get_all_FeatureIds
  Function: Returns an array of all distinct feature hids 
  Returns : @string
  Args    : none

=cut

sub get_all_FeatureIds {
    my ($self) = @_;

    my %idhash;

    foreach my $f ($self->get_all_Features) {
	if (defined($f->hseqname)) {
	    $idhash{$f->hseqname} = 1;
	} else {
	    $self->warn("No sequence name defined for feature. " . $f->seqname . "\n");
	}
    }

    return keys %idhash;
}


=head2 parse_Header

  Title   : parse_Header
  Usage   : my $newid = $self->parse_Header($id);
  Function: Parses different sequence headers
  Returns : string
  Args    : none

=cut

sub parse_Header {
    my ($self,$id) = @_;

    if (!defined($id)) {
	$self->throw("No id input to parse_Header");
    }

    my $newid = $id;

    if ($id =~ /^(.*)\|(.*)\|(.*)/) {
	$newid = $2;
	$newid =~ s/(.*)\..*/$1/;
	
    } elsif ($id =~ /^..\:(.*)/) {
	$newid = $1;
    }
    $newid =~ s/ //g;
    return $newid;
}


sub make_miniseq {
    my ($self,@features) = @_;

    my $strand = $features[0]->strand;
    
    @features = sort {$a->start <=> $b->start} @features;
    
    my $count  = 0;
    my $mingap = $self->minimum_intron;

    my $pairaln  = new Bio::EnsEMBL::Analysis::PairAlign;
    my @genomic_features;

    my $prevend     = 0;
    my $prevcdnaend = 0;

    foreach my $f (@features) {
	print STDERR "Found feature - " . $f->start . "\t" . $f->end . "\t" . $f->strand . "\n";

	if ($f->strand != $strand) {
	    $self->throw("Mixed strands in features set");
	}

	my $start = $f->start;
	my $end   = $f->end;

	if ($start > $end) {
	    my $tmp = $end;
	    $f->end  ($start);
	    $f->start($tmp);;
	}
	
	$start = $f->start - $self->exon_padding;
	$end   = $f->end   + $self->exon_padding;

	my $gap = ($start - $prevend);
	my $cdnagap = abs($f->hstart - $prevcdnaend);
	print STDERR "Feature hstart is " . $f->hstart . "\t" . $prevcdnaend . "\n";
	print STDERR "Padding feature - new start end are $start $end ($cdnagap)\n";


	print STDERR "Count is $count : $mingap " . $gap  . "\n";
	if ($count > 0 && (($gap <  $mingap) || ($cdnagap > 20))) {
	    print(STDERR "Merging exons in " . $f->hseqname . " - resetting end to $end\n");
	    
	    $genomic_features[$#genomic_features]->end($end);
	    $prevend     = $end;
	    $prevcdnaend = $f->hend;
	} else {
	
	    my $newfeature = new Bio::EnsEMBL::SeqFeature;

	    $newfeature->start     ($start);
	    $newfeature->end       ($end);
	    $newfeature->strand    (1);
	    $newfeature->attach_seq($self->genomic_sequence);

	    push(@genomic_features,$newfeature);
	    
	    print(STDERR "Added feature $count: " . $newfeature->start  . "\t"  . 
		                                    $newfeature->end    . "\t " . 
		                                    $newfeature->strand . "\n");

	    $prevend = $end;
	    $prevcdnaend = $f->hend; 
	    print STDERR "New end is " . $f->hend . "\n";

	}
	$count++;
    }

    # Now we make the cDNA features

    my $current_coord = 1;
    
    if ($strand == 1) {
	@genomic_features = sort {$a->start <=> $b->start } @genomic_features;
    } elsif ($strand == -1) {
	print STDERR "Reverse strand - reversing coordinates\n";

	@genomic_features = sort {$b->start <=> $a->start } @genomic_features;
    } else {
	$self->throw("Invalid value for strand [$strand]");
    }

    foreach my $f (@genomic_features) {
	$f->strand(1);
	my $cdna_start = $current_coord;
	my $cdna_end   = $current_coord + ($f->end - $f->start);
	
	my $tmp = new Bio::EnsEMBL::SeqFeature(-start => $cdna_start,
					       -end   => $cdna_end,
					       -strand => $strand);
	
	my $fp  = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
						-feature2 => $tmp);
	
	$pairaln->addFeaturePair($fp);
	
	$self->print_FeaturePair($fp);

	$current_coord = $cdna_end+1;
    }
	
    my $miniseq = new Bio::EnsEMBL::Pipeline::MiniSeq(-id        => 'Genomic',
						      -pairalign => $pairaln);

    my $newgenomic = $miniseq->get_cDNA_sequence;
    print ("New genomic sequence is " . $newgenomic->seq . "\n");
    return $miniseq;

}

sub minimum_intron {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_minimum_intron} = $arg;
    }

    return $self->{_minimum_intron} || 1000;
}

    
sub exon_padding {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_padding} = $arg;
    }

    return $self->{_padding} || 20;
}

sub print_FeaturePair {
    my ($self,$nf) = @_;
    
    print(STDERR "FeaturePair is " . $nf->id    . "\t" . 
	  $nf->start . "\t" . 
	  $nf->end   . "\t(" . 
	  $nf->strand . ")\t" .
	  $nf->hseqname  . "\t" . 
	  $nf->hstart   . "\t" . 
	  $nf->hend     . "\t(" .
	  $nf->hstrand  . ")\n");
}

=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my @ids = get_Sequence(@id)
  Function: Fetches all sequences with ids in the array
  Returns : ref to hash of Bio::PrimarySeq keyed by id
  Args    : none

=cut

sub get_Sequence {
    my ($self,@id) = @_;


ID:    foreach my $id (@id) {

	if (defined($self->{_seq_cache}{$id})) {
	    return $self->{_seq_cache}{$id};
	} 

	my $newid = $self->parse_Header($id);

	next ID unless defined($newid);
	print(STDERR "New id is $newid [$id]\n");

	open(IN,"pfetch -q $newid |") || $self->throw("Error fetching sequence for id [$newid]");

	my $seq;
	
	while (<IN>) {
	    chomp;
	    $seq .= $_;
	}
	
	if (!defined($seq) || $seq eq "no match") {
	    open(IN,"efetch -q $newid |") || $self->throw("Error fetching sequence for id [$newid]");
           
	   $seq = "";
	    
	    while (<IN>) {
		chomp;
		$seq .= $_;
	    }
	}
	
	if (!defined($seq)) {
	    $self->throw("Couldn't find sequence for $newid [$id]");
	}
    
	my $seq = new Bio::PrimarySeq(-id  => $newid,
				      -seq => $seq);
	
	$self->{_seq_cache}{$id} = $seq;

	return $seq;
    }
}

sub get_all_Sequences {
    my ($self,@id) = @_;

    my $seqstr;
    my @newid;

    foreach my $id (@id) {
	my $newid = $self->parse_Header($id);
	push(@newid,$newid);

	print(STDERR "New id is $newid [$id]\n");

	$seqstr .= $newid . " ";
    }

    open(IN,"pfetch -q $seqstr |") || $self->throw("Error fetching sequence for id [$seqstr]");
	
    my $count = 0;
    foreach my $id (@id) {
	my $seq = <IN>;
	chomp($seq);
	if ($seq ne "no match") {
	    $self->{_seq_cache}{$id} = new Bio::PrimarySeq(-seq => $seq,
							   -id  => $newid[$count]);
	}
	$count++;
    }
	
    SEQ: foreach my $id (@id) {
	my $seq   = $self->{_seq_cache}{$id};

	next SEQ unless !defined($seq);
	

	my $newid = $self->parse_Header($id);

	next SEQ unless defined($newid);
	next SEQ if $newid eq "";
	print(STDERR "New id :$newid:$id\n");

	open(IN,"efetch -q $newid |") || $self->throw("Error fetching sequence for id [$newid]");
	    
	$seq = "";
	    
	while (<IN>) {
	    chomp;
	    $seq .= $_;
	}
	
	if (!defined($seq)) {
	    $self->warn("Couldn't find sequence for $newid [$id]");
	}
	
	my $seq = new Bio::PrimarySeq(-id  => $newid,
				      -seq => $seq);

	print("Found seq for $id  $seq\n");

	$self->{_seq_cache}{$id} = $seq;
    }

}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs est2genome on each distinct feature id
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;
    

    my @ids = $self->get_all_FeatureIds;

#    $self->get_all_Sequences(@ids);

    foreach my $id (@ids) {
	my $hseq = $self->get_Sequence(($id));

	if (!defined($hseq)) {
	    $self->throw("Can't fetch sequence for id [$id]\n");
	}

	my $eg = new Bio::EnsEMBL::Pipeline::Runnable::Est2Genome(-genomic => $self->genomic_sequence,
								  -est     => $hseq);

	$eg->run;

	my @f = $eg->output;
	foreach my $f (@f) {
	    print("Aligned output is " . $f->id . "\t" . $f->start . "\t" . $f->end . "\n");
	}
	push(@{$self->{_output}},@f);

    }
}

sub minirun {
    my ($self) = @_;

    my $idhash = $self->get_all_FeaturesById;
    my @ids    = keys %$idhash;

#    $self->get_all_Sequences(@ids);

    ID: foreach my $id (@ids) {

	my $features = $idhash->{$id};
	my @exons;

	print(STDERR "Processing $id\n");
	print(STDERR "Features = " . scalar(@$features) . "\n");

	next ID unless (scalar(@$features) > 1);

	eval {
	    my $miniseq = $self->make_miniseq(@$features);
	    my $hseq    = $self->get_Sequence($id);
	    print("Hseq $id " . $hseq->seq . "\n");
	    if (!defined($hseq)) {
		$self->throw("Can't fetch sequence for id [$id]\n");
	    }
	    
	    my $eg = new Bio::EnsEMBL::Pipeline::Runnable::Est2Genome(-genomic => $miniseq->get_cDNA_sequence,
								      -est     => $hseq);
	    
	    $eg->run;
	    
	    my @f = $eg->output;
	    my @newf;
	    
	    foreach my $f (@f) {
		print(STDERR "Aligned output is " . $f->id    . "\t" . 
		      $f->start . "\t" . 
		      $f->end   . "\t(" . 
		      $f->strand . ")\t" .
		      $f->hseqname  . "\t" . 
		      $f->hstart   . "\t" . 
		      $f->hend     . "\t(" .
		      $f->hstrand  . ")\n");
		
		my @newfeatures = $miniseq->convert_FeaturePair($f);
		
		push(@newf,@newfeatures);
		
		foreach my $nf (@newf) {
		    
		    print(STDERR "Realigned output is " . $nf->id    . "\t" . 
			  $nf->start     . "\t" . 
			  $nf->end       . "\t(" . 
			  $nf->strand    . ")\t" .
			  $nf->hseqname  . "\t" . 
			  $nf->hstart    . "\t" . 
			  $nf->hend      . "\t(" .
			  $nf->hstrand   . ")\n");
		    
		}
	    }
	    
	    push(@{$self->{_output}},@newf);
	};
	if ($@) {
	    print STDERR "Error running est2genome for " . $features->[0]->hseqname . " [$@]\n";
	}
    }
}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_output'}};
}

sub _createfeatures {
    my ($self, $f1score, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id,
        $f1source, $f2source, $f1strand, $f2strand, $f1primary, $f2primary) = @_;
    
    #create analysis object
    my $analysis_obj    = new Bio::EnsEMBL::Analysis
                                (-db              => undef,
                                 -db_version      => undef,
                                 -program         => "est_genome",
                                 -program_version => "unknown",
                                 -gff_source      => $f1source,
                                 -gff_feature     => $f1primary,);
    
    #create features
    my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start =>  $f1start,
                                              -end =>     $f1end,
                                              -seqname =>      $f1id,
                                              -strand =>  $f1strand,
                                              -score =>   $f1score,
                                              -source =>  $f1source,
                                              -primary => $f1primary,
                                              -analysis => $analysis_obj );
 
     my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start =>  $f2start,
                                                -end =>    $f2end,
                                                -seqname =>$f2id,
                                                -strand => $f2strand,
                                                -score =>  undef,
                                                -source => $f2source,
                                                -primary =>$f2primary,
                                                -analysis => $analysis_obj );
    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
                                             -feature2 => $feat2) ;
 
    $self->_growfplist($fp); 
}

sub _growfplist {
    my ($self, $fp) =@_;
    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
}

sub _createfiles {
    my ($self, $genfile, $estfile, $dirname)= @_;
    
    #check for diskspace
    my $spacelimit = 0.1; # 0.1Gb or about 100 MB
    my $dir ="./";
    unless ($self->_diskspace($dir, $spacelimit)) 
    {
        $self->throw("Not enough disk space ($spacelimit Gb required)");
    }
            
    #if names not provided create unique names based on process ID    
    $genfile = $self->_getname("genfile") unless ($genfile);
    $estfile = $self->_getname("estfile") unless ($estfile);    
    #create tmp directory    
    mkdir ($dirname, 0777) or $self->throw ("Cannot make directory '$dirname' ($?)");
    chdir ($dirname) or $self->throw ("Cannot change to directory '$dirname' ($?)"); 
    return ($genfile, $estfile);
}
    

sub _getname {
    my ($self, $typename) = @_;
    return  $typename."_".$$.".fn"; 
}

sub _diskspace {
    my ($self, $dir, $limit) =@_;
    my $block_size; #could be used where block size != 512 ?
    my $Gb = 1024 ** 3;
    
    open DF, "df $dir |" or $self->throw ("Can't open 'du' pipe");
    while (<DF>) 
    {
        if ($block_size) 
        {
            my @L = split;
            my $space_in_Gb = $L[3] * 512 / $Gb;
            return 0 if ($space_in_Gb < $limit);
            return 1;
        } 
        else 
        {
            ($block_size) = /(\d+).+blocks/i
                || $self->throw ("Can't determine block size from:\n$_");
        }
    }
    close DF || $self->throw("Error from 'df' : $!");
}


sub _deletefiles {
    my ($self, $genfile, $estfile, $dirname) = @_;
    unlink ("$genfile") or $self->throw("Cannot remove $genfile ($?)\n");
    unlink ("$estfile") or $self->throw("Cannot remove $estfile ($?)\n");
    chdir ("../");
    rmdir ($dirname) or $self->throw("Cannot remove $dirname \n");
}

1;


