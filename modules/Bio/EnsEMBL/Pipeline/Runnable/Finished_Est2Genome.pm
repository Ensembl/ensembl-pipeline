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

Bio::EnsEMBL::Pipeline::Runnable::Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome->new(
                                             -genomic => $genseq,
                                             -est     => $estseq 
                                             );
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome->new()

=head1 DESCRIPTION

Object to store the details of an est2genome run.
Stores the est2genome matches as an array of Bio::EnsEMBL::FeaturePair

=head2 Methods:

 new,
 genomic_sequence,
 est_sequence,
 run,
 output.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome;

use vars qw(@ISA $verbose);
use strict;
# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Est2Genome);

$verbose = 0;

sub run{
    my ($self, @args) = @_;
    
    # some constant strings
    $self->{'_source_tag'}  = "est2genome";
    my $dirname     = "/tmp/";

    #check inputs
    my $genomicseq = $self->genomic_sequence || $self->throw("Genomic sequence not provided");
    my $estseq     = $self->est_sequence     || $self->throw("EST sequence not provided");

    #extract filenames from args and check/create files and directory
    my ($genname, $estname) = $self->_rearrange(['genomic', 'est'], @args);
    my ($genfile, $estfile) = $self->_createfiles($genname, $estname, $dirname);
    #use appropriate Bio::Seq method to write fasta format files
    {
        my $genOutput = Bio::SeqIO->new(-file => ">$genfile" , '-format' => 'Fasta')
                    or $self->throw("Can't create new Bio::SeqIO from $genfile '$' : $!");
        my $estOutput = Bio::SeqIO->new(-file => ">$estfile" , '-format' => 'Fasta')
                    or $self->throw("Can't create new Bio::SeqIO from $estfile '$' : $!");

        #fill inputs
        $genOutput->write_seq($self->{'_genomic_sequence'}); 
        $estOutput->write_seq($self->{'_est_sequence'});
    }

    my $est_genome_command = $BIN_DIR . "/est2genome -reverse -genome $genfile -est $estfile -space 500000 -out stdout 2>&1 |";
#    my $est_genome_command = $self->est_genome . " -reverse -genome $genfile -est $estfile -space 500000 -out stdout 2>&1 |";
#/nfs/disk100/pubseq/emboss/bin/
        # use -align to get alignment
    print STDERR "\nRunning command $est_genome_command\n\n";

    eval{
	open (ESTGENOME, $est_genome_command) || $self->throw("Can't open pipe from '$est_genome_command' : $!");
	my $firstline = <ESTGENOME>;
	print STDERR "firstline: \t$firstline" if $verbose;
	if ( $firstline =~ m/Align EST and genomic DNA sequences/ ){ 
	    # Catch 'Align EST and genomic DNA sequences'. This comes from STDERR!! [ 2>&1 ]
	    $firstline = <ESTGENOME>;
	    print STDERR "\$firstline (secondline!): \t$firstline" if $verbose;	    
	}

	if( $firstline =~ m/reversed\sest/ ){
	    $self->est_strand(-1);
	} else {
	    $self->est_strand(1);
	}

	if( $firstline =~ m/forward\sgenome/ ){
	    $self->gen_strand(1);
	} else {
	    $self->gen_strand(-1);
	}

	while (<ESTGENOME>) {
	    print STDERR $_ if $verbose;
	    if ($_ =~ /^Segment|^Exon/) {  # We only care about Segments in Exons

		# "gen" = genomic sequence
		my ($primary, $score, $percent_id,
		    $gen_start, $gen_end, $gen_id,
		    $est_start, $est_end, $est_id) = (split)[0,1,2,
							     3,4,5,
							     6,7,8];
		$self->$primary(
		    -score      => $score,
		    -percent_id => $percent_id,
		    -gen_start  => $gen_start,
		    -gen_end    => $gen_end,
		    -gen_id     => $gen_id,
		    -est_start  => $est_start,
		    -est_end    => $est_end,
		    -est_id     => $est_id,
		    -primary    => $primary
		    );
	    }elsif ($_ =~ /Segmentation fault/) {
		$self->warn("Segmentation fault from est_genome\n");
		close (ESTGENOME) or $self->warn("problem closing est_genome: $!\n");
		return(0);
	    }elsif ($_ =~ /(ERROR:.+)/) {
		$self->warn("Error from est_genome: \n<$1>\n");
		close (ESTGENOME) or $self->warn("problem closing est_genome: $!\n");
		return(0);
	    }
	}
	foreach my $seg_array( keys(%{$self->{'_exons'}}) ){
	    my $dnafp = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => $self->{'_exons'}->{$seg_array});
	    $self->add_output($dnafp);
	}
	if(!close(ESTGENOME)){
	    $self->warn("Problems running est_genome when closing pipe: $!\n");
	    return (0);
	}
    };
    $self->_deletefiles($genfile, $estfile);

    if ($@) {
        $self->throw("Error running est_genome:\n$@");
	die $@;
    } else {
        return 1;
    }
}
sub est_strand{
    my ($self,$sign) = @_;
    if(defined($sign) && ($sign eq '1' || $sign eq '-1')){
	$self->{'_est_strand'} = $sign;
    }
    return $self->{'_est_strand'};
}
sub gen_strand{
    my ($self,$sign) = @_;
    if(defined($sign) && ($sign eq '1' || $sign eq '-1')){
	$self->{'_gen_strand'} = $sign;
    }
    return $self->{'_gen_strand'};
}

sub Segment{ # named to match output from est2genome
    my($self) = shift;
    my $p = {
	-score      => undef,
	-percent_id => undef,
	-gen_start  => undef,
	-gen_end    => undef,
	-gen_id     => undef,
	-est_start  => undef,
	-est_end    => undef,
	-est_id     => undef,
	-primary    => undef,
	@_ 
    };
    if( $p->{-gen_end} < $p->{-gen_start} ){
	( $p->{-gen_end}, $p->{-gen_start} ) = ( $p->{-gen_start}, $p->{-gen_end} );
    }
    if( $p->{-est_end} < $p->{-est_start} ){
	( $p->{-est_end}, $p->{-est_start} ) = ( $p->{-est_start}, $p->{-est_end} );
    }
    my $fp = $self->_create_FeaturePair(
					$p->{-score}, $p->{-percent_id},
					$p->{-gen_start}, $p->{-gen_end}, $p->{-gen_id}, 
					$p->{-est_start}, $p->{-est_end}, $p->{-est_id},
					$self->{'_source_tag'}, 
					$self->gen_strand, $self->est_strand,   # est_strand IS stored in the db now!
					$p->{-primary}
					);

    $self->_add_segment_to_exon($fp,$p);
}

sub Exon{ # named to match output from est2genome
    my($self) = shift;
    my $p = {
	-score      => undef,
	-percent_id => undef,
	-gen_start  => undef,
	-gen_end    => undef,
	-gen_id     => undef,
	-est_start  => undef,
	-est_end    => undef,
	-est_id     => undef,
	@_
    };
    #push( @{$self->{'_exon_lines'}} , [ $p->{-gen_start}, $p->{-gen_end}, $p->{-score}, $p->{-percent_id} ] );
    $self->_get_add_exon_lines( [ $p->{-gen_start}, $p->{-gen_end}, $p->{-score}, $p->{-percent_id} ] );
}
sub _add_segment_to_exon{
    my ($self,$fp,$p) = @_;
    my $exon_line_count = scalar(@{$self->{'_exon_lines'}});
    for(my $i = 0; $i < $exon_line_count; $i++){
	# get current __exon data
        my ( $start, $end, $score, $pid ) = ( @{$self->{'_exon_lines'}->[$i]} );
	# check if feature pair belongs to current __exon
        next if ($p->{-gen_start} < $start || $p->{-gen_start} > $end);
	# get last feature pair of current __exon for sanity check
	if( my $prev_fp = $self->_get_last_fp_of_exon($i) ){
	    # sanity check to cater for Bio::EnsEMBL::BaseAlignFeature->_parsefeatures

	    if( $fp->start eq ($prev_fp->end + $fp->strand) 
		|| $fp->end eq ($prev_fp->start + $fp->strand)
		|| $fp->hstart eq ($prev_fp->hend + $fp->hstrand) 
		|| $fp->hend eq ($prev_fp->hstart + $fp->hstrand) ){
		$fp->score($score); 
		$fp->percent_id($pid);
		$self->_add_fp_to_exon($fp,$i);
	    }else{
		print "Splitting the exon\n";
		$self->_change_exon_line($i,[ $start, $prev_fp->end, $score, $pid ]);
		my $new_exon_lines = $self->_get_add_exon_lines([ $fp->start - 1, $end, $score, $pid ]);
 		$fp->score($score); 
 		$fp->percent_id($pid);
		# this adds fp to new exon
		$self->_add_fp_to_exon($fp,(scalar(@{$new_exon_lines}) - 1));
	    }
	}else{
	    $fp->score($score); 
	    $fp->percent_id($pid);
	    $self->_add_fp_to_exon($fp,$i);
	}
    }
}
sub _get_last_fp_of_exon{
    my $self = shift;
    my $exon = shift;
    my $last = 0;
    if(defined($exon)){
	$self->{'_exons'}->{$exon} ||= [];
	$last = scalar( @{$self->{'_exons'}->{$exon}} );
	$last = $last - 1 if $last;
    }   
    return $self->{'_exons'}->{$exon}->[$last] || undef;
}
sub _add_fp_to_exon{
    my $self = shift;
    my $fp = shift;
    my $exon = shift;
    if(defined($exon) && defined($fp)){
	$self->{'_exons'}->{$exon} ||= [];
	push (@{$self->{'_exons'}->{$exon}}, $fp);
    }
}
sub _get_add_exon_lines{
    my $self = shift;
    my $exon_line = shift;
    if(@{$exon_line}){
	push( @{$self->{'_exon_lines'}}, $exon_line ) if scalar(@{$exon_line}) == 4;
    }
    return $self->{'_exon_lines'};
}
sub _change_exon_line{
    my $self = shift;
    my $index = shift;
    my $exon_line = shift;
    if(defined($index) && @{$exon_line}){
	$self->{'_exon_lines'}->[$index] = $exon_line;
    }
    return $self->{'_exon_lines'};
}
sub add_output{
    my ($self, $fp) = @_;
    push (@{$self->{'_output'}}, $fp) if defined ($fp);
}
sub output{
    my ($self) = @_;
    return $self->{'_output'} || [];
}


sub _create_FeaturePair {
    my ( $self, 
	 $f1score, 
	 $f1percent_id, 
	 $f1start, 
	 $f1end, 
	 $f1id, 
	 $f2start, 
	 $f2end, 
	 $f2id, 
	 $f1source, 
	 $f1strand, 
	 $f2strand, 
	 $f1primary ) = @_;
    
    #print "creating feature pair ".$f1primary." ".$f1source." \n";
    my $analysis_obj    = new Bio::EnsEMBL::Analysis
                                (-db              => "none",
                                 -db_version      => "none",
                                 -program         => "est_genome",
                                 -program_version => "none",
                                 -gff_source      => $f1source,
                                 -gff_feature     => $f1primary,);
    #create features
    my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start      =>   $f1start,
                                              -end         =>   $f1end,
                                              -seqname     =>   $f1id,
                                              -strand      =>   $f1strand,
                                              -score       =>   $f1score,
                                              -percent_id  =>   $f1percent_id, 
                                              -source_tag  =>   $f1source,
                                              -primary_tag =>   $f1primary,
                                              -analysis    =>   $analysis_obj );
     
   my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $f2start,
                                              -end         =>   $f2end,
                                              -seqname     =>   $f2id,
                                              -strand      =>   $f2strand,
                                              -score       =>   $f1score,
                                              -percent_id  =>   $f1percent_id, 
                                              -source_tag  =>   $f1source,
                                              -primary_tag =>   $f1primary,
                                              -analysis    =>   $analysis_obj );

    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
                                             -feature2 => $feat2) ;
    return $fp;
}

1;
