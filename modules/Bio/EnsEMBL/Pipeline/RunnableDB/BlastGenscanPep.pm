#!/usr/local/bin/perl -w

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

Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBLoader->new($locator);
my $genscan     = Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$genscan->fetch_input();
$genscan->run();
$genscan->output();
$genscan->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep to add
functionality for reading and writing to databases.
The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::RunnableDB::Blast;
use Bio::EnsEMBL::Analysis::GenscanPeptide;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;

use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::Blast);
=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    
    $self->{'_featurepairs'}= [];
    $self->{'_genseq'}      = undef;
    $self->{'_transcripts'} = [];
    $self->{'_genes'}       = [];
    $self->{'_runnable'}    = [];
    $self->{'_input_id'}    = undef;
    $self->{'_parameters'}  = undef;
    
    my ( $dbobj, $input_id, $analysis, $threshold) = 
            $self->_rearrange (['DBOBJ', 'INPUT_ID', 'ANALYSIS', 'THRESHOLD'], @args);
    
    $self->throw('Need database handle') unless ($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI")  
                unless ($dbobj->isa ('Bio::EnsEMBL::DB::ObjI'));
    $self->dbobj($dbobj);
    
    $self->throw("No input id provided") unless ($input_id);
    $self->input_id($input_id);
    
    $self->throw("Analysis object required") unless ($analysis);
    $self->throw("Analysis object is not Bio::EnsEMBL::Pipeline::Analysis")
                unless ($analysis->isa("Bio::EnsEMBL::Pipeline::Analysis"));
    $self->analysis($analysis);
    
    if ($threshold)
    {
        $self->threshold($threshold);
    }
    else
    {
        $self->threshold(1e-6);
    }   
    return $self;
}

=head2 threshold

    Title   :   threshold
    Usage   :   $obj->threshold($value);
    Function:   Get/set method for threshold p value required for writing
                Feature/FeaturePair to database.
    Args    :   Optional value (p value from Blast)

=cut

sub threshold {
    my ($self, $value) = @_;
    $self->{'_threshold'} = $value if (defined $value);
    return $self->{'_threshold'};
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    print STDERR "Fetching contig $contigid\n";
    my $contig    = $self->dbobj->get_Contig($contigid) 
        or $self->throw("Unable to find contig ($contigid)\n");
    my $genseq    = $contig->primary_seq() 
        or $self->throw("Unable to fetch contig sequence");
    $self->genseq($genseq);
    #need to get features predicted by genscan
    $self->predicted_genes($contig->get_all_PredictionFeatures);
    $self->transcripts($contig->get_genscan_peptides);
}

sub transcripts {
    my ($self, @transcripts) = @_;
    
    if (@transcripts)
    {
        foreach (@transcripts)
        {
            $self->throw("Input $_ is not a Bio::EnsEMBL::Transcript\n")
                unless $_->isa("Bio::EnsEMBL::Transcript");
        }
        push (@{$self->{'_transcripts'}}, @transcripts);
    }
    return @{$self->{'_transcripts'}};
}

sub predicted_genes {
    my ($self, @genes) = @_;
    if (@genes)
    {
        foreach my $gene (@genes)
        {
            $self->throw("Input not a Bio::SeqFeatureI\n")
                unless $gene->isa("Bio::EnsEMBL::SeqFeatureI");
        }
        push (@{$self->{'_genes'}}, @genes);
    }
    return @{$self->{'_genes'}};
}

#converts parameters from string to hash
sub formatted_parameters {
    my ($self) = @_;
     
        #extract parameters into a hash
        my ($parameter_string) = $self->parameters;
        $parameter_string =~ s/\s+//g;
        my @pairs = split (/,/, $parameter_string);
        my %parameters;
        foreach my $pair (@pairs)
        {
            my ($key, $value) = split (/=>/, $pair);
            $parameters{$key} = $value;
        }
        $parameters {'-blast'}  = $self->analysis->program;
        $parameters {'-db'}     = $self->analysis->db;
        
    return %parameters;
}

sub runnable {
    my ($self, @runnable) = @_;
    if (@runnable)
    {
        foreach my $runnable (@runnable)
        {
            $runnable->isa("Bio::EnsEMBL::Pipeline::RunnableI") or
                $self->throw("Input to runnable is not Bio::EnsEMBL::Pipeline::RunnableI");
        }
        push (@{$self->{'_runnable'}}, @runnable);
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;
    #need to pass one peptide at a time
    $self->throw("Input must be fetched before run") unless ($self->genseq);
    foreach my $transcript ($self->transcripts)
    {
        print STDERR "Creating BioPrimarySeq ".$transcript->id."\n";
        my $peptide = Bio::PrimarySeq->new(
                        -id         => $transcript->id,
                        -seq        => $transcript->translate->seq(),
                        -moltype    => 'protein' );
                                        
        my %parameters = $self->formatted_parameters();
        my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(%parameters);
        $runnable->clone($peptide);
        $runnable->threshold($self->threshold());
        $runnable->run();
        $self->runnable($runnable);                                        
    }

    my (@output);
    my @runnable = $self->runnable() 
            or $self->throw("Can't return output - no runnable object(s)");
    #align the output to the contig
    foreach my $run (@runnable)
    {
        foreach my $transcript ($self->transcripts)
        {
            if ($run->clone->id eq $transcript->id)
            {
                print STDERR "MATCHED: ".$run->clone->id." with ".$transcript->id."\n";
                $self->align_to_contig($run, $transcript);
                last;
            }
        }
    }
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    return $self->featurepairs();  
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self) = @_;

    my $db=$self->dbobj();
    my @features = $self->output();
    
    my $contig;
    eval 
    {
        $contig = $db->get_Contig($self->input_id);
    };
    
    if ($@) 
    {
	    print STDERR "Contig not found, skipping writing output to db: $@\n";
    }
    elsif (@features) 
    {
        #should add conditional for evalue here
	    print STDERR "Writing features to database\n";
        foreach my $feature (@features)
        {
            print STDERR ($feature->hseqname()."\t");
        }
        my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	    $feat_Obj->write($contig, @features);
    }
    return 1;
}

# This function creates a hash which is used to map between the exon genomic position
# and a position within the genscan predicted peptide. The hash is then matched
# against blast peptide hits to return a set of featurepairs of exons and blast
# peptides
sub align_to_contig {
    my ($self, $run, $trans) = @_;
    my (%dna_align, @exon_aligns, @featurepairs); #structure for tracking alignment variables
    
    $dna_align {'exons'} = [];
    #exons must be in order of peptide start position
    $trans->sort;
    #calculate boundaries and map exons to translated peptide
    #Note: Exons have an extra 3 bases for the stop codon. Peptides lack this
    foreach my $exon ($trans->each_Exon)
    {
        my %ex_align;
        my $pep = $trans->translate->seq;
        my ($expep) = $exon->translate->seq =~ /[^\*]+/g;
        $self->throw("Exon translation not found in peptide") 
                    unless ($pep =~ /$expep/);
                        
        $ex_align {'name'}      = $exon->id;
        $ex_align {'gen_start'} = $exon->start;
        $ex_align {'gen_end'}   = $exon->end;
        $ex_align {'strand'}    = $exon->strand;
        $ex_align {'phase'}     = $exon->phase;
        $ex_align {'end_phase'} = $exon->end_phase;
        $ex_align {'pep_start'} = index($pep, $expep)+1;
        $ex_align {'pep_end'}   = $ex_align {'pep_start'} + length($expep)-1;
        $ex_align {'trans_start'} = $ex_align{'gen_start'} + ((3- $ex_align{'phase'})%3);
        #This "moves" the dangling phase bases across to the 3 'end of an exon for alignment
        $ex_align {'trans_end'}   = $ex_align{'gen_end'} 
                                    + ((3- $ex_align{'end_phase'})%3);
        #if there is an overhang at the end, a coding codon will be missed from the
        #sequence of exons. Adding 1 covers these blank regions
        $ex_align {'pep_end'}   += 1 if ($ex_align{'end_phase'} != 0);
        #subtract a triplet because the stop codon isn't part of the peptide 
        $ex_align {'trans_end'} -= 3 if ($exon->translate->seq =~ /\*/);
        
        push (@exon_aligns, \%ex_align);
        
        $dna_align {'exon_dna_limit'} += $exon->length;   
    }
    
    $dna_align {'pep_limit'} = $dna_align {'exon_dna_limit'}/3;      
    
    #map each feature to 1 or more exons
    foreach my $fp ($run->output)
    {   
        unless ($fp->start <= $dna_align{'pep_limit'} 
                    && $fp->end <= $dna_align{'pep_limit'})
        {
            $self->throw("Feature coordinates (".$fp->start."-".$fp->end. 
               ") do not fit translated peptide (".$dna_align{'pep_limit'}.")\n");
        }
        #find each matching exon
        my (@aligned_exons);
        foreach my $ex_align (@exon_aligns)
        {
            if ($$ex_align{'pep_end'} >= $fp->start)
            {
                push (@aligned_exons, $ex_align);
                last if ($$ex_align{'pep_end'} >= $fp->end);
            }
        }
        #create sets of featurepairs mapping peptide features to exons
        push (@featurepairs, $self->create_aligned_featurepairs($fp, @aligned_exons));
    }
    return @featurepairs; 
}

# This function takes a blast peptide feature hit and a set of matching exons and
# creates a set of featurepairs aligned to genomic coordinates. It will split
# features if they cross exon boundaries
sub create_aligned_featurepairs {
    my ($self, $fp, @aligned_exons) = @_;
    
    #set the genomic start and end coordinates
    my $first_exon  = $aligned_exons[0];
    my $last_exon   = $aligned_exons[scalar(@aligned_exons)-1];
    my $alignment_start = (($fp->start - $first_exon->{'pep_start'})* 3)  
                            + $first_exon->{'trans_start'};
    my $alignment_end   = (($fp->end - $$last_exon{'pep_start'}+1)* 3)    
                            + $$last_exon{'trans_start'};

    #create featurepairs
    my $prev_pep_end = 0;
    foreach my $ex_align (@aligned_exons)
    {
        my ($ex_start, $ex_end, $pep_start, $pep_end, $start_phase, $end_phase);
        
        #This splits features across multiple exons
        $ex_start   = ($$ex_align{'pep_start'}  <= $fp->start) 
                        ? $alignment_start : $$ex_align{'trans_start'}; 
        $ex_end     = ($$ex_align{'pep_end'}    >= $fp->end) 
                        ? $alignment_end : $$ex_align{'trans_end'} + 1;
        $pep_start  = ($fp->start >= $$ex_align{'pep_start'})
                        ? $fp->hstart : $prev_pep_end + 1; 
        $pep_end    = ($fp->end <= $$ex_align{'pep_end'})
                        ? $fp->hend : $pep_start + (($ex_end-$ex_start)/3)-1 ;    
        
        $prev_pep_end = $pep_end;
        
        $start_phase = $$ex_align{'phase'} + 1;
        $end_phase   = (( 3 - $$ex_align{'end_phase'})%3) + 1;
        
        my $dna_feat = Bio::EnsEMBL::SeqFeature->new (
                                -seqname    => $$ex_align{'name'},
                                -start      => $ex_start, 
                                -end        => $ex_end,
                                -strand     => $$ex_align{'strand'},
                                -score      => $fp->score,
                                -p_value    => $fp->p_value,
                                -percent_id => $fp->percent_id,
                                -analysis   => $fp->analysis,
                                -primary_tag=> $fp->primary_tag,
                                -source_tag => $fp->source_tag );
        
        my $pep_feat = Bio::EnsEMBL::Pep_SeqFeature->new (
                                -seqname    => $fp->hseqname,
                                -start      => $pep_start,
                                -end        => $pep_end,
                                -strand     => 1,
                                -start_frac => $start_phase,
                                -end_frac   => $end_phase,
                                -score      => $fp->score,
                                -p_value    => $fp->p_value,
                                -percent_id => $fp->percent_id,
                                -analysis   => $fp->analysis,
                                -primary_tag=> $fp->primary_tag,
                                -source_tag => $fp->source_tag );
                                    
        my $featurepair = Bio::EnsEMBL::FeaturePair->new (
                                -feature1   => $dna_feat,
                                -feature2   => $pep_feat );
    
        $self->featurepairs($featurepair);    
    }   
}

sub featurepairs {
    my ($self, @fp) = @_;
    
    if (@fp)
    {
        push (@{$self->{'_featurepairs'}}, @fp);
    }
    return @{$self->{'_featurepairs'}};
}
