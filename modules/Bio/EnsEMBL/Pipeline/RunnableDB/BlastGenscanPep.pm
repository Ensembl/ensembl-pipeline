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

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides constructed from 
assembling genscan predicted features to peptide sequence. The resulting blast hits are
written back as FeaturePairs.
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
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
#use Data::Dumper;

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
                -input_id:  Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    
    $self->{'_featurepairs'}= [];
    $self->{'_genseq'}      = undef;
    $self->{'_transcripts'} = [];
    $self->{'_runnable'}    = [];
    $self->{'_input_id'}    = undef;
    $self->{'_parameters'}  = undef;
    
    my ( $dbobj, $input_id, $analysis) = 
            $self->_rearrange (['DBOBJ', 'INPUT_ID', 'ANALYSIS'], @args);
    
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
    
    return $self;
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

#converts parameters from string to hash
sub formatted_parameters {
    my ($self) = @_;
     
        #extract parameters into a hash
        my ($parameter_string) = $self->parameters;
        my %parameters;
        if ($parameter_string)
        {
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs)
            {
                my ($key, $value) = split (/=>/, $pair);
                $key =~ s/\s+//g;
                $parameters{$key} = $value;
            }
        }
        $parameters {'-blast'}  = $self->analysis->program;
        $parameters {'-db'}     = $self->analysis->db_file;
        
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
    print STDERR "Running against ".scalar($self->transcripts)." predictions\n";
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
    my @runnable = $self->runnable();
    #align the output to the contig
    foreach my $run (@runnable)
    {
        foreach my $transcript ($self->transcripts)
        {
            if ($run->clone->id eq $transcript->id)
            {
                print STDERR "MATCHED: ".$run->clone->id." with ".$transcript->id."\n";
                $self->align_hit_to_contig($run, $transcript);
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
        my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	    $feat_Obj->write($contig, @features);
    }
    return 1;
}

# This function creates a hash which is used to map between the exon genomic position
# and a position within the genscan predicted peptide. The hash is then matched
# against blast peptide hits to return a set of featurepairs of exons and blast
# peptides
sub align_hit_to_contig {
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
        $ex_align {'pep_end'}   = ($ex_align {'pep_start'} + length($expep))-1;
        #Codons split across exons due to phase are shared by both exons
        $ex_align {'pep_end'}   += 1 if ($ex_align{'end_phase'} != 0);
        $ex_align {'pep_start'} -= 1 if ($ex_align{'phase'} != 0);
        
        #subtract a triplet because the stop codon isn't part of the peptide 
        #$ex_align {'trans_end'} -= 3 if ($exon->translate->seq =~ /\*/);
        
        push (@exon_aligns, \%ex_align);
        
        $dna_align {'exon_dna_limit'} += $exon->length;   
        #print "Exon: ".$ex_align {'name'}
        #        ." PEP ".$ex_align {'pep_start'}." - ".$ex_align {'pep_end'}
        #        ." GEN ".$ex_align {'gen_start'}." - ".$ex_align {'gen_end'}
        #        ." SPh ".$ex_align {'phase'}." EPh ".$ex_align {'end_phase'}."\n";
    }
    
    $dna_align {'pep_limit'} = $dna_align {'exon_dna_limit'}/3;      
    
    #map each feature to 1 or more exons
    foreach my $fp ($run->output)
    {   
        unless (($fp->end - $fp->start)+1 <= $dna_align{'pep_limit'})
        {
            $self->throw("Feature length (".$fp->start."-".$fp->end. 
               ") is larger than peptide (".$dna_align{'pep_limit'}.")\n");
        }
        #find each matching exon
        my (@aligned_exons);
        foreach my $ex_align (@exon_aligns)
        {
            if ($ex_align->{'pep_end'} >= $fp->start)
            {
                push (@aligned_exons, $ex_align);
                last if ($ex_align->{'pep_end'} >= $fp->end);
            }
        }
        #create sets of featurepairs mapping peptide features to exons
        $self->create_peptide_featurepairs($fp, @aligned_exons);
    } 
}

# This function takes a blast peptide feature hit and a set of matching exons and
# creates a set of featurepairs aligned to genomic coordinates. It will split
# features if they cross exon boundaries
sub create_peptide_featurepairs {
    my ($self, $fp, @aligned_exons) = @_;
    #create featurepairs
    
    foreach my $ex_align (@aligned_exons)
    {
        my ($ex_start, $ex_end, $pep_start, $pep_end, $start_phase, $end_phase);
        #This splits features across multiple exons and records phases
        if ($ex_align->{'pep_start'}  < $fp->start)
        {
            #feature starts inside current exon
            $ex_start   = $ex_align->{'gen_start'}
                            + (($fp->start - $ex_align->{'pep_start'})*3);
            $start_phase= 0;
            $pep_start  = $fp->hstart;
            #print "Start inside exon - "
        }
        else
        {
            #feature starts in a previous exon or absolute start of current exon
            $ex_start   = $ex_align->{'gen_start'};
            $start_phase= $ex_align->{'phase'};
            $pep_start  = $fp->hstart + ($ex_align->{'pep_start'} - $fp->start);
            #print "Start outside or equal to exon -"
        }
        
        if ($$ex_align{'pep_end'}    > $fp->end)
        {
            #feature ends in current exon
            $ex_end     = $ex_align->{'gen_start'}
                            + (($fp->end -  $ex_align->{'pep_start'})*3)+2;
            $end_phase  = 0;
            $pep_end    = $fp->hend;
            #print "End inside exon\n";
        }
        else
        {
            #feature ends in a later exon or absolute end of current exon
            $ex_end     = $ex_align->{'gen_end'};
            $end_phase  = $ex_align->{'end_phase'};
            $pep_end    = $fp->hstart + ($ex_align->{'pep_end'} - $fp->start);
            #print "End outside or equal to exon\n";
        }
        
        my $start_frac = $ex_align->{'phase'} + 1;
        my $end_frac   = (( 3 - $$ex_align{'end_phase'})%3) + 1;
        
        my $dna_feat = Bio::EnsEMBL::SeqFeature->new (
                                -seqname    =>  $ex_align->{'name'},
                                -start      =>  $ex_start, 
                                -end        =>  $ex_end,
                                -strand     =>  $ex_align->{'strand'},
                                -score      =>  $fp->score,
                                -p_value    =>  $fp->p_value,
                                -percent_id =>  $fp->percent_id,
                                -analysis   =>  $fp->analysis,
                                -primary_tag=>  $fp->primary_tag,
                                -source_tag =>  $fp->source_tag, 
                                -phase      =>  $start_phase,  
                                -end_phase  =>  $end_phase );
        
        my $pep_feat = Bio::EnsEMBL::Pep_SeqFeature->new (
                                -seqname    =>  $fp->hseqname,
                                -start      =>  $pep_start,
                                -end        =>  $pep_end,
                                -strand     =>  $fp->hstrand,
                                -start_frac =>  $start_frac,
                                -end_frac   =>  $end_frac,
                                -score      =>  $fp->score,
                                -p_value    =>  $fp->p_value,
                                -percent_id =>  $fp->percent_id,
                                -analysis   =>  $fp->analysis,
                                -primary_tag=>  $fp->primary_tag,
                                -source_tag =>  $fp->source_tag );
                                    
        my $featurepair = Bio::EnsEMBL::FeaturePair->new (
                                -feature1   => $dna_feat,
                                -feature2   => $pep_feat );
    
        $self->featurepairs($featurepair);    
    }   
}

sub featurepairs {
    my ($self, $fp) = @_;
    if ($fp)
    {
        $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair") 
                unless $fp->isa("Bio::EnsEMBL::FeaturePairI");
        push (@{$self->{'_featurepairs'}}, $fp);
    }
    #print STDERR   "FEATURES: ".(@{$self->{'_featurepairs'}})."\n";
    return @{$self->{'_featurepairs'}};
}
