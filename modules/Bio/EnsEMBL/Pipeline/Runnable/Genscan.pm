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

Bio::EnsEMBL::Pipeline::Runnable::Genscan

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $clonefile = '/nfs/disk65/mq2/temp/bA151E14.seq'; 
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::Genscan object
  my $genscan = Bio::EnsEMBL::Pipeline::Runnable::Genscan->new (-CLONE => $seq);
  $genscan->workdir($workdir);
  $genscan->run();
  my $seqfeature = $genscan->output();
  my @exons = $genscan->output_exons();
  my @genes = $genscan->output_genes();

=head1 DESCRIPTION

This package is based on Genscan2ace.
Genscan takes a Bio::Seq object and runs Genscan on it. The
resulting output is parsed to produce a set of Bio::SeqFeatures. 

=head2 Methods:

=over 4

=item new($seq_obj)

=item    Genscan($path_to_Genscan)

=item    workdir($directory_name)

=item    run()

=item    output()

=item    output_exons()

=item    output_genes()

=back

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Genscan;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object );

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Genscan->new (-CLONE => $seq);
    Function:   Initialises Genscan object
    Returns :   a Genscan Object
    Args    :   A Bio::Seq object 
                (Genscan location and parameter file location optional)

=cut

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
    $self->{_seqfeature} = undef;   #Bio::SeqFeature container for Subseqfeatures
    $self->{_genes} = [];           #an array of arrays of SeqFeatures
    $self->{_clone} = undef;        #location of Bio::Seq object
    $self->{_genscan} = undef;      #location of Genscan script
    $self->{_workdir} = undef;      #location of temp directory
    $self->{_filename} =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;     #file to store results of genscan
    $self->{_protected} =[];        #a list of file suffixes protected from deletion
    $self->{_parameters} =undef;    #location of parameters for genscan
    my($clonefile, $genscan, $parameters) = $self->_rearrange(['CLONE', 'GENSCAN', 'PARAM'], @args);
    
    $self->clone($clonefile) if ($clonefile);
    if ($genscan)       
    {$self->genscan($genscan) ;}
    else                
    { 
        eval 
        {  $self->genscan($self->locate_runnable('genscan')); }; 
        if ($@) 
        {  $self->genscan('/nfs/disk100/humpub/OSFbin/genscan'); }
    }
    if ($parameters)    
    { $self->clone($parameters) ; }
    else                
    {$self->parameters('/nfs/disk100/humpub/OSFbin/HumanIso.smat'); }     
    return $self; # success - we hope!
}


###################
#get/set methods
###################

sub clone {
    my ($self, $seq) = @_;
    if ($seq)
    {
        $seq->isa("Bio::Seq") || $self->throw("Input isn't a Bio::Seq");
        $self->{_clone} = $seq ;
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".genscan");
        $self->seqfeature($seq);
    }
    return $self->{_clone};
}

sub seqfeature {
    my ($self, $seq) = @_;
    if ($seq)
    {
        #new Bio::SeqFeature to hold exon information
        #create analysis object
        my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => 'genscan',
                            -program_version => "unknown",
                            -gff_source      => 'genscan',
                            -gff_feature     => 'exon');
    
        #create and fill Bio::EnsEMBL::Seqfeature objects   
        my $sequence = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $self->clone->id(),
                            -start   => 1,
                            -end     => $self->clone->length(),
                            -strand  => 1,
                            -score   => undef,
                            -frame   => undef,
                            -source  => undef,
                            -primary => $self->clone->moltype(),
                            -analysis => $analysis_obj);
        $self->{_seqfeature} = $sequence;
    }
    return $self->{_seqfeature};
}

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

=head2 genscan

    Title   :   genscan
    Usage   :   $obj->genscan('/nfs/disk100/humpub/OSFbin/genscan');
    Function:   Get/set method for the location of genscan
    Args    :   File path (optional)

=cut

sub genscan {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Genscan not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_genscan} = $location ;
    }
    return $self->{_genscan};
}

=head2 parameters

    Title   :   parameters
    Usage   :   $obj->parameters('/nfs/disk100/humpub/OSFbin/HumanIso.smat');
    Function:   Get/set method for the location of genscan parameters
    Args    :   File path (optional)

=cut

sub parameters {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Genscan parameters not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_parameters} = $location ;
    }
    return $self->{_parameters};
}

=head2 workdir

    Title   :   workdir
    Usage   :   $obj->wordir('~humpub/temp');
    Function:   Get/set method for the location of a directory to contain temp files
    Args    :   File path (optional)

=cut

sub add_exon {
    my ($self, $exon) =@_;
    if ($exon)
    {
        $exon->isa("Bio::EnsEMBL::SeqFeature") 
            || $self->throw("Input isn't a Bio::EnsEMBL::SeqFeature");
        #BUG: 'Use of uninitialized value at SeqFeature.pm line 644'
        #FIX: Pass empty string to add_sub_SeqFeature
        $self->seqfeature->add_sub_SeqFeature($exon, '');
    }
}

sub add_gene {
    my ($self, $gene) =@_;
    $gene->isa("Bio::EnsEMBL::SeqFeature") 
            || $self->throw("Input isn't a Bio::EnsEMBL::SeqFeature");
    push(@{$self->{'_genes'}}, $gene);
}

###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs genscan and creates array of sub-seqfeatures
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir) = @_;
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Genscan\n");
    #set directory if provided
    $self->workdir('genscan') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run genscan       
    $self->run_genscan();
    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
}

sub run_genscan {
    my ($self) = @_;
    print "Running genscan on ".$self->filename."\n";
    open (OUTPUT, $self->genscan.' '.$self->parameters.' '.$self->filename.'|')
            or $self->throw("Couldn't open pipe to genscan: $!\n");  
    open (RESULTS, ">".$self->results)
            or $self->throw("Couldn't create file for genscan results: $!\n");
    print RESULTS <OUTPUT>;
    close OUTPUT;
    close RESULTS;
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses Genscan output to give a set of seqfeatures
    Returns :   none
    Args    :   optional filename

=cut


sub parse_results {
    my ($self) = @_;
    my %exon_type = ('Sngl', 'Single Exon',
                     'Init', 'Initial Exon',
                     'Intr', 'Internal Exon',
                     'Term', 'Terminal Exon');
    open (GENSCAN, "<".$self->results)
        or $self->throw ("Couldn't open file ".$self->results.": $!\n");
    if (<GENSCAN> =~ m|NO EXONS/GENES PREDICTED IN SEQUENCE| )
    {
        print "No genes predicted\n";
        return;
    }
    while (<GENSCAN>)
    {
        # Last line before predictions contains nothing
        # but spaces and dashes
        if (/^\s*-[-\s]+$/) 
        { 
            while (<GENSCAN>) 
            {
                my @element = split;
                my %feature;
                
                if ($element[1] && $element[1] =~ /init|term|sngl|intr/i)
                {
                   $feature {name} = $element[0];
                   #arrange numbers so that start is always < end
                   if ($element[2] eq '+')
                   {
                        $feature {start} = $element[3];
                        $feature {end} = $element[4];
                        $feature {strand} = 1;
                   }
                   elsif ($element[2] eq '-')
                   {
                        $feature {start} = $element[4];
                        $feature {end} = $element[3];
                        $feature {strand} = -1;
                   }
                   
                   $feature {score} = $element[11];
                   $feature {frame} = $element[6];
                   $feature {type} = $exon_type{$element[1]};
                   $feature {primary} = 'exon';
                   $feature {source} = 'genscan 1.0';
                   $self->create_exon(\%feature);
                }
                #data ends with line 'predicted peptide sequence(s)'
                elsif ($element[0] && $element[0] =~ /predicted/i)
                {
                    close GENSCAN;
                    return;
                }
            }
        }
    }
}

sub create_exon {
    my ($self, $feat) = @_;
    #$self->create_gene();
    
    #create analysis object
    my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => $feat->{source},
                            -program_version => "unknown",
                            -gff_source      => $feat->{source},
                            -gff_feature     => $feat->{primary});
    
    #create and fill Bio::EnsEMBL::Seqfeature objects   
    my $exon = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $feat->{name},
                            -start   => $feat->{start},
                            -end     => $feat->{end},
                            -strand  => $feat->{strand},
                            -score   => $feat->{score},
                            -frame   => $feat->{frame},
                            -source_tag  => $feat->{source},
                            -primary_tag => $feat->{type},
                            -analysis => $analysis_obj);  
    $self->add_exon($exon);
}

#only called by output_genes()
sub create_genes {
    my ($self) = @_;
    my @exons = $self->output_exons();
    my @ordered_exons = sort { $a->seqname <=> $b->seqname } @exons;
    for (my $ex_index = 0; $ex_index < scalar(@ordered_exons); $ex_index++)
    {
        my $ex = $ordered_exons[$ex_index];
        if ($ex->primary_tag =~ /Single|Initial/) #lone exon or intial exon
        {
            $ex->seqname =~ /\./;     #extract number before '.'
            my $gene_number = $`;     #this is the gene number
            my @sub_exons;    #an array of exons to be loaded
            #create gene as new SeqFeature
            my $gene = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $gene_number,
                            -strand  => $ex->strand,
                            -score   => $ex->score,
                            -frame   => $ex->frame,
                            -source  => $ex->source_tag,
                            -primary => $ex->primary_tag,
                            -analysis => $ex->analysis);
            #start/end point depends on orientation                 
            my $gene_start = $ex->start;
            my $gene_end = $ex->end;
            $gene->start($gene_start) if ($ex->strand eq '1');
            $gene->end($gene_end) if ($ex->strand eq '-1');
            
            for (my $gene_index = $ex_index; 
                    $gene_index < scalar(@ordered_exons);
                    $gene_index++)
            {
                my $current_exon = $ordered_exons[$gene_index];
                $current_exon->seqname =~ /\./; #extract gene num
                if ($gene_number == $`)
                {
                    push (@sub_exons, $current_exon);
                    $gene_start = $current_exon->start;
                    $gene_end = $current_exon->end;
                }
            }
            $gene->start($gene_start) if ($ex->strand eq '-1');
            $gene->end ($gene_end) if ($ex->strand eq '1');
            #now that the bounds are set can load sub_seqfeatures
            foreach my $sub_exon (@sub_exons)
            {
                $gene->add_sub_SeqFeature($sub_exon, ''); #add exon as sub-seqfeature
            }
            $self->add_gene($gene); #add gene to main object
        }
    }
}

##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns a single SeqFeature with exons as sub-SeqFeatures
    Returns :   A single SeqFeature with exons as sub-SeqFeatures
    Args    :   none

=cut

sub output {
my ($self) = @_;
    return $self->{'_seqfeature'};
}

=head2 output_exons

    Title   :   output_exons
    Usage   :   obj->output_exons()
    Function:   Returns an array of SeqFeatures corresponding to exons
    Returns :   An array of SeqFeatures corresponding to exons
    Args    :   none

=cut

sub output_exons {
    my ($self) = @_;
    my @exons = $self->seqfeature->sub_SeqFeature();
    return @exons;
}

=head2 output_genes

    Title   :   output_genes
    Usage   :   obj->output_genes()
    Function:   Returns an array of SeqFeatures representing predicted genes
    Returns :   An array of SeqFeatures (genes) containing sub-seqfeatures (exons)
    Args    :   none

=cut

sub output_genes {
    my ($self) = @_;
    #Decided to create genes here to save uneccesary storage - genes are just reorganised exons
    $self->create_genes();
    return @{$self->{'_genes'}};
}
