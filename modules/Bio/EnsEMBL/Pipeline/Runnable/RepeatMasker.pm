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

Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker

=head1 SYNOPSIS


    
=head1 DESCRIPTION

=head2 Methods:


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{_fplist} = [];        #an array of feature pairs
    $self->{_clone} = undef;        #location of Bio::Seq object
    $self->{_repeatmasker} = undef; #location of RepeatMaskerHum script
    $self->{_workdir} = undef;      #location of temp directory
    $self->{_filename} =undef;      #file to store Bio::Seq object
    
    $self->repeatmasker('/tmp_mnt/nfs/disk100/humpub/scripts/RepeatMaskerHum');
    my( $clonefile ) = $self->_rearrange(['CLONE'], @args);
    $self->clone($clonefile) if ($clonefile);       
    return $self; # success - we hope!
}

#################
# get/set methods 
#################

sub clone {
    my ($self, $seq) = @_;
    if ($seq)
    {
        $seq->isa("Bio::Seq") || $self->throw("Input isn't a Bio::Seq");
        $self->{_clone} = $seq ;
        $self->filename($self->clone->id."$$.seq");
    }
    return $self->{_clone};
}

sub filename {
    my ($self, $filename) = @_;
    $self->{_filename} = $filename if ($filename);
    return $self->{_filename};
}

sub repeatmasker {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("RepeatMaskerHum not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_repeatmasker} = $location ;
    }
    return $self->{_repeatmasker};
}

sub workdir {
    my ($self, $directory) = @_;
    if ($directory)
    {
        #mkdir ($directory) unless (-d $directory);
        print "$directory doesn't exist\n" unless (-d $directory);
        $self->{_workdir} = $directory;
    }
    return $self->{_workdir};
}

###########
# Analysis methods
##########
    
sub run {
    my ($self, $dir, $args) = @_;
    #set arguments for repeatmasker
    $args = '-low' unless ($args);
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for RepeatMasker\n");
    #set directory if provided
    $self->workdir('RepeatMaskerOutput') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile();        
    $self->run_repeatmasker($args);
    #parse output of repeat masker 
    $self->parse_repmask();
    #$self->deletefiles();
}

sub run_repeatmasker {
    my ($self, $args) = @_;
    #run RepeatMaskerHum
    print "Running RepeatMasker\n";
    system ($self->repeatmasker." $args ".$self->filename); 
    #or $self->throw("Error running RepeatMasker: $!\n")
    #open repeat predictions
    open (REPOUT, "<".$self->filename.".RepMask.out")
            or $self->throw($self->filename.".RepMask.out not created by RepeatMaskerHum\n");   
    close REPOUT;
}

sub parse_repmask {
    my ($self) = @_;
    print "Parsing output\n";
    open (REPOUT, "<".$self->filename.".RepMask.out")
        or $self->throw("Error opening ".$self->filename."RepMask.out\n");
    #check if no repeats found
    if (<REPOUT> =~ /no repetitive sequences detected/)
    {
        print "RepeatMaskerHum didn't find any repetitive sequences\n";
        close REPOUT;
        return;
    }
    #extract values
    my @output = <REPOUT>;
    for (my $index = 2; $index < scalar(@output); $index++) #loop from 3rd line
    {  
        my @element = split (/\s+/, $output[$index]);  
        my (%feat1, %feat2);
        $feat1 {name} = $element[5];
        $feat1 {score} = $element[1];
        $feat1 {start} = $element[6];
        $feat1 {end} = $element[7];
        #set strand ('+' = 1 and 'c' = -1)
        $feat1 {strand} = 1 if ($element[9] eq '+');
        $feat1 {strand} = -1 if ($element[9] eq 'C');
        $feat2 {name} = $element[10];
        $feat2 {score} = $element[1];
        $feat2 {start} = ($element[12] =~ s/(|)//); #strip away parentheses from $element[11]
        $feat2 {end} = $element[13];
        $feat2 {strand} = 1;
        $self->createfeaturepair(\%feat1, \%feat2); #may need to use references
    }
    close REPOUT;   
}

sub createfeaturepair {
    my ($self, $feat1, $feat2) = @_;
    #some contant strings
    my $source = 'RepeatMaskerHum';
    my $primary = 'similarity';
    
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => $source,
                            -program_version => "unknown",
                            -gff_source      => $source,
                            -gff_feature     => $primary,);
    
    #create and fill Bio::EnsEMBL::Seqfeature objects
    my $seqfeature1 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat1->{name},
                            -start   => $feat1->{start},
                            -end     => $feat1->{end},
                            -strand  => $feat1->{strand},
                            -score   => $feat1->{score},
                            -source  => $source,
                            -primary => $primary,
                            -analysis => $analysis_obj);
    
    my $seqfeature2 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat2->{name},
                            -start   => $feat2->{start},
                            -end     => $feat2->{end},
                            -strand  => $feat2->{strand},
                            -score   => $feat2->{score},
                            -source  => $source,
                            -primary => $primary,
                            -analysis => $analysis_obj);
    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $seqfeature1,
                                             -feature2 => $seqfeature2) ;
    $self->_growfplist($fp);                             
}

sub _growfplist {
    my ($self, $fp) =@_;    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
}

##############
# input/output methods
#############

sub output {
    my ($self) = @_;
    return @{$self->{'_fplist'}};
}


sub writefile {
    my ($self) = @_;
    print "Writing sequence to ".$self->filename."\n";
    #create Bio::SeqIO object and save to file
    my $clone_out = Bio::SeqIO->new(-file => ">".$self->filename , '-format' => 'Fasta')
            or $self->throw("Can't create new Bio::SeqIO from ".$self->filename.":$!\n");
    $clone_out->write_seq($self->clone) 
            or $self->throw("Couldn't write to file ".$self->filename.":$!");
}

sub deletefiles {
    my ($self) = @_;
    #delete all analysis files? probably need to 'glob' or something
    unlink ($self->filename);
}

sub checkdir {
    my ($self) = @_;
    #check for disk space
    my $spacelimit = 0.01;
    $self->throw("Not enough disk space ($spacelimit required):$!\n") 
                        unless ($self->diskspace('./', $spacelimit));
    my $dir = $self->workdir();
    chdir ($dir) or $self->throw("Cannot change to directory $dir ($!)\n");
}

sub diskspace {
    my ($self, $dir, $limit) =@_;
    my $block_size; #could be used where block size != 512 ?
    my $Gb = 1024 ** 3;
    
    open DF, "df $dir |" or $self->throw ("Can't open 'du' pipe ($!)\n");
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
    close DF || $self->throw("Error from 'df' : $!\n");
}
