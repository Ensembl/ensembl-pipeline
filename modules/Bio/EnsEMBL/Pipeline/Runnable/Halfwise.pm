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

Bio::EnsEMBL::Pipeline::Runnable::Halfwise

=head1 SYNOPSIS

    
=head1 DESCRIPTION

=head2 Methods:


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Halfwise;

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
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Halfwise->new (-CLONE => $seq);
    Function:   Initialises Halfwise object
    Returns :   a HAlfwise Object
    Args    :   A Bio::Seq object (-CLONE) 
                (Halfwise location and arguments optional)

=cut
sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
    $self->{_seqfeature} = undef;   #Bio::SeqFeature
    $self->{_clone} = undef;        #location of Bio::Seq object
    $self->{_halfwise} = undef;     #location of Genscan script
    $self->{_workdir} = undef;      #location of temp directory
    $self->{_filename} =undef;      #file to store Bio::Seq object
    $self->{_results}   =undef;     #file to store results of Halfwise
    $self->{_protected} =[];        #a list of file suffixes protected from deletion
    $self->{_arguments} =undef;     #parameters for halfwise
    my($clonefile, $halfwise, $parameters) 
            = $self->_rearrange(['CLONE', 'HALFW', 'ARGS'], @args);
    
    $self->clone($clonefile) if ($clonefile);
    if ($halfwise)      {$self->halfwise($halfwise) ;}
    else                {$self->halfwise('/usr/local/pubseq/scripts/halfwise'); }
    if ($parameters)    {$self->clone($parameters) ; }
    else                {$self->arguments('-init wing -pseudo -caceh -cut 25 -aln 200 -quiet'); }     
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
        $self->results($self->filename.'.halfwise');
    }
    return $self->{_clone};
}

sub filename {
    my ($self, $filename) = @_;
    $self->{_filename} = $filename if ($filename);
    return $self->{_filename};
}

sub results {
    my ($self, $results) = @_;
    $self->{_results} = $results if ($results);
    return $self->{_results};
}

=head2 protect
    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes
    
=cut
sub protect {
    my ($self, @filename) =@_;
    push (@{$self->{_protected}}, @filename) if (@filename);
    return @{$self->{_protected}};
}

=head2 halfwise
    Title   :   halfwise
    Usage   :   $obj->halfwise('/nfs/disk100/humpub/OSFbin/halfwise');
    Function:   Get/set method for the location of halfwise
    Args    :   File path (optional)
    
=cut
sub halfwise {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("halfwise not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{_halfwise} = $location ;
    }
    return $self->{_halfwise};
}

=head2 arguments
    Title   :   arguments
    Usage   :   $obj->arguments('-init wing -pseudo -caceh -cut 25 -aln 200 -quiet');
    Function:   Get/set method for halfwise arguments arguments
    Args    :   File path (optional)
    
=cut
sub arguments {
    my ($self, $args) = @_;
    if ($args)
    {
        $self->{_arguments} = $args ;
    }
    return $self->{_arguments};
}

=head2 workdir
    Title   :   workdir
    Usage   :   $obj->wordir('~humpub/temp');
    Function:   Get/set method for the location of a directory to contain temp files
    Args    :   File path (optional)
    
=cut
sub workdir {
    my ($self, $directory) = @_;
    if ($directory)
    {
        mkdir ($directory, '777') unless (-d $directory);
        $self->throw ("$directory doesn't exist\n") unless (-d $directory);
        $self->{_workdir} = $directory;
    }
    return $self->{_workdir};
}

sub growfplist {
    my ($self, $fp) =@_;    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
}

###########
# Analysis methods
##########
=head2 run
    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs halfwise and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut
sub run {
    my ($self, $dir) = @_;
    #check clone
    my $seq = $self->clone() || $self->throw("Clone required for Halfwise\n");
    #set directory if provided
    $self->workdir('halfwise') unless ($self->workdir($dir));
    $self->checkdir();
    #write sequence to file
    $self->writefile(); 
    #run genscan       
    $self->run_halfwise();
    #parse output and create features
    $self->parse_halfwise();
    $self->deletefiles();
}

sub run_halfwise {
    my ($self) = @_;
    print "Running halfwise on ".$self->filename."\n";
    open (OUTPUT, $self->halfwise.' '.$self->parameters.' '.$self->filename.'|')
            or $self->throw("Couldn't open pipe to halfwise: $!\n");  
    open (RESULTS, ">".$self->results)
            or $self->throw("Couldn't create file for halfwise results: $!\n");
    print RESULTS <OUTPUT>;
    close OUTPUT;
    close RESULTS;
}

=head2 parsefile
Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses halfwise ace output to give a set of feature pairs
    Returns :   none
    Args    :   optional filename

=cut
sub parsefile {
    my ($self, $filename) = @_;
    $self->results($filename) if ($filename);
    $self->parse_halfwise();
}

sub parse_halfwise {
    my ($self) = @_;
    #file is in ace format. 
    open (HALFWISE, "<".$self->results)
        or $self->throw ("Couldn't open file ".$self->results.": $!\n");
    if (<HALFWISE> =~ /#No hits found in the blast search/ )
    {
        print "No hit found with halfwise\n";
        return;
    }
    
    while (<HALFWISE>)
    {
        if (/#Complete Analysis/)
        {
            my (%feat1, %feat2);
            my @halfwise = <HALFWISE>;
            for (my $index=0; $index < scalar(@halfwise); 
                             $index++, $_ = $halfwise[$index])
            {
                #find second sequence declaration (Feature 2)
                if ($feat2 {name} && $_ =~ /$feat2{name}/)
                {
                    for (; ($index < scalar(@halfwise)) && ($_ =~ /\b\S/);
                            $index++, $_ = $halfwise[$index])
                    {
                        chomp; #remove trailing \n
                        if (/Method/)#ACE: method
                        {
                           my @element = split;
                           $feat1 {source} = $element[1];
                        }
                        elsif (/Database/)#ACE: database
                        {
                           my @element = split;
                           $feat2 {db} = $element[1];
                           $feat2 {db_version} = $element[2];
                        }
                        elsif (/CDS\b|Pseudogene/)#ACE: type
                        {
                           $feat2 {primary} = $_;
                        }
                        elsif (/CDS_predicted_by/)#ACE: CDS_predicted_by genwise score
                        {
                           my @element = split;
                           $feat2 {program} = $element[1];
                           $feat1{score} = $element[2];
                           $feat2{score} = $element[2];
                        }
                        elsif (/Source_Exons/)
                        {
                           my @element = split;
                           $feat2{start} = $element[1];
                           $feat2{end} =  $element[2];
                        }
                        $feat1 {source} = 'Halfwise';
                        $feat1 {primary} = 'similarity';
                        $feat2 {source} = 'Halfwise';
                    } 
                    #Should have fully defined feat1 and feat2 at his point 
                    $self->createfeaturepair(\%feat1, \%feat2); 
                }               
                #find first sequence declaration (Feature 1)
                elsif (/\bSequence/)
                {
                    chomp;
                    my @element = split;
                    $feat1 {name} = $element[1];
                    $feat1 {strand} = 1;
                    for (; ($index < scalar(@halfwise)) && ($_ =~ /\S/ );
                            $index++, $_ = $halfwise[$index])
                    {
                        
                        chomp;
                        if (/subsequence/i)
                        {
                            my @element = split;
                            $feat2 {name} = $element[1];
                            #check orientation
                            if ($element[2] < $element [3])
                            {
                                $feat1 {start} = $element[2];
                                $feat1 {end} = $element[3];
                                $feat2 {strand} = 1;
                            }
                            else
                            {
                                $feat1 {start} = $element[3];
                                $feat1 {end} = $element[2];
                                $feat2 {strand} = -1;
                            }
                        }
                    }
                }
                #search next line    
            }
        }     
    }
}

sub createfeaturepair {
    my ($self, $feat1, $feat2) = @_;
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => $feat2->{db},
                            -db_version      => $feat2->{db_version},
                            -program         => $feat1->{source},
                            -program_version => 'unknown',
                            -gff_source      => $feat1->{source},
                            -gff_feature     => $feat1->{primary});
    
    #create and fill Bio::EnsEMBL::Seqfeature objects
    my $seqfeature1 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat1->{name},
                            -start   => $feat1->{start},
                            -end     => $feat1->{end},
                            -strand  => $feat1->{strand},
                            -score   => $feat1->{score},
                            -source_tag  => $feat1->{source},
                            -primary_tag => $feat1->{primary},
                            -analysis => $analysis_obj);
    
    my $seqfeature2 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname => $feat2->{name},
                            -start   => $feat2->{start},
                            -end     => $feat2->{end},
                            -strand  => $feat2->{strand},
                            -score   => $feat2->{score},
                            -source_tag  => $feat2->{source},
                            -primary_tag => $feat2->{primary},
                            -analysis => $analysis_obj);
    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $seqfeature1,
                                             -feature2 => $seqfeature2) ;
    $self->growfplist($fp);                             
}

##############
# input/output methods
#############
=head2 output
    Title   :   output
    Usage   :   obj->output()
    Function:   A list of Feature pairs 
    Returns :   A list of Feature pairs
    Args    :   none

=cut
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
    #delete all analysis files 
    my @list = glob($self->filename."*");
    foreach my $result (@list)
    {
        my $protected = undef; #flag for match found in $protected
        foreach my $suffix ($self->protect)
        {        
            $protected = 'true' if ($result eq $self->filename.$suffix);
        }
        unless ($protected)
        {
            unlink ($result) or $self->throw ("Couldn't delete $result :$!");    
        }
    }
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
