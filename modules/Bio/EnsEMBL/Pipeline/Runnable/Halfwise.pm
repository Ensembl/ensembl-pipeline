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

  #create and fill Bio::Seq object
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $seqfile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();

  #create Bio::EnsEMBL::Pipeline::Runnable::Genscan object
  my $halfwise = Bio::EnsEMBL::Pipeline::Runnable::Halfwise->new (-QUERY => $seq);
  $halfwise->workdir($workdir);
  #$halfwise->run();
  $halfwise->parsefile('/nfs/disk65/mq2/temp/halfwise.output');
  my $featurepairs = $halfwise->output();

=head1 DESCRIPTION

Takes a Bio::Seq (or Bio::PrimarySeq) object and runs halfwise on it. The results are returned as an
array of Bio::EnsEMBL::FeaturePairs. The location of halfwise can be set using
the halfwise() method

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
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Halfwise->new (-QUERY => $seq);
    Function:   Initialises Halfwise object
    Returns :   a HAlfwise Object
    Args    :   A Bio::Seq object (-QUERY) 
                (Halfwise location and arguments optional)

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    $self->{'_seqfeature'} = undef;  # Bio::SeqFeature
    $self->{'_query'}      = undef;  # location of Bio::Seq object
    $self->{'_halfwise'}   = undef;  # location of Genscan script
    $self->{'_workdir'}    = undef;  # location of temp directory
    $self->{'_filename'}   = undef;  # file to store Bio::Seq object
    $self->{'_results'}    = undef;  # file to store results of Halfwise
    $self->{'_protected'}  = [];     # a list of file suffixes protected 
                                     # from deletion
    $self->{'_arguments'}    = undef;     #parameters for halfwise
    my($query, $halfwise, $arguments) 
	= $self->_rearrange([qw(QUERY HALFW ARGS)], @args);
    
    $self->query($query) if ($query);
    if ($halfwise)  {$self->halfwise($halfwise) ;}
    else                
    {   
        eval 
        { $self->halfwise($self->locate_executable('halfwise')); };
        if ($@)
        {  $self->halfwise ('/usr/local/pubseq/scripts/halfwise'); }
    }
    if ($arguments)    
    {$self->arguments($arguments) ; }
    else                
    {$self->arguments('-init wing -pseudo -caceh -cut 25 -aln 200 -quiet'); }     
    return $self; 
}


###################
#get/set methods
###################

sub query {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeq") || $seq->isa("Bio::Seq")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_query'} = $seq ;
        
        $self->queryname($self->query->id);
        $self->filename($self->query->id.".$$.seq");
        $self->results($self->filename.'.halfwise');
    }
    return $self->{'_query'};
}

=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

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

=head2 queryname

    Title   :   queryname
    Usage   :   $obj->queryname('AC00074');
    Function:   Get/set method for query name. 
                This must be set manually when a file or pipe is parsed and the queryname is 
                not present in the executable output
    Args    :   File suffixes

=cut

=head2 workdir

    Title   :   workdir
    Usage   :   $obj->wordir('~humpub/temp');
    Function:   Get/set method for the location of a directory to contain temp files
    Args    :   File path (optional)

=cut

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
    my ($self) = @_;
    #check seq
    my $seq = $self->query() || $self->throw("Seq required for Halfwise\n");
    #set directory if provided
    $self->workdir('/tmp') unless $self->workdir();
    $self->checkdir();
    #write sequence to file
    $self->writefile();
    #run genscan       
    $self->run_halfwise();
    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
}

sub run_halfwise {
    my ($self) = @_;
    print "Running halfwise on ".$self->filename."\n";
    open (OUTPUT, $self->halfwise.' '.$self->filename.' '.$self->arguments.'|')
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
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut

#New and improved! takes filenames and handles, therefore pipe compliant!
sub parse_results {
    my ($self) = @_;
    my $filehandle;
    if (ref ($self->results) !~ /GLOB/)
    {
        #file is in ace format.
        open (HALFWISE, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*HALFWISE;
    }
    else
    {
        $filehandle = $self->results;
    }
    
    if (<$filehandle> =~ /#No hits found in the blast search/ )
    {
        print "No hit found with halfwise\n";
        return;
    }
    
    while (<$filehandle>)
    {
        if (/#Complete Analysis/)
        {
            my (%feat1, %feat2);
            my @halfwise = <$filehandle>;
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
                        $feat2 {primary} = 'similarity';
                        $feat2 {db} = 'unknown';
                        $feat2 {db_version} = 'unknown';
                    } 
                    #Should have fully defined feat1 and feat2 at his point 
                    $self->create_FeaturePair(\%feat1, \%feat2); 
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

##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   returns a list of Feature pairs 
    Returns :   A list of Feature pairs
    Args    :   none

=cut

sub output {
my ($self) = @_;
    return @{$self->{'_fplist'}};
}
