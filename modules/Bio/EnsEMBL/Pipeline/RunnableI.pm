#
# Interface for running programs
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

Bio::EnsEMBL::Pipeline::RunnableI

=head1 SYNOPSIS

=head1 DESCRIPTION

Interface for running external programs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableI;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Analysis::Programs;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);

=head1 ABSTRACT METHODS

These methods need to be implemented in any
module which implements
C<Bio::EnsEMBL::Pipeline::RunnableI>.

=head2 run

    $self->run();

Actually runs the analysis programs.  If the
analysis has fails, it should throw an
exception.  It should also remove any temporary
files created (before throwing the exception!).

=head2 output

    @output = $self->output();

Return a list of objects created by the analysis
run (eg: C<Bio::EnsEMBL::FeaturePair> objects).

=cut

sub run {
    my ($self) = @_;

    $self->throw("run not implemented");
}

sub output {
    my ($self) = @_;

    $self->throw("output not implemented");
}

#########################
# Added by MAQ 
# functions used by Runnable modules replacing hp.pl functions
# These aren't really interfaces more to do with background jobs.
#########################

sub locate_executable {
    my ($self, $runnable) = @_;
    if ($runnable)
    {
        Bio::EnsEMBL::Analysis::Programs->import($runnable);
        return $Bio::EnsEMBL::Analysis::Programs::Program_Paths { $runnable };
    }
}

#disk io methods

#parsefile is a public method
sub parsefile {
    my ($self, $filename) = @_;
    $self->results($filename) if ($filename);
    $self->parse_results();
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

sub writefile {
    my ($self, $seqobj, $seqfilename) = @_;
    unless ($seqobj)
    {
        print "Writing sequence to ".$self->filename."\n";
        #create Bio::SeqIO object and save to file
        my $clone_out = Bio::SeqIO->new(-file => ">".$self->filename , '-format' => 'Fasta')
               or $self->throw("Can't create new Bio::SeqIO from ".$self->filename.":$!\n");
        $clone_out->write_seq($self->clone) 
                or $self->throw("Couldn't write to file ".$self->filename.":$!");
    }
    else
    {
        $seqfilename = 'filename' unless ($seqfilename);
        print "Writing sequence to ".$self->$seqfilename()."\n";
        #create Bio::SeqIO object and save to file
        my $clone_out = Bio::SeqIO->new(-file => ">".$self->$seqfilename(). '-format' => 'Fasta')
               or $self->throw("Can't create new Bio::SeqIO from ".$self->$seqfilename().":$!\n");
        $clone_out->write_seq($self->$seqobj()) 
                or $self->throw("Couldn't write to file ".$self->$seqfilename().":$!");
    
    }
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

# FeaturePair methods:
sub createfeaturepair {
    my ($self, $feat1, $feat2) = @_;
    #some contant strings
    
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => $feat2->{db},
                            -db_version      => $feat2->{db_version},
                            -program         => $feat2->{program},
                            -program_version => $feat2->{p_version},
                            -gff_source      => $feat2->{source},
                            -gff_feature     => $feat2->{primary});
    
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

sub growfplist {
    my ($self, $fp) =@_;    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
}
