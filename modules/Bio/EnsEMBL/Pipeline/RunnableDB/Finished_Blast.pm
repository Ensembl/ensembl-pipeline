#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Blast

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $blast   = Bio::EnsEMBL::Pipeline::RunnableDB::Blast->new ( 
                                                    -db         => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
$blast->fetch_input();
$blast->run();
$blast->output();
$blast->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Blast to add
functionality for reading and writing to databases.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Finished_Blast;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Pipeline::Config::Blast;
use Bio::EnsEMBL::Pipeline::Config::General;
use Hum::Sort 'ace_sort';
use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

my %UNGAPPED;
my %UNMASKED;
my %NO_DESC;

foreach my $db (@$DB_CONFIG) {
    my (   $name,         $ungapped,         $unmasked,       $pfetch)
    = ($db->{'name'}, $db->{'ungapped'}, $db->{min_unmasked}, $db->{no_desc});
    
    if($db && $name){
        $UNGAPPED{$name} = $ungapped;
        $UNMASKED{$name} = $unmasked;
        $NO_DESC{$name}  = $pfetch;
    }else{
        my($p, $f, $l) = caller;
        warn("either db ".$db." or name ".$name." isn't defined so can't work $f:$l\n");
    }
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
    my $genseq;
   
    $self->throw("No input id") unless defined($self->input_id);
    print STDERR "INPUT ID: " . $self->input_id . "\n";    

    my $contig = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
    $genseq    = $contig->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING,$SOFT_MASKING) or $self->throw("Unable to fetch contig");
    $self->query($genseq || $contig);
    
    my $seq = $self->query->seq;
    
    if ($seq =~ /[CATG]{3}/) {
        $self->input_is_void(0);
        #$self->check_with_seg($self->query);
    } else {
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }
    
    my $ungapped;
    if($UNGAPPED{$self->analysis->db_file}){
        $ungapped = 1;
    } else {
        $ungapped = undef;
    }
    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast->new(-query          => $self->query,
        -database       => $self->analysis->db_file,
        -program        => $self->analysis->program,
        -options        => $self->analysis->parameters,
        -threshold_type => 'PVALUE',
        -threshold      => 1,
        -ungapped       => $ungapped,
    );
    
    $self->runnable($runnable);

    return 1;
}

sub check_with_seg{
    my ($self, $seqObj_to_test) = @_;

    warn "need a Bio::Seq Obj" unless $seqObj_to_test;

    my ($filename) = $self->_createfiles('/tmp',[qw(seg_checking)]);
    my $file = Bio::SeqIO->new(-file   => ">$filename", 
                               -format => 'Fasta') 
        or $self->throw("Can't create Bio::SeqIO $filename $!");
    my $translated = $seqObj_to_test->translate;
    #warn "************************************************************";
    #warn $translated->seq;
    #warn "************************************************************";
    $file->write_seq($translated);
    
    my $seg_cmd = "nseg $filename -x";
    my $seg = Bio::SeqIO->new(-file   => "$seg_cmd |",
                              -format => 'Fasta')
        or $self->throw("Can't create Bio::SeqIO $seg_cmd $!");
    my $seq;
    eval{
        $seq = $seg->next_seq->seq;
    };
    unlink($filename);
    if($@){
        $self->throw("There was a problem with SEG masking.\nI tried to '$seg_cmd'");
    }
    #warn "************************************************************";
    #warn $seq;
    #warn "************************************************************";
    if($seq =~ /[CATG]{3}/i){
        $self->input_is_void(0);
    }else{
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides after SEG filtering");
    }
    
}
sub _createfiles {
    my ($self, $dirname, $filenames) = @_;
    
    my $unique = {};
    $unique    = { map { $_, $unique->{$_}++ } @$filenames };
    my @files  = ();

    $dirname ||= '/tmp';
    $dirname   =~ s!(\S+)/$!$1!;

    foreach my $file(@$filenames){
        if($unique->{$file}){
            #name not unique add random
            $file .= ".$$.".int(rand(200));
            push(@files, "$dirname/$file");
        }else{
            #name was unique just add it
            push(@files, "$dirname/$file.$$");
        }
    }

    return @files;
}
=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::xxxx->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;
    
    foreach my $runnable ($self->runnable) {
        
        $self->throw("Runnable module not set") unless ($runnable);
        
        # Not sure about this
        $self->throw("Input not fetched")       unless ($self->query);
        eval{
            $runnable->run();
        };
        if($@){
            chomp $@;
            $self->failing_job_status($1) 
                if $@ =~ /^\"([A-Z_]{1,40})\"$/i; # only match '"ABC_DEFGH"' and not all possible throws
            $self->throw("$@");
        }
        my $db_version = $runnable->get_db_version if $runnable->can('get_db_version');
        $self->db_version_searched($db_version); # make sure we set this here
        # $self->updateAnalysis_with_db_version(); # alter the analysis.
        # trying to move this from Job.pm
        if ( (my @output = $runnable->output) && !($NO_DESC{$self->analysis->db_file})) {
            my $dbobj      = $self->db;
            my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
            my %ids        = map { $_->hseqname, 1 } @output;
            $seqfetcher->write_descriptions( $dbobj, keys(%ids) );
        }else{
            warn "Either didn't find anything or Don't fetch descriptions was true\n";
        }
    }
    return 1;
}

=head2 db_version_searched

    Title   :  db_version_searched
               [ distinguished from Runnable::*::get_db_version() ]
    Useage  :  $self->db_version_searched('version string')
               $obj->db_version_searched()
    Function:  Get/Set a blast database version that was searched
               The actual look up is done in Runnable::Finished_Blast
               This is just a holding place for the string in this
               module
    Returns :  String or undef
    Args    :  String
    Caller  :  $self::run()
               Job::run_module()

=cut

sub db_version_searched{
    my ($self, $arg) = @_;
    
    $self->{'_db_version_searched'} = $arg if $arg;

    return $self->{'_db_version_searched'};
}

sub updateAnalysis_with_db_version{
    my ($self) = @_;
    my $db_version  = $self->db_version_searched();
    $self->throw(" *** No db version *** " . 
                 " Probably BlastableVersion's fault." . 
                 " check out debugging that.") unless $db_version;
    my $analysis    = $self->analysis();
    my $current_dbv = $analysis->db_version();
    unless($current_dbv){
        # got to store it if nothing already there.
        # set the analysis to have the searched db_version
        $analysis->db_version($db_version); 
        # and store it
        $analysis->store();
        # and return
        return 1;
    }
    if ($db_version ne $current_dbv){
        # this temporarily sets the analysis obj to be
        # the db version of the searched blastable.
        # later on we should do a $analysis->store here
        $analysis->db_version($db_version);
        return 1;
    }

    # method to check most up to date not used
    # got to check which is the most up to date
    my @dbsearched  = ($db_version  =~ /(\d+)-(\w+)-(\d+)(.+)?/);
    my @current_dbv = ($current_dbv =~ /(\d+)-(\w+)-(\d+)(.+)?/);
    if(@dbsearched >= 3 && @current_dbv >= 3){
        # comparing 01-Jan-04
        my ($day, $mon, $yr, $ver)     = @dbsearched;
        my ($cday, $cmon, $cyr, $cver) = @current_dbv;
        $ver  =~ s/\D+//g if $ver;  # only interested in numbers
        $cver =~ s/\D+//g if $cver; # only interested in numbers

        my $lt = 0;
        # this means that the info had a version
        # and the current db version was lower than what it was searched
        $lt &= 8 if ($ver && $cver && $cver < $ver);
        $lt &= 1 if ($cyr  < $yr);
        $lt &= 2 if ($cmon < $mon);
        $lt &= 4 if ($cday < $day);

        return 1 unless($lt);  # current analysis->db_version is up-to-date
        
    }
}

1;
