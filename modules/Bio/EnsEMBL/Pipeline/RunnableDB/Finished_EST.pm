
### Bio::EnsEMBL::Pipeline::RunnableDB::Finished_EST

package Bio::EnsEMBL::Pipeline::RunnableDB::Finished_EST;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_EST;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::PrimarySeq;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ( $new, @args ) = @_;
    my $self = $new->SUPER::new(@args);
    
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)
    #print STDERR Dumper($self);
    return $self;
}

sub fetch_input {
    my ($self) = @_;
    
    my @fps;
    
    $self->throw("No input id") unless defined( $self->input_id );
    
    my $contigid = $self->input_id;
    my $rawContigAdaptor = $self->db->get_RawContigAdaptor();
    my $contig   = $rawContigAdaptor->fetch_by_name($contigid);
    
    my $masked   = $contig->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING,$SOFT_MASKING) or $self->throw("Unable to fetch contig");
    # make a Bio::PrimarySeq obj to remove some db intensive queries later on
    my $unmasked = Bio::PrimarySeq->new(-display_id => $contig->display_id(),
					-id         => $contig->id(),
					-seq        => $contig->seq()
					);
    my $seq = $masked->seq;
    if( scalar($seq =~ s/([CATG])/$1/g) > 3 ){
        $self->input_is_void(0);
        $self->check_with_seg($masked);
    }else{
        $self->input_is_void(1);
        $self->warn("Need at least 3 nucleotides");
    }

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Finished_EST->new(
        '-query'      => $masked,
        '-unmasked'   => $unmasked,
        '-analysis'   => $self->analysis,
    );
    $self->runnable($runnable);
}
sub check_with_seg{
    my ($self, $seqObj_to_test) = @_;

    warn "need a Bio::Seq Obj" unless $seqObj_to_test;

    my ($filename) = $self->_createfiles('/tmp',[qw(seg_checking)]);
    my $file = Bio::SeqIO->new(-file   => ">$filename", 
                               -format => 'Fasta') 
        or $self->throw("Can't create Bio::SeqIO $filename $!");
    $file->write_seq($seqObj_to_test);

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
sub runnable {
    my ( $self, $arg ) = @_;
    
    if ( defined($arg) ) {
        $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI")
        unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");
        
        $self->{_runnable} = $arg;
    }
    
    return $self->{_runnable};
}

sub run {
    my ($self) = @_;
    
    my $runnable = $self->runnable;
    $runnable || $self->throw("Can't run - no runnable object");
    
    $runnable->run;
    
    my $db_version = $runnable->db_version_searched if $runnable->can('db_version_searched');
    $self->db_version_searched($db_version);
    if ( my @output = $runnable->output ) {
	        my $dbobj      = $self->db;
	        my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
	        my %ids        = map { $_->hseqname, 1 } @output;
	        $seqfetcher->write_descriptions( $dbobj, keys(%ids) );        
    }
    
    return 1;
    
}
sub db_version_searched{
    my $self = shift;
    $self->{'_db_version_searched'} = shift if @_;
    return $self->{'_db_version_searched'};
}


sub output {
    my ($self) = @_;
    my @runnable = $self->runnable;
    my @results;
    foreach my $runnable (@runnable) {
        print STDERR "runnable = " . $runnable[0] . "\n";
        push ( @results, $runnable->output );
    }
    return @results;
}
1;

=head1 NAME - Bio::EnsEMBL::Pipeline::RunnableDB::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

