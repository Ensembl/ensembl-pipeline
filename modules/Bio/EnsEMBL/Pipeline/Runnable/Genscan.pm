

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
  my @genes = $genscan->output();
  my @exons = $genscan->output_exons();
  my $seqfeature = $genscan->output_singlefeature();

=head1 DESCRIPTION

This package is based on Genscan2ace.
Genscan takes a Bio::Seq (or Bio::PrimarySeq) object and runs Genscan on it. The
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

=head1 SEE ALSO

=over 4

=item B<Bio::EnsEMBL::Pipeline::RunnableI>

=item B<Bio::EnsEMBL::Pipeline::RunnableDB::Genscan> 

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
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis; 
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::Utils;
use Bio::Root::RootI;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Genscan->new (-CLONE => $seq);
    Function:   Initialises Genscan object
    Returns :   a Genscan Object
    Args    :   A Bio::Seq object 
                (Genscan location and matrix file location optional)

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    #$self->{_seqfeature} = undef;   #Bio::SeqFeature container for Subseqfeatures
    $self->{'_exons'}     = [];    # an array of Bio::Seqfeatures (exons)
    $self->{'_genes'}     = [];    # an array of arrays of SeqFeatures
    #$self->{_feature_data} = [];  # an array of data for SeqFeatures
    $self->{'_peptides'}  = [];    # genscan predicted peptide (used for phase)
    $self->{'_clone'}     = undef; # location of Bio::Seq object
    $self->{'_genscan'}   = undef; # location of Genscan script
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of genscan
    $self->{'_protected'} = [];    # a list of file suffixes protected from deletion
    $self->{'_parameters'} =undef; #location of parameters for genscan

    my($clone, $genscan, $parameters, $matrix) = 
        $self->_rearrange([qw(CLONE GENSCAN PARAM MATRIX)], @args);


    $self->clone($clone);
    $genscan = 'genscan'       unless ($genscan);
    $matrix  = 'HumanIso.smat' unless ($matrix);

    $self->genscan($self->find_executable($genscan));
    $self->matrix ($self->find_file($matrix));


    if ($parameters)    
    { $self->parameters($parameters) ; }
    else                
    {$self->parameters(''); }     

    return $self;
}


###################
#get/set methods
###################

sub clone {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_clone'} = $seq ;
        $self->filename($self->clone->id.".$$.seq");
        $self->results($self->filename.".genscan");
    }
    return $self->{'_clone'};
}


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
        $self->{'_genscan'} = $location ;
    }
    return $self->{'_genscan'};
}

=head2 matrix

    Title   :   matrix
    Usage   :   $obj->matrix('/nfs/disk100/humpub/OSFbin/HumanIso.smat');
    Function:   Get/set method for the location of genscan matrix
    Args    :   File path (optional)

=cut

sub matrix {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Genscan matrix not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_matrix'} = $location ;
    }
    return $self->{'_matrix'};
}


=head2 parameters

    Title   :   parameters
    Usage   :   $obj->parameters('parameters');
    Function:   Get/set method for the location of genscan parameters
    Args    :   File path (optional)

=cut

sub parameters {
    my ($self, $param) = @_;
    if ($param)
    {
        $self->{'_parameters'} = $param;
    }
    return $self->{'_parameters'};
}

sub exons {
    my ($self, $exon) =@_;
    if ($exon)
    {
        $exon->isa("Bio::EnsEMBL::SeqFeature") 
                || $self->throw("Input isn't a Bio::EnsEMBL::SeqFeature");
        push(@{$self->{'_exons'}}, $exon);
    }
    return @{$self->{'_exons'}};
}

#empties exon array after they have been converted to genes
sub clear_exons {
    my ($self) = @_;
    $self->{'_exons'} = [];
}

sub genscan_genes {
    my ($self, $gene) =@_;
    if ($gene)
    {
        $gene->isa("Bio::EnsEMBL::SeqFeature") 
                || $self->throw("Input isn't a Bio::EnsEMBL::SeqFeature");
        push(@{$self->{'_genes'}}, $gene);
        @{$self->{'_genes'}} = sort { $a->seqname <=> $b->seqname } @{$self->{'_genes'}};
    }
    return @{$self->{'_genes'}};
}

sub genscan_peptides {
    my ($self, $peptide) = @_;
    push (@{$self->{'_peptides'}}, $peptide) if ($peptide);
    return @{$self->{'_peptides'}};
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

    my ($self) = @_;
    #check seq
    my $seq = $self->clone() || $self->throw("Seq required for Genscan\n");
    #set directory if provided
    $self->workdir('/tmp') unless $self->workdir();
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
    #print "Running genscan on ".$self->filename."\n";
    system ($self->genscan.' '.$self->matrix.' '.$self->filename.' > '.$self->results);
    $self->throw($self->results." not created by Genscan\n") unless (-e $self->results);
}

=head2 parsefile

    Title   :  parsefile
    Usage   :   $obj->parsefile($filename)
    Function:   Parses Genscan output to give a set of seqfeatures
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
                NOTE: genscan can not assign phases to exons from the output
                file unless the sequence is supplied as a Bio::Seq object.
    Returns :   none
    Args    :   optional filename

=cut

sub parse_results {
    my ($self) = @_;

    my %exon_type = ('Sngl', 'Single Exon',
                     'Init', 'Initial Exon',
                     'Intr', 'Internal Exon',
                     'Term', 'Terminal Exon');

    my $filehandle;
    if (ref ($self->results) !~ /GLOB/)
    {
        open (GENSCAN, "<".$self->results)
            or $self->throw ("Couldn't open file ".$self->results.": $!\n");
        $filehandle = \*GENSCAN;
    }
    else
    {
        $filehandle = $self->results;
    }

    if (<$filehandle> =~ m|NO EXONS/GENES PREDICTED IN SEQUENCE|i )
    {
        print STDERR "No genes predicted\n";
        return;
    }

    #The big parsing loop - parses exons and predicted peptides
    while (<$filehandle>)
    {
        # Last line before predictions contains nothing but spaces and dashes
        if (/^\s*-[-\s]+$/) 
        { 
            while (<$filehandle>) 
            {
                my %feature; 
                if (/init|term|sngl|intr/i)
                {
                   my @element = split;
                   $self->throw("Unable to parse Genscan ouput (".scalar(@element).") Line: $_\n") 
                                unless (scalar(@element) == 13); 
                                
                   my ($gene, $exon) = split (/\./, $element[0]); 

                   $feature {name} = $gene + ($exon/1000); #name must be a number
                   #arrange numbers so that start is always < end
                   if ($element[2] eq '+')
                   {
                        $feature {'start'} = $element[3];
                        $feature {'end'} = $element[4];
                        $feature {'strand'} = 1;
                   }
                   elsif ($element[2] eq '-')
                   {
                        $feature {'start'} = $element[4];
                        $feature {'end'} = $element[3];
                        $feature {'strand'} = -1;
                   }
                   $feature {'score'} = $element[12];
                   $feature {'p'}     = $element[11];
                   $feature {'type'} = $exon_type{$element[1]};
                   $feature {'program'} = 'Genscan';
                   $feature {'program_version'} = '1.0';
                   $feature {'primary'} = 'prediction';
                   $feature {'source'} = 'genscan';

                   $self->create_feature(\%feature);
                }
                #gene/exon data ends with line 'predicted peptide sequence(s)'
                elsif (/predicted peptide/i)
                {   
                    last;   
                }
            }
            #begin parsing of peptide data
            my $peptide_string;
            my $reading_peptide =0;
            while (<$filehandle>)
            {
                if (/\A[\W]\Z/)
                {
                    $self->genscan_peptides($peptide_string);
                    $peptide_string = undef;
                    $reading_peptide = 0;
                }
                elsif ($reading_peptide)
                {
                    $_ =~ s/\n//g;
                    $peptide_string .= $_;
                }
                elsif (/^>/)
                {
                    $reading_peptide = 1;
                }
            }
            $self->genscan_peptides($peptide_string); #final peptide 
        }
    }
    #end of big loop. Now build up genes
    $self->create_genes();

    unless ($self->clone)
    {
        print STDERR "Can't calculate phases if Bio::Seq isn't supplied\n";
        return;
    }
    $self->calculate_and_set_phases_new();
    $self->clear_exons(); #free up unecessary storage 
}


sub calculate_and_set_phases_new {
    my ($self) = @_;

    my @genes       = $self->genscan_genes();
    my @peptides    = $self->genscan_peptides();

    $self->throw("Mismatch in number of genes (".scalar(@genes).
                 ") and peptides ("             .scalar(@peptides).
                 ") parsed from file") unless (scalar(@genes) == scalar (@peptides));

    my $i = 0;
    my $count = 1;
  GENE: while ($i < scalar(@genes)) {

#      print STDERR "Gene " . $genes[$i]->seqname . "\n";

      my $peptide = $peptides[$i];
      my @exons   = $genes[$i]->sub_SeqFeature();
#      print STDERR "Exons are $#exons\n";
      my @newtran = Bio::EnsEMBL::DBSQL::Utils::fset2transcript_3frame($genes[$i],$self->clone);

#      print STDERR "\nPeptide is " . $peptides[$i] . "\n";

      my $translation_found = 0;

      foreach my $tran (@newtran) {

        my $temp_tran = $tran->translate->seq;

#       print STDERR "Translation is " . $temp->tran . "\n";

        # clean the translated sequence

        #genscan translated partial genes correctly whilst exon translation begin with M
        #$temp_tran =~ s/^M//i; #remove initial M from exon

        # remove any initial X's from the translation
        $temp_tran =~ s/^x//i;
        
        # remove any terminal X's from the translation
        $temp_tran =~ s/x$//i;
	my $genscan = $peptides[$i];
	my $x = 0;
	while (($x = index($temp_tran, 'X', $x)) != -1) {
	  #print STDERR "Found an 'X' at ", $i + 1, "\n";
	  substr($genscan, $x, 1) = 'X';
	  $x++;
	}
    
	#print "\nafter\ngenscan: $genscan\nensembl: $ensembl\n";
        if (index($genscan ,$temp_tran) >= 0) {
#         print STDERR $tran->temporary_id . " " . $tran->translate->seq . "\n";

          $translation_found = 1;
          foreach my $exon ($tran->get_all_Exons) {
           #print "exon ".$exon->seqname." ".$exon->start . " " . $exon->end . " " . $exon->phase . " " . $exon->end_phase . " " .$exon->strand . "\n";
          }
          $tran->temporary_id($self->clone->id . "." . $count);
          $count++;
	  #print "translation = ".$temp_tran."\n";
	  if($temp_tran =~/\*/){
	    #print "translation ".$tran->temporary_id." = ".$temp_tran."\n";
	    #print "transcript contains stop codons!";
	  } else {
	   #print "adding translation ".$tran->temporary_id." \n";
	    $self->add_Transcript($tran);
	  }
	  $i++;
          next GENE;
        }
      }
      
      unless ($translation_found) {
        $self->throw("[Genscan.pm] Unable to match Genscan peptide ".$peptides[$i].
                     " in a translation\n");
      }

#      print "\n";
      $i++;
    }
  }

sub add_Transcript {
  my ($self,$transcript) = @_;

  if (defined($transcript)) {
    if (!defined($self->{_transcripts})) {
      $self->{_transcripts} = [];
    }
    push(@{$self->{_transcripts}},$transcript);
  }
}

sub each_Transcript {
  my ($self) = @_;

  if (!defined($self->{_transcripts})) {
    $self->{_transcripts} = [];
  }
  return @{$self->{_transcripts}};
}


sub create_feature {
    my ($self, $feat) = @_;
    #$self->create_gene();

    #create analysis object
    my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => $feat->{'program'},
                            -program_version => $feat->{'program_version'},
                            -gff_source      => $feat->{'source'},
                            -gff_feature     => $feat->{'primary'});

    #create and fill Bio::EnsEMBL::Seqfeature objects   
    my $exon = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname => $feat->{'name'},
                            -start   => $feat->{'start'},
                            -end     => $feat->{'end'},
                            -strand  => $feat->{'strand'},
                            -score   => $feat->{'score'},
                            -frame   => $feat->{'frame'},
                            -p_value => $feat->{'p'},
                            -source_tag  => $feat->{'source'},
                            -primary_tag => $feat->{'type'},
                            -analysis => $analysis_obj);  
    $self->exons($exon);
}

#creates groups of exons as subseqfeatures.
#relies on seqname of exons being in genscan format 3.01, 3.02 etc
sub create_genes {
    my ($self) = @_;
    #print "creating genes \n";
    my (%genes, %gene_start, %gene_end, %gene_score, %gene_p,
        %gene_strand, %gene_source, %gene_primary, %gene_analysis);

    my @ordered_exons = sort { $a->seqname <=> $b->seqname } $self->exons();

    #sort exons into hash by initial numbers of seqname (genes)
    foreach my $exon (@ordered_exons)
    {
        my ($group_number) = ($exon->seqname =~ /(\d+)\./);
	#print "seqname =  ".$exon->seqname."\n";
        #intialise values for new gene
        unless (defined ($genes {$group_number}))
        {
            $genes          {$group_number} = [];
            $gene_start     {$group_number} = $exon->start;
            $gene_end       {$group_number} = $exon->end;
            $gene_score     {$group_number} = 0 ;
            $gene_strand    {$group_number} = $exon->strand;
            $gene_source    {$group_number} = $exon->source_tag ;
            $gene_primary   {$group_number} = $exon->primary_tag;
            $gene_analysis  {$group_number} = $exon->analysis;
            $gene_p         {$group_number} = 0 ;
        }
        #fill array of exons
        push (@{$genes {$group_number}}, $exon);
        #calculate gene boundaries and total score
        $gene_start {$group_number} = $exon->start() 
            if ($exon->start() < $gene_start{$group_number});
        $gene_end   {$group_number} = $exon->end() 
            if ($exon->end() > $gene_end{$group_number});
        $gene_score {$group_number} += $exon->score();
        $gene_p     {$group_number} += $exon->p_value();
    }

    #create Bio::SeqFeature objects (genes) with SubSeqFeatures (exons)
    foreach my $gene_number (keys(%genes))
    {
        my $gene = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname     => $gene_number,
                            -strand      => $gene_strand   {$gene_number},
                            -score       => $gene_score    {$gene_number}
                                            /(scalar @{$genes {$gene_number}}),
                            -p_value     => $gene_p         {$gene_number}
                                            /(scalar @{$genes {$gene_number}}),
                            -start       => $gene_start    {$gene_number},
                            -end         => $gene_end      {$gene_number},
                            -source_tag  => $gene_source   {$gene_number},
                            -primary_tag => $gene_primary  {$gene_number},
                            -analysis    => $gene_analysis {$gene_number}, )
                    or $self->throw("Couldn't create Bio::EnsEMBL::SeqFeature object");

        foreach my $exon (@{$genes {$gene_number}})
        {
          $gene->add_sub_SeqFeature($exon, '');
        }
        $self->genscan_genes($gene); #add gene to main object
      }
}

##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of SeqFeatures representing predicted genes 
                with exons stored as SubSeqFeatures.
    Returns :   An array of SeqFeatures (genes) containing sub-seqfeatures (exons)
    Args    :   none

=cut

sub output {

  my ($self) = @_;

  my @feat;

  my $analysis = Bio::EnsEMBL::Analysis->new(   -db              => undef,
                                                -db_version      => undef,
                                                -program         => 'genscan',
                                                -program_version => 1,
                                                -gff_source      => 'genscan',
                                                -gff_feature     => 'prediction',
                                                -logic_name      => 'genscan',
                                                );

  
  foreach my $transcript ($self->each_Transcript) {
    my @exons = $transcript->get_all_Exons;

    if ($exons[0]->strand == 1) {
      @exons = sort {$a->start <=> $b->start } @exons;
    } else {
      @exons = sort {$b->start <=> $a->start } @exons;
    }
#    print STDERR "\n" .$transcript->temporary_id . "\n";
#    print "\ntranscript ".$transcript->temporary_id." translates to ".$transcript->translate->seq."\n\n";
    foreach my $exon (@exons) {
      my $f = new Bio::EnsEMBL::SeqFeature(-seqname => $self->clone->id.".".$exon->seqname,
                                           -start   => $exon->start,
                                           -end     => $exon->end,
                                           -strand  => $exon->strand,
                                           -phase   => $exon->phase,
                                           -end_phase => $exon->end_phase,
                                           -score   => $exon->score,
                                           -p_value   => $exon->p_value,
                                           -source_tag => 'genscan',
                                           -primary_tag => 'prediction',
                                           -analysis     => $analysis);
      my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $self->clone->id.".".$exon->seqname,
                                            -start   => $exon->start,
                                            -end     => $exon->end,
                                            -strand  => $exon->strand,
                                            -phase   => $exon->phase,
                                            -end_phase => $exon->end_phase,
                                            -score   => $exon->score,
                                            -p_value   => $exon->p_value,
                                            -source_tag => 'genscan',
                                            -primary_tag => 'prediction',
                                            -analysis     => $analysis);

#      print STDERR $exon->start . " " . $exon->end . " " . $exon->phase . " " . $exon->strand . "\n";

      $f->analysis($analysis);
      $f2->analysis($analysis);

      my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
                                             -feature2 => $f2);

      push(@feat,$fp);
    }

  }
  return @feat;
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
    my @exons;

    foreach my $gene ($self->genscan_genes)
    {
        push (@exons, $gene->sub_SeqFeature);
    }
    print STDERR "No exons predicted\n" unless (@exons);
    @exons = sort { $a->seqname <=> $b->seqname } @exons;
    return @exons;
}

=head2 output_singlefeature

    Title   :   output_singlefeature
    Usage   :   obj->output_singlefeature()
    Function:   Returns a single SeqFeature with exons as sub-SeqFeatures
    Returns :   A single SeqFeature with exons as sub-SeqFeatures
    Args    :   none

=cut

sub output_singlefeature {
    my ($self) = @_;
    my ($start, $end, $score, $analysis, $primary, $source, $p_value);

    my (@genes) = $self->genes();
    print STDERR "No exons predicted\n" unless (@genes);
    #calculate boundaries and aggregate values
    foreach my $gene (@genes)
    {
        $start      =  $gene->start()  if (!defined($start) || $gene->start() < $start);
        $end        =  $gene->end()    if (!defined($end)   || $gene->end()   > $end);
        $score      += $gene->score();
        $p_value    += $gene->p_value();
        $analysis   =  $gene->analysis();
        $primary    =  $gene->primary_tag();
        $source     =  $gene->source_tag();
    }
    $score  = $score/scalar(@genes); #average score

    my $single = Bio::EnsEMBL::SeqFeature->new
                        (   -seqname        => 'genscan',
                            -strand         => 1,
                            -score          => $score,
                            -start          => $start,
                            -end            => $end,
                            -source_tag     => $source,
                            -primary_tag    => $primary,
                            -p_value        => $p_value,
                            -analysis       => $analysis )
                    or $self->throw("Couldn't create Bio::EnsEMBL::SeqFeature object");

    foreach my $gene (@genes)
    {
        foreach my $exon ($gene->sub_SeqFeature())
        {
            $single->add_sub_SeqFeature($exon, '');
        }    
    }
    return $single;
}

1;
