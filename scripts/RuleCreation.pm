#RuleCreation is a package for creating ensembl-pipeline rules on the
#basis of a config file or writing those rules to a config file once
#retrived from a database. It has the same basis as AnalysisCreation
#and is based on an idea from Glenn Proctor
# Cared for by ensembl
#
# Copyright ensembl
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

RuleCreation

=head1 SYNOPSIS

This class will parse config files to produce rule objects and store 
them and will take rule tables and write a config file

[RepeatMask]
condition=SubmitContig

[Pmatch]
condition=SubmitChromosome

[Pmatch_Wait]
condition=Pmatch

[BestPmatch]
condition=Pmatch_Wait
condition=SubmitGenome
#comment lines can be made if they start with a # symbol
=head1 DESCRIPTION



=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

the class itself obviously doesn't' need to be instantiated but the 
either the script which uses it should be in the same directory as 
it or the directory which contains it should be in you PERL5LIB

the rule_setup script which should be found in the directory can
perform both functions for you if you have the appropriate database
and config files


an example config file can be found in this directory called 
example_rule.conf
=cut

package RuleCreation;

use strict;
use warnings;
use Bio::EnsEMBL::Pipeline::Rule;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(parse_files write_into_db read_db write_file 
                 create_rules);

verbose('WARNING');


=head2 parse_files

  Arg [1]   : array of filenames
  Function  : parse a rule config file and produce rule objects
  Returntype: an hash ref of rule goal/conditions
  Exceptions: if file doesn't exist
              if config format is in correct
              if key already exists for a particular header'
  Caller    : 
  Example   : my $rule_hash = parse_files($file);

=cut



sub parse_files {

  my @files = shift;
  if(@files == 0){
    throw("Can't parse files if we don't have any @files");
  }

  my %config;
  # read each file
  
  foreach my $file (@files) {
    
    if (! -e $file) {
      throw("rule file $file not found\n");
    }
    my $header = "";
    open (FILE, $file) or throw "Couldn't open file $file";
    while (<FILE>) {
      chomp();
      # Comment or blank line
      next if (/^\s$/ || /^\#/);
      # [HEADER]
      if (/^\[(.*)\]$/) {  # $1 will be the header name, without the []
        $header = $1;
        #print "Reading stanza $header\n";
        if(!$config{$header}){
          $config{$header} = [];
        }else{
          throw("Can't have two different rules with the same goal ".
                "header");
        }
      }
      # key=value
      if (/^([^=\s]+)\s*=\s*(.+)/) {   # $1 = key, $2 = value
        my $key = lc($1);  # keys stored as all lowercase, 
        #values have case preserved
        my $value = $2;
        if (length($header) == 0) {
          throw("Found key/value pair $key/$value outside stanza");
        }
        if($key eq 'condition'){
          push(@{$config{$header}}, $value);
        }else{
          warn("Don't recognise ".$key." so doing nothing with ".
               $value);
        }
      }
    } # while <FILE>
    close FILE;
  }
  foreach my $goal(keys(%config)){
    if($config{$goal} == 0){
      throw("Can't create a rule for ".$goal.
            " without any conditions");
    }
  }
  return \%config;
}


=head2 create_rules

  Arg [1]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Arg [2]   : hashref for hash keyed on goal analysis logic_name
  containing condition logic_names
  Function  : create rules on the basis of the hash
  Returntype: array ref of Bio::EnsEMBL::Pipeline::Rule's
  Exceptions: throws if not passed a dbadaptor or it is of the wrong 
  type or if not passed a rule hash or it is empty
  Example   : my @rules = @{create_rules($db, $rule_hash)}'

=cut


sub create_rules{
  my ($db, $config) = @_;

  if(!$db || !$db->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
    throw("Must define a db and ".$db." must be a ".
          "Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
  }
  if(!$config || keys(%$config) == 0){
    throw("Not going to create any rules as rule hash ".
          $config." is empty ");
  }
  my @rules;
  my $analysis_adaptor = $db->get_AnalysisAdaptor;
  my $rule_adaptor = $db->get_RuleAdaptor;
  foreach my $goal(keys(%$config)){
    my $analysis = $analysis_adaptor->fetch_by_logic_name($goal);
    if(!$analysis){
      throw("Can't have rule with goal ".$goal." if no analysis ".
            "of that name exists");
    }
    my @conditions = @{$config->{$goal}};
    my %checked;
    foreach my $condition(@conditions) {
      throw "$condition is not found in the database\n" 
        unless scalar($analysis_adaptor->fetch_by_logic_name($condition));
      $checked{$condition} = 1;
    }
    throw("Sorry, you can't insert a rule if the goal is the same as ".
          "one of the conditions") if exists $checked{$goal};
    my $rule = Bio::EnsEMBL::Pipeline::Rule->new
      (
       -goalanalysis => $analysis
      );
    foreach my $condition(keys %checked) {
      $rule->add_condition($condition);
    }
    push(@rules, $rule);
  }
  return \@rules;
}

=head2 write_into_db

  Arg [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
      [2]   : Ref to an array of rule objects
  Function  : Write the rule objects into the database
  Returntype: N/A
  Exceptions: if dbadaptor is the wrong type of object
  Caller    : 
  Example   : &write_into_db($db, \@rules);

=cut

sub write_into_db{
  my ($db, $rules) = @_;

  if(!$db || !($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
    throw("need a Pipeline::DBAdaptor not ".$db);
  }
  if(!$rules || @$rules == 0){
    throw("Need rules to write to the database");
  }
  my $ruleadaptor = $db->get_RuleAdaptor;
  my @existing_rules = $ruleadaptor->fetch_all;
  my %existing;
  foreach my $exist(@existing_rules){
    $existing{$exist->goalAnalysis->logic_name} = $exist;
  }
  foreach my $rule(@$rules){
    if(!$existing{$rule->goalAnalysis->logic_name}){
      $ruleadaptor->store($rule);
    }else{
      warn("Not storing rule ".$rule->goalAnalysis->logic_name.
           " rule with same goal already exists ".
           $existing{$rule->goalAnalysis->logic_name}->dbID);
    }
  }
}

=head2 read_db

  Arg [1]   : Bio::EnsEMBL::DBAdaptor
  Function  : Read the rule objects from the database
  Returntype: array ref of rule objects
  Exceptions: if db isn't the correct type'
  Caller    : 
  Example   : my $rules = &read_db($db);

=cut
sub read_db{
  my $db = shift;

  if(!($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
    throw("need a DBAdaptor not ".$db);
  }
 
  my $rule_adaptor = $db->get_RuleAdaptor;

  my @rules = $rule_adaptor->fetch_all;
  return \@rules;
}


=head2 write_file

  Arg [1]   : filename
      [2]   : arrayref of rule objects
  Function  : write a config file for the objects given
  Returntype: N/A
  Exceptions: if file doesnt exist
  Caller    : 
  Example   : &write_file($file, $rules);

=cut

sub write_file{
  my ($file, $rules) = @_;
  #print STDERR "Opening ".$file."\n";
  open (FH, '>'.$file) or throw ("couldn't open $file to write to");
  foreach my $rule(@$rules){
    
    print FH "[".$rule->goalAnalysis->logic_name."]\n";
    my @conditions = @{$rule->list_conditions};
    foreach my $condition(@conditions){
      print FH "condition=".$condition."\n";
    }
  }
}
