package Bio::EnsEMBL::Pipeline::Monitor;

use strict;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Root);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db);
    Function:   creates a pipeline monitoring object if given a pipeline
                database.  Can be used to query and control a running
                pipeline.
    Returns :   
    Args    :   -db:     A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($db) =   $self->_rearrange([qw(DB)],@args);

    if ($db) {
      $self->db($db);
    } else {
      $self->throw("No database object input");
    }

    return $self;
}


sub db {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    if (! $arg->isa("Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor")) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{_db} = $arg;
  }
  return $self->{_db};
}

sub print_header {
  my ($self,$str) = @_;
  
  print "\n$str [" . $self->db->host . " " . $self->db->dbname . "]\n\n";
}

sub show_current_status {
  my ($self) = @_;
  #print STDERR "running current status\n";
  #Show running/failed jobs grouped by status and analysis name.
  my @times = times;
  #print STDERR "1, @times\n";
  my $jobadp = $self->db->get_JobAdaptor;
  
  my $results = $jobadp->list_current_status;
  @times = times;
  #print STDERR "2, @times\n";
  my %status;
  my %tasks;
 JOB:foreach my $result(@$results){
    my ($job_id, $t, $input_id, $s, $timestamp) = @$result;
    #print STDERR "status = ".$s."\n";
    if($s eq 'SUCCESSFUL'){
      next JOB;
    }
    #print STDERR "have task ".$t." status ".$s."\n";
    if(!$tasks{$t}){
      $tasks{$t} = {};
      
    }
    if(!$tasks{$t}{$s}){
      $tasks{$t}{$s} = 0;
     
    }
    $tasks{$t}{$s}++;
   
  }
  @times = times;
  #print STDERR "3, @times\n";
  my $maxcount = 0;
  my $maxstatus = 0;
  my $maxname = 0;
  foreach my $s(keys(%tasks)){
    my %hash = %{$tasks{$s}};
    my @tasks = keys(%hash);
    foreach my $t(@tasks){
      my $count = $status{$s}{$t};
      my $status = $t;
      my $name = $s;
      if(length($count) > $maxcount){
	$maxcount = length($count);
      }
      if (length($status) > $maxstatus) {
	$maxstatus = length($status);
      }
      if (length($name) > $maxname) {
	$maxname = length($name);
      }
    }
  }
  $maxcount++;
  $maxstatus++;
  $maxname++;
  @times = times;
  #print STDERR "4, @times\n";
  $self->print_header("Pipeline current status");
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n","Name","Status","Count");
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n","----","------","-----");

  foreach my $s(keys(%tasks)){
    my %hash = %{$tasks{$s}};
    my @tasks = keys(%hash);
    foreach my $t(@tasks){
      printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n",$s,$t,$tasks{$s}{$t});
    }
  }
  print "\n\n";
  @times = times;
  #print STDERR "5 @times\n\n";
  
  
}




 
sub show_current_status_summary {
  my ($self) = @_;

 my $jobadp = $self->db->get_JobAdaptor;
  
  my $results = $jobadp->list_current_status;
  
  my %status;
  foreach my $result(@$results){
    my ($job_id, $taskname, $input_id, $s, $timestamp) = @$result;
    if(!$status{$s}){
      $status{$s} = 0;
    }
    
    $status{$s}++;
  }
  
  foreach my $s(keys(%status)){
  
  }
  my $maxcount = 0;
  my $maxstatus = 0;
  foreach my $s(keys(%status)){
    my $count = $status{$s};
    my $status = $s;
    if(length($count) > $maxcount){
      $maxcount = length($count);
    }
    if (length($status) > $maxstatus) {
      $maxstatus = length($status);
    }
  }
  $maxcount++;
  $maxstatus++;
  $self->print_header("Pipeline current status summary");
   printf("%-${maxstatus}s %-${maxcount}s\n","Status","Count");
  printf("%-${maxstatus}s %-${maxcount}s\n","------","-----");
  foreach my $s(keys(%status)){
   printf("%-${maxstatus}s %-${maxcount}s\n", $s, $status{$s});
  }
  print "\n\n";
		

}
  

# show finished jobs grouped by analysisId

#

#show running/failed jobs grouped by status

sub show_finished {
  my ($self) = @_;
  #print STDERR "running current status\n";
  #Show running/failed jobs grouped by status and analysis name.
  my @times = times;
  #print STDERR "1, @times\n";
  my $jobadp = $self->db->get_JobAdaptor;
  
  my $results = $jobadp->list_current_status;
  @times = times;
  #print STDERR "2, @times\n";
  my %status;
  my %tasks;
 JOB:foreach my $result(@$results){
    my ($job_id, $t, $input_id, $s, $timestamp) = @$result;
    #print STDERR "status = ".$s."\n";
    if($s ne 'SUCCESSFUL'){
      next JOB;
    }
    #print STDERR "have task ".$t." status ".$s."\n";
    if(!$tasks{$t}){
      $tasks{$t} = {};
      
    }
    if(!$tasks{$t}{$s}){
      $tasks{$t}{$s} = 0;
     
    }
    $tasks{$t}{$s}++;
   
  }
  @times = times;
  #print STDERR "3, @times\n";
  my $maxcount = 0;
  my $maxstatus = 0;
  my $maxname = 0;
  foreach my $s(keys(%tasks)){
    my %hash = %{$tasks{$s}};
    my @tasks = keys(%hash);
    foreach my $t(@tasks){
      my $count = $status{$s}{$t};
      my $status = $t;
      my $name = $s;
      if(length($count) > $maxcount){
	$maxcount = length($count);
      }
      if (length($status) > $maxstatus) {
	$maxstatus = length($status);
      }
      if (length($name) > $maxname) {
	$maxname = length($name);
      }
    }
  }
  $maxcount++;
  $maxstatus++;
  $maxname++;
  @times = times;
  #print STDERR "4, @times\n";
  $self->print_header("Pipeline current status");
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n","Name","Status","Count");
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n","----","------","-----");

  foreach my $s(keys(%tasks)){
    my %hash = %{$tasks{$s}};
    my @tasks = keys(%hash);
    foreach my $t(@tasks){
      printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n",$s,$t,$tasks{$s}{$t});
    }
  }
  print "\n\n";
  @times = times;
  #print STDERR "5 @times\n\n";
  
  
} 
# show analysis processes
sub show_analysis {
  my ($self) = @_;

  my $sth = $self->db->prepare("select analysis_id,logic_name,db,program,parameters,module from analysis");
  my $res = $sth->execute;

  my $maxname;
  my $maxid;
  my $maxdb;
  my $maxprog;
  my $maxparam;
  my $maxmodule;

  my @ids;
  my @names;
  my @dbs;
  my @progs;
  my @params;
  my @modules;

  while (my $ref = $sth->fetchrow_hashref) {
    my $id      = $ref->{'analysis_id'};
    my $name    = $ref->{'logic_name'};
    my $db      = $ref->{'db'};
    my $prog    = $ref->{'program'};
    my $param   = $ref->{'parameters'};
    my $module  = $ref->{'module'};

    if (length($id) > $maxid) {
      $maxid = length($id);
    }
    if (length($name) > $maxname) {
      $maxname = length($name);
    }
    if (length($db) > $maxdb) {
      $maxdb = length($db);
    }
    if (length($prog) > $maxprog) {
      $maxprog = length($prog);
    }
    if (length($module) > $maxmodule) {
      $maxmodule = length($module);
    }
    if (length($param) > $maxparam) {
      $maxparam = length($param);
    }

    push(@names,$name);
    push(@ids,$id);
    push(@dbs,$db);
    push(@progs,$prog);
    push(@modules,$module);
    push(@params,$param);
  }

  $maxname++;
  $maxid++;
  $maxprog++;
  $maxparam++;
  $maxmodule++;
  $maxdb++;
  
  $self->print_header("Analysisprocess");

  printf("%-${maxname}s %-${maxid}s %-${maxdb}s %-${maxprog}s %-${maxparam}s %-${maxmodule}s \n","Name","Id","db","Program","Params","Module");
  printf("%-${maxname}s %-${maxid}s %-${maxdb}s %-${maxprog}s %-${maxparam}s %-${maxmodule}s \n","----","--","--","-------","------","------");

  while (my $name = shift(@names)) {
    my $id   = shift @ids;
    my $db   = shift @dbs;
    my $prog = shift @progs;
    my $param = shift @params;
    my $module = shift @modules;

    printf("%-${maxname}s %-${maxid}s %-${maxdb}s %-${maxprog}s %-${maxparam}s %-${maxmodule}s \n",$name,$id,$db,$prog,$param,$module);
  }

  print "\n";
}




1;

