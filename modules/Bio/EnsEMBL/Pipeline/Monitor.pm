package Bio::EnsEMBL::Pipeline::Monitor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use strict;

use vars qw(@ISA);

@ISA = qw ();

=pod

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db);
    Function:   creates a pipeline monitoring object if given a pipeline
                database.  Can be used to query and control a running
                pipeline.
    Returns :   
    Args    :   -dbobj:     A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {},$class;
    my ($dbobj) =   rearrange([qw(DBOBJ)],@args);

    if ($dbobj) {
      $self->dbobj($dbobj);
    } else {
      throw("No database object input");
    }

    return $self;
}


sub dbobj {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    if (! $arg->isa("Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor")) {
      throw("[$arg] is not a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{_dbobj} = $arg;
  }
  return $self->{_dbobj};
}

sub print_header {
  my ($self,$str) = @_;
  unless($self->{'_header'}){
  	$self->{'_header'} = sprintf("[%s:%s %s]\n\n", $self->dbobj->dbc->host, $self->dbobj->dbc->port, 
                                                       $self->dbobj->dbc->dbname);
  }
  print "\n$str " . $self->{'_header'};
}

sub show_current_status {
  my ($self) = @_;

  #Show running/failed jobs grouped by status and analysis name.

  my $sth = $self->dbobj->prepare("select count(*), js.status, a.logic_name ,a.analysis_id from analysis a, job_status js, job j where j.job_id = js.job_id and a.analysis_id = j.analysis_id and js.is_current = 'y' group by a.logic_name, js.status");

  my $res = $sth->execute;

  my $maxcount = undef;
  my $maxstatus;
  my $maxname;
  my $maxaid;
  my @counts;
  my @status;
  my @names;
  my @aid;

  while (my $ref = $sth->fetchrow_hashref) {
    my $count  = $ref->{'count(*)'};
    my $status = $ref->{'status'};
    my $name   = $ref->{'logic_name'};
    my $aid    = $ref->{'analysis_id'};

    if (!defined($maxcount) || length($count) > $maxcount) {
      $maxcount = length($count);
    }
    if (!defined($maxstatus) || length($status) > $maxstatus) {
      $maxstatus = length($status);
    }
    if (!defined($maxname) || length($name) > $maxname) {
      $maxname = length($name);
    }
    if (!defined($maxaid) || length($aid) > $maxaid) {
      $maxaid = length($aid);
    }

    push(@counts,$count);
    push(@status,$status);
    push(@names,$name);
    push(@aid,$aid);

  }
    $maxcount++;
    $maxstatus++;
    $maxname++;
    $maxaid++;

  $self->print_header("Pipeline current status");

  
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s   %-${maxaid}s\n","Name","Status","Count","Analysis-id");
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s   %-${maxaid}s\n","----","------","-----","-----------");
  
  while (my $count = shift(@counts)) {
    my $status = shift @status;
    my $name   = shift @names;
    my $aid    = shift @aid; 

    printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s     %-${maxaid}s \n",$name,$status,$count,$aid);
  }

  print("\n");

}


  # Show running/failed jobs grouped by status and analysisId
#  my $sth = $self->dbobj->prepare("select count(*),js.status,job.analysis_id from current_status js,job where job.job_id = js.job_id group by job.analysis_id";


# show running/failed jobs grouped by status
sub show_current_status_summary {
  my ($self) = @_;

  my $sth = $self->dbobj->prepare("select count(*), status from job_status where is_current = 'y' group by status");

  my $res = $sth->execute;

  my $maxcount = 0;
  my $maxstatus = 0;

  my @counts;
  my @status;

  while (my $ref = $sth->fetchrow_hashref) {
    my $count  = $ref->{'count(*)'};
    my $status = $ref->{'status'};

    if (length($count) > $maxcount) {
      $maxcount = length($count);
    }
    if (length($status) > $maxstatus) {
      $maxstatus = length($status);
    }

    push(@counts,$count);
    push(@status,$status);

  }
  $maxcount++;
  $maxstatus++;
  
  $self->print_header("Pipeline status summary");

  printf("%-${maxstatus}s %-${maxcount}s\n","Status","Count");
  printf("%-${maxstatus}s %-${maxcount}s\n","------","-----");

  while (my $count = shift(@counts)) {
    my $status = shift @status;

    printf("%-${maxstatus}s %-${maxcount}s\n",$status,$count);
  }

  print "\n";

}
  

# show finished jobs grouped by analysisId

#

#show running/failed jobs grouped by status
sub show_finished_summary {
  my ($self, $no_submit, $show_percent) = @_;

  my $sth = $self->dbobj->prepare("select count(*),a.logic_name,a.analysis_id from input_id_analysis i, analysis  a where a.analysis_id = i.analysis_id group by a.analysis_id");

  my $res = $sth->execute;

  my $maxcount = 0;
  my $maxname  = 0;
  my $maxid    = 0;
  my $max_nums = undef;
  my @counts;
  my @names;
  my @ids;

  while (my $ref = $sth->fetchrow_hashref) {
    my $count  = $ref->{'count(*)'};
    my $name = $ref->{'logic_name'};
    my $id = $ref->{'analysis_id'};

    if ($name =~ /Submit/){
      next if($no_submit);
    }

    if (length($count) > $maxcount) {
      $maxcount = length($count);
    }

    if (length($name) > $maxname) {
      $maxname = length($name);
    }

    if (length($id) > $maxid) {
      $maxid = length($id);
    }

    push(@counts,$count);
    push(@names,$name);
    push(@ids,$id);
  }

  $maxcount++;
  $maxname++;
  $maxid++;

  if($show_percent){
    $max_nums = percent_finished($self, 0);
    $self->print_header("Finished job summary");
    printf("%-${maxcount}s %-${maxname}s %-${maxid}s%s\n","Count","Name","Id","%Done");
    printf("%-${maxcount}s %-${maxname}s %-${maxid}s%s\n","-----","----","--","-----");

    while (my $count = shift(@counts)) {
      my $name  = shift @names;
      my $id    = shift @ids;
      my $total = 0;
      my $add   = " ";
      if($count && $max_nums->{$id}){
	$total = (100*$count/($max_nums->{$id}));
	#$total = $count."/".$max_nums->{$id};
	$add = '%';
      }
      printf("%-${maxcount}s %-${maxname}s %-${maxid}s %3d%s\n",$count,$name,$id,$total,$add);
    }
  }
  else{
    $self->print_header("Finished job summary");

    printf("%-${maxcount}s %-${maxname}s %-${maxid}s\n","Count","Name","Id");
    printf("%-${maxcount}s %-${maxname}s %-${maxid}s\n","-----","----","--");

    while (my $count = shift(@counts)) {
      my $name = shift @names;
      my $id   = shift @ids;

      printf("%-${maxcount}s %-${maxname}s %-${maxid}s\n",$count,$name,$id);
    }
  }
  print "\n";

}

# show analysis processes
sub show_analysisprocess {
  my ($self) = @_;

  my $sth = $self->dbobj->prepare("select analysis_id,logic_name,db,program,parameters,module from analysis");
  my $res = $sth->execute;

  my $maxname = 0;
  my $maxid = 0;
  my $maxdb = 0;
  my $maxprog = 0;
  my $maxparam = 0;
  my $maxmodule = 0;

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
  
  $self->print_header("Analysis");

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



# show rules
sub show_Rules {
  my ($self) = @_;

  my $sth = $self->dbobj->prepare("select a.logic_name,rg.rule_id from rule_goal rg, analysis a where a.analysis_id = rg.goal");

  my $res = $sth->execute;

  my @names;
  my @ids;
 
  my $maxname = 0;
  my $maxid = 0;

  while (my $ref = $sth->fetchrow_hashref) {
    my $id = $ref->{'rule_id'};
    my $name = $ref->{'logic_name'};

    if (length($id) > $maxid) { $maxid = length($id);}
    if (length($name) > $maxname) {$maxname = length($name);}

    push(@ids,$id);
    push(@names,$name);
  }
  
  $maxname++;
  $maxid++;

  $self->print_header("Pipeline rules");

  printf("%-${maxname}s %-${maxid}s\n","Name","Id");
  printf("%-${maxname}s %-${maxid}s\n","----","--");

  while (my $name = shift @names) {
    my $id = shift @ids;
    printf("%-${maxname}s %-${maxid}s\n",$name,$id);
  }

  print "\n";
}

sub show_Rules_and_Conditions {
  my ($self) = @_;

  my $sql = "select a.logic_name,rg.rule_id,rc.rule_condition from rule_conditions rc,rule_goal rg, analysis a where a.analysis_id = rg.goal and rg.rule_id = rc.rule_id";
  
  my $sth = $self->dbobj->prepare($sql);

  my $res = $sth->execute;

  my @names;
  my @ids;
  my @conds;

  my $maxname = 0;
  my $maxid = 0;
  my $maxcond = 0;

  while (my $ref = $sth->fetchrow_hashref) {
    my $id = $ref->{'rule_id'};
    my $name = $ref->{'logic_name'};
    my $cond = $ref->{'rule_condition'};

    if (length($id) > $maxid) { $maxid = length($id);}
    if (length($name) > $maxname) {$maxname = length($name);}
    if (length($cond) > $maxcond) {$maxcond = length($cond);}

    push(@ids,$id);
    push(@names,$name);
    push(@conds,$cond);
  }
  
  $maxname++;
  $maxid++;
  $maxcond++;

  $self->print_header("Rules and Conditions");

  printf("%-${maxname}s %-${maxid}s %-${maxcond}s \n","Name","Id","Condition");
  printf("%-${maxname}s %-${maxid}s %-${maxcond}s \n","----","--","---------");
  

  while (my $name = shift @names) {
    my $id = shift @ids;
    my $cond = shift @conds;

    printf("%-${maxname}s %-${maxid}s %-${maxcond}s \n",$name,$id,$cond);
  }

  print "\n";
}

sub show_jobs_by_status_and_analysis {
  my ($self,$status,$analysis) = @_;

  if (!defined($status) || !defined($analysis)) {
    throw("No status and/or analysis input\n");
  }

#  my $sth = $self->dbobj->prepare("select job.* from job_status js,job,analysis a where a.analysis_id = job.analysis_id and job.job_id = js.job_id and js.status = '$status' and a.logic_name = '$analysis'");
  my $sth = $self->dbobj->prepare("
  	SELECT job.* FROM job_status js,job,analysis a 
  	WHERE a.analysis_id = job.analysis_id 
  		&& job.job_id = js.job_id 
  		&& js.status = '$status' 
  		&& a.logic_name = '$analysis'
  		&& js.is_current = 'y'");
  my $res = $sth->execute;

  my @jobIds;
  my @input_ids;
  my @lsfs;
  my @out;
  my @err;
  my @retry;
  my @exhost;

  my $maxjobid = 0;
  my $maxinputid = 0;
  my $maxlsf = 0;
  my $maxout = 0;
  my $maxerr = 0;
  my $maxretry = 0;
  my $maxexhost = 0;

  while (my $ref = $sth->fetchrow_hashref) {
    my $jobId         = $ref->{'job_id'};
    my $input_id      = $ref->{'input_id'};
    my $submission_id = $ref->{'submission_id'};
    my $out           = $ref->{'stdout_file'};
    my $err           = $ref->{'stderr_file'};
    my $retry         = $ref->{'retry_count'};
    my $exhost        = $ref->{'exec_host'};

    if (length($jobId) > $maxjobid) {$maxjobid = length($jobId);}
    if (length($input_id) > $maxinputid) {$maxinputid = length($input_id);}
    if (length($submission_id) > $maxlsf) {$maxlsf = length($submission_id);}
    if (length($out) > $maxout) {$maxout = length($out);}
    if (length($err) > $maxerr) {$maxerr = length($err);}
    if (length($retry) > $maxretry) {$maxretry = length($retry);}
    if (length($exhost) > $maxexhost) {$maxexhost = length($exhost);}

    push(@jobIds,$jobId);
    push(@input_ids,$input_id);
    push(@lsfs,$submission_id);
    push(@out,$out);
    push(@err,$err);
    push(@retry,$retry);
    push(@exhost,$exhost);

  }

  $self->print_header("Jobs by status $status and analysis $analysis");
  
  printf("%-${maxinputid}s  %-${maxjobid}s  %-${maxlsf}s  %-${maxretry}s  %-${maxout}s  %-${maxerr}s  %-${maxexhost}s\n","Input id","Job","submission id","retry","output file","error file","exhost");
  printf("%-${maxinputid}s  %-${maxjobid}s  %-${maxlsf}s  %-${maxretry}s  %-${maxout}s  %-${maxerr}s  %-${maxexhost}s\n","--------","---","-------------","-----","-----------","----------","------");

  while (my $jobid = shift @jobIds) {
    my $inputid = shift @input_ids;
    my $lsf     = shift @lsfs;
    my $out     = shift @out;
    my $err     = shift @err;
    my $retry   = shift @retry;
    my $exhost  = shift @exhost;
    
    printf("%-${maxinputid}s  %-${maxjobid}s  %-${maxlsf}s  %-${maxretry}s  %-${maxout}s  %-${maxerr}s  %-${maxexhost}s\n",$inputid,$jobid,$lsf,$retry,$out,$err,$exhost);
  }

  print "\n";
}


# ==================================================== #

sub rules_cache{
    my $self = shift;
    unless($self->{'_rules_cache'}){
            my $rule_adaptor = $self->dbobj->get_RuleAdaptor();
            my @rules = $rule_adaptor->fetch_all;
            foreach my $rule (@rules)  {
                $self->{'_rules_cache'}->{$rule->goalAnalysis->dbID} = 
		  $rule->goalAnalysis->logic_name if ($rule->list_conditions()->[0]);
            }
    }
    return $self->{'_rules_cache'};
}

=pod

=head2 get_unfinished_analyses***

The following methods return an anonymous list of arrays:
 get_unfinished_analyses_for_input_id( $input_id )
 get_unfinished_analyses_for_assembly_type( $assembly_type )
 get_unfinished_analyses

with the following data spec:
 [
   [ input_id, logic_name, analysis_id ]
   [ input_id, logic_name, analysis_id ]
   ...
 ]
 
They all take an optional $print (Boolean) which prints the arrays for you.
 
You can grow you're own look up hash to see if a $input_id, $logic_name combination has
finished by doing.
 
 my $unfin = $monitor_obj->get_unfinished_analyses_for_assembly_type($assembly_type);
 my $hash;
 map { $hash->{$_->[0]}->{$_->[1]} = $_->[2] } @$unfin;
 
Now B<$hash-E<gt>{$input_id}-E<gt>{$logic_name} = $analysis_id>.

=head2 get_no_hit_contigs_for_analysis

 Takes two arguments, and finds all the finished input_ids for which the given analysis
did NOT find any hits.
 Returns data structure as above.
 
 get_no_hit_contigs_for_analysis($feature_table, $analysis_id)
'
=cut

sub get_unfinished_analyses_for_input_id{
    my ($self, $contig_name, $print, $unfinished) = @_;
    my $rules_cache = $self->rules_cache();
    my $sth = $self->dbobj->prepare(qq{SELECT c.name, a.analysis_id AS a_id, a.logic_name 
                                                    FROM contig c STRAIGHT_JOIN analysis a 
                                                    LEFT JOIN input_id_analysis i ON c.name = i.input_id 
				                    && a.analysis_id = i.analysis_id 
                                                    WHERE i.input_id IS NULL && c.name = ?});
    $sth->execute($contig_name);
    while(my $row  = $sth->fetchrow_hashref){
            push(@{$unfinished}, [ $row->{'name'}, $row->{'logic_name'}, $row->{'a_id'} ])
                            if $rules_cache->{$row->{'a_id'}};
    }
    map {print join("\t: ", @$_) . "\n"} @$unfinished if $print;
    return $unfinished;
}

sub get_unfinished_analyses_for_assembly_type{
    my ($self, $assembly_type, $print, $unfinished) = @_;
    my $rules_cache = $self->rules_cache();
    my $sth = $self->dbobj->prepare(qq{SELECT c.name, a.analysis_id AS a_id, a.logic_name, b.type 
                                         FROM assembly b, contig c STRAIGHT_JOIN analysis a 
                                    LEFT JOIN input_id_analysis i ON c.name = i.input_id
					   && a.analysis_id = i.analysis_id
                                        WHERE i.input_id IS NULL 
					   && b.type = ?
					   && b.contig_id = c.contig_id});
    $sth->execute($assembly_type);
    while(my $row  = $sth->fetchrow_hashref){
            push(@{$unfinished}, [ $row->{'name'}, $row->{'logic_name'}, $row->{'a_id'} ])
                            if $rules_cache->{$row->{'a_id'}};
    }
    map {print join("\t: ", @$_) . "\n"} @$unfinished if $print;
    return $unfinished || [];
}

sub get_unfinished_analyses{
    my ($self, $print, $unfinished) = @_;
    my $rules_cache = $self->rules_cache();
    my $sth = $self->dbobj->prepare(qq{SELECT c.name, a.analysis_id AS a_id, a.logic_name 
                                                    FROM contig c STRAIGHT_JOIN analysis a 
                                                    LEFT JOIN input_id_analysis i ON c.name = i.input_id && a.analysis_id = i.analysis_id 
                                                    WHERE i.input_id IS NULL ORDER BY c.name});
    $sth->execute();
    while(my $row  = $sth->fetchrow_hashref){
            push(@$unfinished, [ $row->{'name'}, $row->{'logic_name'}, $row->{'a_id'} ])
                            if $rules_cache->{$row->{'a_id'}};
    }
    map {print join("\t: ", @$_) . "\n"} @$unfinished if $print;
    print "Waiting for " . scalar(@$unfinished) . " to complete.\n" if $print;
    return $unfinished;
}

sub get_no_hit_contigs_for_analysis{
    my $self = shift;
    my ($feature_table, $analysis_id, $print, $no_hits) = @_;
    my $sth = $self->dbobj->prepare(qq{SELECT c.contig_id, c.name
                                                    FROM contig c STRAIGHT_JOIN input_id_analysis i 
                                                    LEFT JOIN $feature_table f ON f.analysis_id = i.analysis_id && f.contig_id = c.contig_id  
                                                    WHERE f.contig_id IS NULL && f.analysis_id IS NULL && i.analysis_id = ? && c.name = i.input_id});
    $sth->execute($analysis_id);
    $no_hits = [];
    while(my $row = $sth->fetchrow_hashref){
            push(@{$no_hits}, [ $row->{'name'}, undef, $analysis_id ] );
    }
    map {print join("\t: ", @$_) . "\n"} @$no_hits if $print;
    return $no_hits;
}

sub lock_status{
    my ($self,$print) = @_;
    my $str = "Locked By USER: %s, HOST: %s, PID: %s, STARTED: %s (%s) \n";
    my @a = ();
    my ($dbuser, $dbhost, $dbport, $dbname) = ($self->dbobj->dbc->username,
                                               $self->dbobj->dbc->host,
                                               $self->dbobj->dbc->port,
                                               $self->dbobj->dbc->dbname);
    if (my $lock_str = $self->dbobj->pipeline_lock) {
            my($user, $host, $pid, $started) = $lock_str =~ /(\w+)@(\S+):(\d+):(\d+)/;
            $self->print_header("This pipeline is LOCKED") if $print;
            @a = ($user, $host, $pid, scalar(localtime($started)), $started);
            printf($str, @a) if $print;
    }else{
            $self->print_header("This pipeline is FREE") if $print;
    }
    return @a;
}


#get the maximum numbers of jobs for every analysis
#for calculating the percent-done numbers
sub percent_finished{
  my ($self, $print) = @_;

  my $sql0 = "select a.analysis_id, a.logic_name ".
             "from analysis a;";
  my $sth0 = $self->dbobj->prepare($sql0) or die "cant prepare: $sql0";
  my $sql1 = "select rg.goal, rc.rule_condition ".
             "from rule_conditions rc,rule_goal rg ".
	     "where rg.rule_id = rc.rule_id;";
  my $sth1 = $self->dbobj->prepare($sql1) or die "cant prepare: $sql1";
  my $sql2 = "select count(*) as count from input_id_analysis ".
             "where analysis_id = ?;";
  my $sth2 = $self->dbobj->prepare($sql2) or die "cant prepare: $sql2";

  my %analysis_id        = ();
  my %analysis_name      = ();
  my %analysis_maxnum    = ();
  my %analysis_condition = ();
  my $ref3;

  $sth0->execute();
  while (my $ref0 = $sth0->fetchrow_hashref) {
    $analysis_name{$ref0->{'analysis_id'}}      = $ref0->{'logic_name'};
    $analysis_id{$ref0->{'logic_name'}}         = $ref0->{'analysis_id'};
    if($ref0->{'logic_name'} =~ /Wait/i){
      $analysis_maxnum{$ref0->{'analysis_id'}}  = 1;
    }
    else{
      $analysis_maxnum{$ref0->{'analysis_id'}}  = 0;
    }
  }
  $sth1->execute();
  while (my $ref1 = $sth1->fetchrow_hashref) {
    next if($ref1->{'rule_condition'} =~ m/.+wait/i);
    $analysis_condition{$ref1->{'goal'}} = $ref1->{'rule_condition'};
  }

  foreach my $analysis_id (keys %analysis_name){
    next if($analysis_name{$analysis_id} =~ /Wait/i);
    my $initial_id = $analysis_id;
    $sth2->execute($initial_id);
    $ref3 = $sth2->fetchrow_hashref;
    $analysis_maxnum{$initial_id} = $ref3->{'count'};

    while(defined($analysis_condition{$analysis_id})){
      $analysis_id = $analysis_id{ $analysis_condition{$analysis_id} };
      $sth2->execute($analysis_id);
      $ref3 = $sth2->fetchrow_hashref;
      if($ref3->{'count'} > $analysis_maxnum{$initial_id}){
	$analysis_maxnum{$initial_id} = $ref3->{'count'};
      }
      #print "\n\t-> ".$analysis_id.": ".$analysis_maxnum{$initial_id};
    }
  }

  #print out for testing
  if($print){
    print "\n\n\nanalysis_id\tlogic_name\tmax_num\n";
    foreach my $analysis_id (sort{$a<=>$b} keys %analysis_name){
      print $analysis_name{$analysis_id}."\t\t".
	    $analysis_id."\t".
	    $analysis_maxnum{$analysis_id}."\n";
    }
  }

  return(\%analysis_maxnum);
}


1;
