package Bio::EnsEMBL::Pipeline::Monitor;

use strict;

use vars qw(@ISA);

@ISA = qw (Bio::Root::RootI);

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
    my $self = $class->SUPER::new(@args);

    my ($dbobj) =   $self->_rearrange([qw(DBOBJ)],@args);

    if ($dbobj) {
      $self->dbobj($dbobj);
    } else {
      $self->throw("No database object input");
    }

    return $self;
}


sub dbobj {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    if (! $arg->isa("Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor")) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{_dbobj} = $arg;
  }
  return $self->{_dbobj};
}

sub print_header {
  my ($self,$str) = @_;

  print "\n$str [" . $self->dbobj->host . " " . $self->dbobj->dbname . "]\n\n";
}

sub show_current_status {
  my ($self) = @_;

  #Show running/failed jobs grouped by status and analysis name.

  my $sth = $self->dbobj->prepare("select count(*),js.status,ap.logic_name from analysisprocess ap,current_status js,job where job.jobId = js.jobId and ap.analysisId = job.analysisId group by ap.logic_name,js.status");

  my $res = $sth->execute;

  my $maxcount;
  my $maxstatus;
  my $maxname;

  my @counts;
  my @status;
  my @names;

  while (my $ref = $sth->fetchrow_hashref) {
    my $count  = $ref->{'count(*)'};
    my $status = $ref->{'status'};
    my $name   = $ref->{'logic_name'};

    if (length($count) > $maxcount) {
      $maxcount = length($count);
    }
    if (length($status) > $maxstatus) {
      $maxstatus = length($status);
    }
    if (length($name) > $maxname) {
      $maxname = length($name);
    }

    push(@counts,$count);
    push(@status,$status);
    push(@names,$name);

  }
    $maxcount++;
    $maxstatus++;
    $maxname++;

  $self->print_header("Pipeline current status");

  
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n","Name","Status","Count");
  printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n","----","------","-----");

  while (my $count = shift(@counts)) {
    my $status = shift @status;
    my $name   = shift @names;


    printf("%-${maxname}s %-${maxstatus}s %-${maxcount}s\n",$name,$status,$count);
  }

  print("\n");

}


  # Show running/failed jobs grouped by status and analysisId
#  my $sth = $self->dbobj->prepare("select count(*),js.status,job.analysisId from current_status js,job where job.jobId = js.jobId group by job.analysisId";


# show running/failed jobs grouped by status
sub show_current_status_summary {
  my ($self) = @_;

  my $sth = $self->dbobj->prepare("select count(*),status from current_status group by status");

  my $res = $sth->execute;

  my $maxcount;
  my $maxstatus;

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
  my ($self) = @_;

  my $sth = $self->dbobj->prepare("select count(*),ap.logic_name,ap.analysisId from InputIdAnalysis i, analysisprocess ap where ap.analysisId = i.analysisId group by ap.analysisId;");
  
  my $res = $sth->execute;

  my $maxcount;
  my $maxname;
  my $maxid;

  my @counts;
  my @names;
  my @ids;

  while (my $ref = $sth->fetchrow_hashref) {
    my $count  = $ref->{'count(*)'};
    my $name = $ref->{'logic_name'};
    my $id = $ref->{'analysisId'};

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

  $self->print_header("Finished job summary");

  printf("%-${maxcount}s %-${maxname}s %-${maxid}s\n","Count","Name","Id");
  printf("%-${maxcount}s %-${maxname}s %-${maxid}s\n","-----","----","--");

  while (my $count = shift(@counts)) {
    my $name = shift @names;
    my $id   = shift @ids;

    printf("%-${maxcount}s %-${maxname}s %-${maxid}s\n",$count,$name,$id);
  }
  print "\n";

}
  
# show analysis processes
sub show_analysisprocess {
  my ($self) = @_;

  my $sth = $self->dbobj->prepare("select analysisId,logic_name,db,program,parameters,module from analysisprocess");
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
    my $id      = $ref->{'analysisId'};
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



# show rules
sub show_Rules {
  my ($self) = @_;

  my $sth = $self->dbobj->prepare("select ap.logic_name,rg.ruleid from RuleGoal rg, analysisprocess ap where ap.analysisId = rg.goalAnalysisId");

  my $res = $sth->execute;

  my @names;
  my @ids;
 
  my $maxname;
  my $maxid;

  while (my $ref = $sth->fetchrow_hashref) {
    my $id = $ref->{'ruleid'};
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

  my $sth = $self->dbobj->prepare("select ap.logic_name,rg.ruleid,rc.conditionLiteral from RuleConditions rc,RuleGoal rg, analysisprocess ap where ap.analysisId = rg.goalAnalysisId and rg.ruleId = rc.ruleId");

  my $res = $sth->execute;

  my @names;
  my @ids;
  my @conds;

  my $maxname;
  my $maxid;
  my $maxcond;

  while (my $ref = $sth->fetchrow_hashref) {
    my $id = $ref->{'ruleid'};
    my $name = $ref->{'logic_name'};
    my $cond = $ref->{'conditionLiteral'};

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
    $self->throw("No status and/or analysis input\n");
  }
	       

  my $sth = $self->dbobj->prepare("select job.* from jobstatus js,job,analysisprocess ap where ap.analysisId = job.analysisId and job.jobId = js.jobId and js.status = '$status' and ap.logic_name = '$analysis'");

  my $res = $sth->execute;

  my @jobIds;
  my @input_ids;
  my @lsfs;
  my @out;
  my @err;
  my @retry;

  my $maxjobid;
  my $maxinputid;
  my $maxlsf;
  my $maxout;
  my $maxerr;
  my $maxretry;

  while (my $ref = $sth->fetchrow_hashref) {
    my $jobId = $ref->{'jobId'};
    my $input_id = $ref->{'input_id'};
    my $LSF_id = $ref->{'LSF_id'};
    my $out    = $ref->{'stdout_file'};
    my $err    = $ref->{'stderr_file'};
    my $retry  = $ref->{'retry_count'};

    if (length($jobId) > $maxjobid) {$maxjobid = length($jobId);}
    if (length($input_id) > $maxinputid) {$maxinputid = length($input_id);}
    if (length($LSF_id) > $maxlsf) {$maxlsf = length($LSF_id);}
    if (length($out) > $maxout) {$maxout = length($out);}
    if (length($err) > $maxerr) {$maxerr = length($err);}
    if (length($retry) > $maxretry) {$maxretry = length($retry);}

    push(@jobIds,$jobId);
    push(@input_ids,$input_id);
    push(@lsfs,$LSF_id);
    push(@out,$out);
    push(@err,$err);
    push(@retry,$retry);

  }

  $self->print_header("Jobs by status $status and analysis $analysis");
  
  printf("%-${maxinputid}s  %-${maxjobid}s  %-${maxlsf}s  %-${maxretry}s  %-${maxout}s  %-${maxerr}s\n","Input id","Job","LSF id","retry","output file","error file");
  printf("%-${maxinputid}s  %-${maxjobid}s  %-${maxlsf}s  %-${maxretry}s  %-${maxout}s  %-${maxerr}s\n","--------","---","------","-----","-----------","----------");

  while (my $jobid = shift @jobIds) {
    my $inputid = shift @input_ids;
    my $lsf     = shift @lsfs;
    my $out     = shift @out;
    my $err     = shift @err;
    my $retry   = shift @retry;
    
    printf("%-${maxinputid}s  %-${maxjobid}s  %-${maxlsf}s  %-${maxretry}s  %-${maxout}s  %-${maxerr}s\n",$inputid,$jobid,$lsf,$retry,$out,$err);
  }

  print "\n";
}
1;

