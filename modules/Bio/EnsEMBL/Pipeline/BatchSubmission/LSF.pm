package Bio::EnsEMBL::Pipeline::BatchSubmission::LSF;


use Bio::EnsEMBL::Pipeline::BatchSubmission;
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Pipeline::BatchSubmission);


sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'bsub'} = undef;
  
  return $self;
 
}


##################
#accessor methods#
##################

sub bsub{
  my($self, $arg) = @_;

  if(defined($arg)){
    $self->{'bsub'} = $arg;
  }

  return $self->{'bsub'};

}

##other accessor are in base class##

######################
#command line methods#
######################

sub construct_command_line{
  my($self, $command, $stdout, $stderr) = @_; 
#command must be the first argument then if stdout or stderr aren't definhed the objects own can be used
  if(!$command){
    $self->throw("cannot create bsub if nothing to submit to it : $!\n");
  }
  my $bsub_line;
  $self->command($command);
  if($stdout){
    $bsub_line = "bsub -o ".$stdout;
  }else{
    $bsub_line = "bsub -o ".$self->stdout_file;
  }
  if($self->nodes){
    my $nodes = $self->nodes;
    # $nodes needs to be a space-delimited list
    $nodes =~ s/,/ /;
    $nodes =~ s/ +/ /;
    # undef $nodes unless $nodes =~ m{(\w+\ )*\w};
    $bsub_line .= " -m '".$nodes."' ";
  }
  $bsub_line .= " -q ".$self->queue    if defined $self->queue;
  $bsub_line .= " -J ".$self->jobname  if defined $self->jobname;
  $bsub_line .= " ".$self->parameters." "  if defined $self->parameters;
  if($stderr){
    $bsub_line .= " -e ".$stderr;
  }else{
    $bsub_line .= " -e ".$self->stderr_file;
  }
  $bsub_line .= " -E \"".$self->pre_exec."\"" if defined $self->pre_exec; 
  ## must ensure the prexec is in quotes ##
  $bsub_line .= " ".$command;
  $self->bsub($bsub_line);
  #print "have command line\n";
}



sub open_command_line{
  my ($self)= @_;

  open(SUB, $self->bsub." 2>&1 |");
  my $lsf;
  while(<SUB>){
    if (/Job <(\d+)>/) {
      $lsf = $1;
    }
  }
  $self->id($lsf);
  close(SUB);
}


sub get_pending_jobs {
  my($self, %args) = @_;

  my ($user)  = $args{'-user'}  || $args{'-USER'}  || undef;
  my ($queue) = $args{'-queue'} || $args{'-QUEUE'} || undef;

  my $cmd = "bjobs -p";
  $cmd   .= " -q $queue" if $queue;
  $cmd   .= " -u $user"  if $user;
  $cmd   .= " | grep PEND";

  open CMD, "$cmd 2>&1 |" or do {
    return undef;
  };

  my @lines = <CMD>;

  close CMD or do {
    return undef;
  };

  return 0 if $lines[0] =~ /No pending job found/;
  return scalar @lines;
}

1;
