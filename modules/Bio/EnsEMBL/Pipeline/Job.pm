use strict;
use warnings;




sub new {

}


sub dbID {
  my $self = shift;

	if(@_) {
		$self->{'dbID'} = shift;
	}

	return $self->{'dbID'};
}


sub task_name {
	my $self = shift;

	if(@_) {
		$self->{'task_name'} = shift;
	}

	return $self->{'task_name'};
}


sub adaptor {
	my $self = shift;


}


sub input_id {

}




=head2 module

  Arg [1]    : string $module
               A fully qualified module name
  Example    : $job->module('Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker');
  Description: Getter setter for the fully qualified module name that does the
               work associated with this job.  The module is required to
               implement the methods fetch_input, run, write_output
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub module {
	my $self = shift;

	if(@_) {
		$self->{'module'} = shift;
	}

	return $self->{'module'};
}



=head2 run

  Arg [1]    : none
  Example    : $job->run();
  Description: Runs this job locally. This will require the module
               specified by the module() subroutine and instantiate it with the
               parameters from the parameters() subroutine.  Upon instantiation
               the fetch_input(), run(), and write_output() subroutines will
               be called, and the job status table updated accordingly.
  Returntype : none
  Exceptions : none
  Caller     : pipeLineManager, runner.pl script (on farm)

=cut

sub run {
  my $self = shift;

	my $module = $self->module();

  eval {

  #require the module

  # instantiate with parameter constructor args

  #set job status to reading

  #call $module->fetch_input()

  #set job status to running

  #call $module->run()

  #set job status to writing

  #call $module->write_output();

  };

  if($@) {
   #set job status to failed
  }

  #set job status to successful

  return;
}





=head2 parameters

  Arg [1]    : 
  Example    : 
  Description: Getter/Setter for Additional parameters apart from the input id
               that have to
               be passed into the running module.  This is a free form
               string retrieved from the job table.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub parameters {
	my $self = shift;
}


=head2 get_status

  Arg [1]    : none
  Example    :
    foreach $hashref (@{$job->get_status}) { 
        print $hashref->{'status'}, ' ', $hashref->{'time'}; 
    }
  Description: Retrieves listref of hashrefs containing two keys 'status' and
               'time'.  The first value of the listref is the most recent
               status, the last value is the first status (CREATED entry).
  Returntype : listref of hashrefs
  Exceptions : none
  Caller     : general

=cut

sub get_status {
	my $self = shift;

  return $self->{'status'};
}


#
# sets the current status, and time, must still be stored in DB
#
sub set_current_status {
  my $self = shift;
  my $status = shift;

  my $time = time();
}





1;

