#
# Mickeymouse implementation of RunnableI
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

Bio::EnsEMBL::Pipeline::Runnable::ProcessList

=head1 SYNOPSIS

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::ProcessList();

    my $status   = $runnable->run;
    
    my @output   = $runnable->output;

=head1 DESCRIPTION

Mickeymouse implementation of RunnableI that returns the process list

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::ProcessList;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    return $self;
}


=head2 run
  Title   : run
  Usage   : $self->run
  Function: Runs a ps command and stores the
            output in an array
  Returns : nothing
  Args    : none

=cut

sub run {
    my ($self) = @_;

    open(IN,"ps -ealf |");

    my @proc;
    my $count = 0;

    while ((my $line = <IN>) && ($count < 10)) {
	chomp($line);
	my ($f,$s,$name,$pid,$ppid,$c,$pri,$ni,$addr,$sz,$wchan,$stime,$tty,$time,$cmd) = split(' ',$line,15);

	my $tmp;
	$tmp->{f}     = $f;
	$tmp->{s}     = $s;
	$tmp->{name}  = $name;
	$tmp->{pid}   = $pid;
	$tmp->{ppid}  = $ppid;
	$tmp->{c}     = $c;
	$tmp->{pri}   = $pri;
	$tmp->{ni}    = $ni;
	$tmp->{addr}  = $addr;
	$tmp->{sz}    = $sz;
	$tmp->{wchan} = $wchan;
	$tmp->{stime} = $stime;
	$tmp->{tty}   = $tty;
	$tmp->{time}  = $time;
	$tmp->{cmd}   = $cmd;

	push(@proc,$tmp);
	$count++;
    }

    close(IN);
    $self->{_proc} = [];
    push(@{$self->{_proc}},@proc);

    
}


=head2 output
  Title   : output
  Usage   : my @out = $self->output
  Function: Returns the output from the ps
            command in an array of hashes
            Each element of the array contains
            details of one process
  Returns : @HASH
  Args    : none

=cut

sub output {
    my ($self) = @_;

    return @{$self->{_proc}};

}


1;
