use strict;
use warnings;


use Bio::EnsEMBL::Pipeline::Job;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::Job);

sub new {
	my $caller = shift;
	my $class = ref($caller) || $caller;

	#call superclass constructor

	#process additional named args
}


sub stdout {


}


sub stderr {

}


sub submission_id {
	
}



1;
