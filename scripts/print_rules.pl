use strict;

use Getopt::Long;

use Bio::EnsEMBL::Pipeline::DBSQL::Obj;

my $host     = 'localhost';
my $port     = '410000';
my $dbname   = 'simon_dec12';
my $dbuser   = 'ensro';
my $pass     =  undef;


&GetOptions(
	    'host:s'     => \$host,
	    'port:n'     => \$port,
	    'dbname:s'   => \$dbname,
	    'dbuser:s'   => \$dbuser,
	    'pass:s'     => \$pass,
	    );



my $db  = new Bio::EnsEMBL::Pipeline::DBSQL::Obj(-host   => $host,
						 -port   => $port,
						 -dbname => $dbname,
						 -user   => $dbuser,
						 -pass   => $pass);

my $sth=$db->prepare("select rg.ruleId,a.program,rc.conditionLiteral from RuleGoal rg, RuleConditions rc, analysis a where rc.ruleId = rg.ruleId and a.id = rg.goalAnalysisId");
$sth->execute();
while (my ($rule,$goal,$condition)= $sth->fetchrow_array) {
    print STDERR "Rule $rule: Do $goal only if $condition\n";
}


