use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <input table> <IDlist>";
}

open $id,$ARGV[1] || die "cannot open :$!";

while(<$id>){
	chomp;
	$ID =(split /\s+/)[0];
	$hash{$ID}  = 1;
}
close $id;

open $f, $ARGV[0] || die "cannot open :$!";

while(<$f>){
	chomp;
	$ID =(split /\s+/)[0];
	print $_,"\n" if exists $hash{$ID} ;

}

close $f;
