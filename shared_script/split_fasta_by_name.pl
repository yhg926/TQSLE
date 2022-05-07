use warnings;
use diagnostics;
if (@ARGV!=1){
	die "*.pl <.fasta>";
}

$/='>';
open $fas,$ARGV[0] || die "can not open $ARGV[0]:$!";
<$fas>;

while(<$fas>){
	chomp;
	$head=(split /\s+/)[0];
	open ($out,'>', $head.".fasta") or die $!;
	print $out '>'.$_ ;
}

