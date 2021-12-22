use warnings;
use diagnostics;
if(@ARGV!=1){
	die "*.pl <*.fasta>";
}

open $fa,$ARGV[0] || die "$!";
$/='>';
<$fa>;
while(<$fa>){
	chomp;
	@lines=(split /\n+/);
	shift @lines;
	
	$seq = join '', @lines;
	
	print length $seq,"\n";
}
