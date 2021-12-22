use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <*.fasta> <min_length>"
}

$len=$ARGV[1];
open $fa,$ARGV[0]||die "$!";
$/='>';
<$fa>;
while(<$fa>){
	chomp;
	@lines=(split /\n+/);
	$head=shift @lines;
	$seq = join '',@lines;
	if (length $seq > $len){
		print '>',$head,"\n",$seq,"\n";
	}

}
