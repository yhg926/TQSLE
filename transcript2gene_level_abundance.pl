use warnings;
use diagnostics;
if(@ARGV!=2){
	die "Usage: *.pl <abundance.out> <transcriptID2geneID.tab>";
}

open $t2g,$ARGV[1] || die "cannot open $ARGV[1]:$!";

while(<$t2g>){
	chomp;
	($tid,$gid) = (split /\s+/)[0,1];
	$hash{$tid} = $gid;
}
close $t2g;

open $ab,$ARGV[0] || die "cannot open $ARGV[0]:$!";

while(<$ab>){
	chomp;
	($tid,$tpm)=(split /\s+/)[0,-1];
	next if !exists $hash{$tid};

	$gene{$hash{$tid}} += $tpm;

}
close $ab;

foreach $gid (sort keys %gene){
	print $gid,"\t",$gene{$gid},"\n";

}

