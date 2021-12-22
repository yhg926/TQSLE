#!/usr/bin/perl
use warnings;
use diagnostics;
if(@ARGV!=2){
    die "*.pl <Transcript abundance> <allTranscript2GeneID>";
}
open $t2g,$ARGV[1] || die "cannot open $ARGV[1]" ;

while(<$t2g>){
	chomp;
	($tid, $gid)=(split /\s+/)[0,1];
	$hash{$tid} = $gid;

}
close $t2g;

#input format: Transcript_ID\tAbundace\n
open $f,$ARGV[0] || die "cannot open $ARGV[0]" ;

while(<$f>){
	chomp;
	($tid,$x)=(split /\s+/)[0,1];

	$gid = $hash{$tid};
	if(!exists $sum{$gid}) { $sum{$gid} = 0 };
	$sum{$gid} += $x;
}
close($f);

foreach $gid (sort keys %sum){

	print $gid,"\t",$sum{$gid},"\n";

}



