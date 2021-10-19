#input format: Transcript_ID\tGene_ID\tAbundace\n
#!/usr/bin/perl
use warnings;
use diagnostics;
if(@ARGV!=1){
	die "*.pl <Transcript_ID2Gene_ID4Abundace>";
}

open $f,$ARGV[0] || die "cannot open $ARGV[0]" ;

while(<$f>){
	chomp;
	($gid,$val)=(split /\s+/)[1,2] ;
	$hash{$gid} += $val ;
}

foreach $ele (keys %hash){
	print $ele,"\t",$hash{$ele},"\n";
	
}
