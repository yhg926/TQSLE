use warnings;
use diagnostics;

if(@ARGV != 2){
	die "*.pl <flux.fastq> <K>";
}

$K = $ARGV[1];
open $fq, $ARGV[0] || die "cannot open $ARGV[0]";

while( $head = <$fq> ){

	 $read = <$fq>;  
	 <$fq>;
	 <$fq>;

	$rl = length($read);
	($tid,$tlen,$fst,$fend,$ort) = (split /:/, $head)[2,4,5,6,7];

	if (!exists $hash{$tid}) {
		${$hash{$tid}}[$tlen - $K] = 0; 
	}


	if ($ort =~ /1/){
		$kst = $fst;
		$kend = $fst + $rl - $K;  			
	}
	else {
		$kst = $fend - $rl; 
		$kend = $fend - $K; 
	}

	for( $i = $kst; $i <= $kend; $i++ ){

		${$hash{$tid}}[$i]++;
	}

}
close $fq;


foreach $ele (keys %hash){

	$klen = @{$hash{$ele}};

	for($n = 0 ; $n < $klen; $n++){
		
		${$hash{$ele}}[$n] = 0 if !defined ${$hash{$ele}}[$n] ;
		print $ele,"\t", $klen,"\t", $n, "\t", ${$hash{$ele}}[$n],"\n";

	}

}

