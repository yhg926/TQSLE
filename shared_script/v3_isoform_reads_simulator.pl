use warnings;
use diagnostics;
if(@ARGV !=4){
	die "*.pl <transcript.fa> <readlen> <true count list> <depth per melocular>";
}

open $ct,$ARGV[2] || die "cannot open $ARGV[2]:$!";
$ind=0;
while(<$ct>){
    chomp;
    if (!/^-?\d+$/){
        print "not a integer\n";
        exit(1);
    }
    else {
        $trexp[$ind++] = $_ ;
    }
}
close $ct;

open $tr,$ARGV[0] || die "cannot open $ARGV[0]:$!";

open ($out,'>',"perfectreads.fastq") or die $!;

$rlen=$ARGV[1];
$dp= $ARGV[3];
$/='>';
<$tr>;
while(<$tr>){		
	$ln++;
	last if $ln > $ind; # transcript num > truth count num
	chomp;
	next if length $_ < 31;
	@lines=(split /\n+/);
	shift @lines;
	$seq = join '',@lines;
	$tlen = length $seq;
	print $tlen,"\n";

	for ($frgn=0; $frgn < $tlen*$dp*$trexp[$ln-1]/$rlen ; $frgn++) {#for ($frgn=0; $frgn < $tlen*$dp*$trexp[$ln-1]/$rlen +1; $frgn++) {
			$rp1 = int(rand($tlen - $rlen));
			$rp2 =int( rand($tlen - $rlen));
			$read1 = substr($seq, $rp1,$rlen);
			$read2 = substr($seq, $rp2,$rlen);
			$rid++;

			printf $out "\@read%08d_%04dp frag%04d ENSMUST%011d gene_id:ENSMUSG%011d gene_name:gene%05d\n%s\n\+\n",$rid, $rp1,$frgn,$ln,$ln,$ln,$read1;
            print $out 'I' x $rlen,"\n";
            printf $out "\@read%08d_%04dp frag%04d ENSMUST%011d gene_id:ENSMUSG%011d gene_name:gene%05d\n%s\n\+\n",$rid,$rp2,$frgn,$ln,$ln,$ln,$read2;
            print $out 'I' x $rlen,"\n";

		}

}

close $tr;
close $out;
