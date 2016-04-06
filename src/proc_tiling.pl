#!/usr/bin/env perl
#
use strict;
use warnings;
#>chrX 156040895 bases
#22904	32912	213229	10009	100.00	95.81	+	011957F
#246142	253576	8706	7435	100.00	94.34	-	012847F
#262283	262644	15303	362	100.00	98.10	-	016256F
#ref_start ref_end gap	contig_len cov  idt	strand  contigID
#0	1	2	3	4	5	6	7
my (%chr,$current);
die "Usage: $0 <1.tiling 2.tiling ...>\n" unless @ARGV >= 1;
for my $i(@ARGV) {
    open IN,'<',$i or die "open($i): $!\n";
    while(<IN>) {
	chomp;
	if(/^>/) {
	    if(/^>(\w+)\s+(\d+)\s+bases/) {
		$current = $1;
		$chr{$1} = {};
		$chr{$1}->{len} = $2;
	    } else {
		die "unrecognized header:$_\n";
	    }
	} else {
	    my @f = split;
	    if($f[1] > $chr{$current}->{len}) {
		#sometimes alignment end may go beyond this chromosome due to query
		#contig is too long
		$f[1] = $chr{$current}->{len};
	    }
	    if($f[0] < 0) {
		$f[0] = 0;
	    }
	    if(defined $chr{$current}->{cover}) {
		$chr{$current}->{cover} += $f[1]-$f[0];
	    } else {
		$chr{$current}->{cover} = $f[1]-$f[0];
	    }
	}
    }
    close IN;
}

#output
print "ID\tLENGTH\tCOVER_LENGTH\tCOVERAGE\n";
for my $i(sort keys %chr) {
    my $len = $chr{$i}->{len};
    my $cover = $chr{$i}->{cover};
    my $cov_ratio = $cover>$len? 1:$cover/$len;
    print join("\t",$i,$len,$cover,$cov_ratio),"\n";
}
