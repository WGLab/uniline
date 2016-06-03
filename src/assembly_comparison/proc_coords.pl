#!/usr/bin/env perl

use strict;
use warnings;

#    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
#===============================================================================================================================
#15007271 15007689      32423    32005        419      419      86.64   133797422    51437       0.00     0.81   chr10	m150320_100742_42199_c100794652550000001823158109091525_s1_p0/54078/7931_59368
#0		1	2	3		4	5	6	7		8	  9	  10	  11	12
die "Usage: $0 <1.coords 2.coords ...>\n" unless @ARGV >= 1;
my %fa=();
for my $i(@ARGV) {
    open IN,'<',$i or die "open($i): $!\n";
    while(<IN>){
	s/\|//g;
	s/^[\t ]+//;
	next if /^[=\/\s\[]|^NUCMER/;
	chomp;
	my @f=split;
	my $id = $f[12];
	my $chr = $f[11];
	my $len = $f[8];
	my $mapped_len = $f[5];
	#reverse ref and query
	#my $id = $f[11];
	#my $chr = $f[12];
	#my $len = $f[7];
	#my $mapped_len = $f[4];

	if(defined $fa{$id} ) {
	    push @{$fa{$id}->{chr}},$chr;
	    $fa{$id}->{mapped_len} += $mapped_len;
	} else {
	    $fa{$id} = {};
	    $fa{$id}->{chr} = [$chr];
	    $fa{$id}->{mapped_len} = $mapped_len;
	    $fa{$id}->{len} = $len;
	}
    }
    close IN;
}
print join("\t","FASTA_ID","Mapped_length","Total_length","Mapping_ratio","Mapped_chr"),"\n";
for my $i(keys %fa) {
    my @chr = @{$fa{$i}->{chr}};
    print join("\t",$i,
	$fa{$i}->{mapped_len},
	$fa{$i}->{len},
	#sprintf("%.2f",$fa{$i}->{mapped_len}/$fa{$i}->{len}*100)."%",
	sprintf("%.2f",$fa{$i}->{mapped_len}/$fa{$i}->{len}),
	join(" ",@chr)
    ), "\n";
}
