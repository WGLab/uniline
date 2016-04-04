#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std qw/getopts/;

my %seen;
my %opts;
getopts('l:',\%opts);
die "Usage: $0 <GFA output BED>\n".
" -l	min length of original gap\n" unless @ARGV==1;
my $input = shift @ARGV;
open IN,'<',$input or die "open($input): $!\n";
while(<IN>) {
    #example line
    #000267F-010-01	0	432	chr1_297968_347968_Gap3	45	-	extension_partial
    #0			1	2	3			4	5	6
    my @f=split;
    my $l=$f[2]-$f[1]; #length of gap on source assembly
    my @g=($f[3]=~/(.*?)_(\d+)_(\d+)_(\w+\d+)/); #interpret original gap
    my @orig=@g;
    $g[2]=$g[1]+$l unless /full/; #if not fully closed, then add partial closed length to start point to get end point
    if(defined $seen{$f[3]}) {
	$seen{$f[3]} = 0; #the gap is not uniquely closed, skip it
    } else {
	$seen{$f[3]} = 0; #the gap is not uniquely closed, skip it
	if(defined $opts{l}) {
	    next unless $orig[2]-$orig[1] >= $opts{l};
	} 
	$seen{$f[3]} = join("\t",@g[0..2],$f[4],$g[3],$f[6]);
    }
}
close IN;
my ($total,$outputCount);
for my $i(keys %seen) {
    $total++;
    if($seen{$i} ne '0') {
	$outputCount++;
	print $seen{$i},"\n";
    }
}
warn "Output $outputCount/$total results\n";
warn "All done\n";
