#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/basename/;
use Getopt::Std;

my %options;
getopts("l:",\%options);
die "Usage: $0 [options] <bed> <indexed genome fasta>\n".
" -l <INT>	min length for output,default 10,000\n"
unless @ARGV == 2;
my $minLen = $options{l}||10_000; #min length for unmapped region
my $bedtools = "bedtools";
my $bed = shift;
my $genome = "/tmp/".rand($$).".genome";
my $sorted_bed = "/tmp/".rand($$).".sorted.bed";
my $fa = shift;
my $fai = "$fa.fai";
my $outFa = basename $bed.".complement.$minLen.fasta";
my $complementBED = basename $bed.".complement.$minLen.bed";

!system("cut -f 1,2 $fai | sort -k 1,1 -k2,2n > $genome") or die "cut $fai: $!\n";
!system("sort -k 1,1 -k2,2n $bed > $sorted_bed") or die "bedtools sort: $!\n";
!system("$bedtools complement -i $sorted_bed -g $genome | perl -ne '\@f=split;print if \$f[2]-\$f[1] >= $minLen' > $complementBED") or die "bedtools complement: $!\n";
!system("$bedtools getfasta -fi $fa -bed $complementBED -fo $outFa") or die "bedtools getfasta: $!\n";
