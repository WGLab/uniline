#!/usr/bin/env perl
use strict;
use Getopt::Std;
use warnings;


die "Usage: $0 [-r] <bwa|blasr> <read length> <interval size> <prefix> <reference.fa> <query.fa>\n".
"-r	random interval of specified size\n" unless @ARGV >= 6;
my %opts;
getopts('r',\%opts);
my $aligner = shift @ARGV;
my $readLen = shift @ARGV;
my $interval = shift @ARGV;
my $prefix = shift @ARGV;
my $ref = shift @ARGV;
my $query = shift @ARGV;
warn "NOTICE: processing $query against $ref\n";
my $queryGenome = "$query.fai";
my $queryFai = "$query.genome";
die "interval length must be larger than 2*read length\n" unless $interval > 2*$readLen;
unless (-e $queryFai) {
    !system("samtools faidx $query") or die "samtools: $!\n";
}
unless (-e $queryGenome) {
    !system("cut -f 1,2 $queryFai > $queryGenome") or die "cut: $!\n";
}
###based on genome size file
##generate bed file specifying 5Mb windows
#start of window is start of first 1000bp read, end of window is start of second 1000bp read
$prefix .= "_interval${interval}_read${readLen}_aligner${aligner}";
my $out = "${prefix}misjoiningRate.txt";
my ($windowBed, $readBed, $readFasta, $bam) = ($prefix."windows.bed", $prefix."read.bed", $prefix."read.fasta", $prefix.".bam");
if ($opts{r}) {
    !system("bedtools random -l $interval -n 100000  -g $queryGenome > $windowBed") or die;
    !system("perl -ane 'next unless \$F[2] - \$F[1] == $interval;print join(\"\\t\",\$F[0],\$F[1],\$F[1]+$readLen,\"\$F[0]window\$F[3]read1\"),\"\\n\";print join(\"\\t\",\$F[0],\$F[2] - $readLen,\$F[2],\"\$F[0]window\$F[3]read2\"),\"\\n\";' $windowBed > $readBed") or die;
} else {
    !system("bedtools makewindows -g $queryGenome -w ".($interval - $readLen)." -s $interval -i winnum > $windowBed") or die;
    !system("perl -ane 'next unless \$F[2] - \$F[1] == $interval - $readLen;print join(\"\\t\",\$F[0],\$F[1],\$F[1]+$readLen,\"\$F[0]window\$F[3]read1\"),\"\\n\";print join(\"\\t\",\$F[0],\$F[2],\$F[2] + $readLen,\"\$F[0]window\$F[3]read2\"),\"\\n\";' $windowBed > $readBed") or die;
}

###based on window coordinates, create BED for read extraction
###extract sequences to create artificial reads
!system("bedtools getfasta -fi $query -bed $readBed -fo $readFasta -name") or die;
#map with BWA
if ($aligner eq 'bwa') {
    !system("bwa mem -t 6 -v 1 $ref $readFasta 2>/dev/null | samtools view -@ 2 -bS - | samtools sort -@ 2 -f -n - $bam") or die; 
} elsif ($aligner eq 'blasr') {
    !system("blasr -nproc 6 -sam $readFasta $ref 2>/dev/null | samtools view -bS - | samtools sort -f -n - $bam") or die; 
} else {
    die "unknown aligner: $aligner\n";
}
#calculate number of inconsistent reads
#allow 100kb discrepancy
#source ./python_virtualenv/bin/activate
!system("calcMisjoinedScaffold.py ".($interval - 0.1*$interval)." ".($interval + 0.1 * $interval)." $bam 2>$out") or die "$!";
