#!/usr/bin/env perl

use strict;
use warnings;
use lib "/home/rcf-02/yunfeigu/perl5";
use lib "/home/yunfeiguo/projects/SeqMule/lib";
use SeqMule::Utils;
use SeqMule::Parallel;
use Getopt::Std;
use Data::Dumper;

die "Usage: $0 [options] <1.coords 2.coords ...>\n".
" -g <INT>	gap allowance for merging near regions, default:100\n".
" -i <INT>	identity, default 90\n".
" -l <INT>	min len2(length of mapped query), default 50\n".
" -r		reverse chromosome and query order\n".
" -f 		treat files as list of coords file paths\n"
unless @ARGV >= 1;
my %options;
getopts("rg:i:l:f",\%options);
#PARAMETERS
my $gap = $options{g} || 100; #max gap allowed for two alignments to be stiched together
my $tmpdir = "/tmp";
my $bedtools = "bedtools";
my $debug = 1;
my $minIdt = $options{i} || 90; #mininum identity between two sequences in a mapping, in percentage
my $minLen2 = $options{l} || 50; #min length of alignment
my $hpc_qsub = "qsub -S /bin/bash -V -l walltime=1:0:0 -l nodes=1:ppn=1 -l mem=2GB -A lc_kw -q laird ";
my $bcl_qsub = "qsub -S /bin/bash -V -l h_vmem=2G -cwd ";
my $qsub = $bcl_qsub;
my $cwd = $ENV{PWD};

my $steps="
* convert all coords to bed, specifying which part of query is mapped to which part of ref
* merge all these bed files with gap allowance of 1kb
(the following can be done by getComplementFasta.pl)
* get complement of all these bed files
* output sequence larger than 1kb
";

#INPUT
my @cleanQ; #files to be removed
my $outputFile = File::Spec->catfile($tmpdir,rand($$).".mummer.unmerged.bed");
open (my $fh,">",$outputFile) or die "open($outputFile): $!\n";
warn "\@ARGV: @ARGV\n" if $debug == 1;
print $fh "#CONTIG\tSTART\tEND\tCHR\tSTART\tEND\tQUERY_LENGTH\tREF_LENGTH\tIDENTITY\n";
if ($options{f}) {
    my @fileList = `cat $ARGV[0]`;
    chomp @fileList;
    warn "No file\n" and exit unless @fileList;
    while(my $coord = shift @fileList) {
	&coord2bed($fh,$coord);
    }
} else {
    while(my $coord = shift @ARGV) {
	&coord2bed($fh,$coord);
    }
}
close $fh;
warn "reading coords done\n";
warn "Output written to ",&mergeBED($outputFile,$gap,"all.mapped.bed"),"\n";
warn "Clean up ...\n";
&cleanup();
warn "All done\n";

###################################################################
sub coord2bed {
	#read coords file, store parsed result in a hash
	my $fh = shift;
	my $coord = shift;
	warn "will parse $coord\n" if $debug == 1;
	return unless $coord;
	open IN,'<',$coord or die "open($coord): $!\n";
	while(<IN>){
		s/\|//g;
		s/^[\t ]+//;
		next if /^[=\/\s\[]|^NUCMER/;
		chomp;
		my $parsedLine = &parseLine($_);
		my $id = $parsedLine->{id};
		next unless &passFilter($parsedLine);
		($parsedLine->{query_start},$parsedLine->{query_end}) = &smallerFirst($parsedLine->{query_start},$parsedLine->{query_end});
		($parsedLine->{ref_start},$parsedLine->{ref_end}) = &smallerFirst($parsedLine->{ref_start},$parsedLine->{ref_end});
		$parsedLine->{query_start} -= 1; #convert to 0-based start
		$parsedLine->{ref_start} -= 1;
		#we need to make sure all regions that are accepted as final mapped regions come from
		#a continuous sequence on the chr
		#store original data, only for splitting
		warn "Parsing line $.: $_\n" if $debug == 3;
		warn Dumper($parsedLine) if $debug == 3;
        print $fh join("\t",
            $parsedLine->{id},$parsedLine->{query_start},$parsedLine->{query_end},
            $parsedLine->{chr},$parsedLine->{ref_start},$parsedLine->{ref_end},
	    $parsedLine->{len1},$parsedLine->{len2},$parsedLine->{idt},),"\n";
	}
	close IN;
}
sub passFilter {
	#filter line by global filters
	my $parsedLine = shift;
	my $result = 1;
	$result = 0 if $parsedLine->{idt} < $minIdt;
	$result = 0 if $parsedLine->{len2} < $minLen2;

	return $result;
}
sub parseLine {
	#parse a line in *.coords
	my $line = shift;
	my $result = {};
	my @f=split ' ',$line;
	die "ERROR: expected 13 fields: $line\n" unless @f==13;
	warn "DEBUG:@f\n" if $debug == 3;
#    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
#===============================================================================================================================
#15007271 15007689      32423    32005        419      419      86.64   133797422    51437       0.00     0.81   chr10	m150320_100742_42199_c100794652550000001823158109091525_s1_p0/54078/7931_59368
#0		1	2	3		4	5	6	7		8	  9	  10	  11	12
	#here we only care how much of the query is aligned, regardless of alignment location
	if($options{r}) {
	    $result->{query_start} = $f[0];
	    $result->{query_end} = $f[1];
	    $result->{ref_start} = $f[2];
	    $result->{ref_end} = $f[3];
	    $result->{id} = $f[11];
	    $result->{chr} = $f[12];
	    $result->{idt} = $f[6];
	    $result->{len1} = $f[5]; #alignment length on ref
	    $result->{len2} = $f[4]; #alignment length on query
	} else {
	    $result->{query_start} = $f[2];
	    $result->{query_end} = $f[3];
	    $result->{ref_start} = $f[0];
	    $result->{ref_end} = $f[1];
	    $result->{id} = $f[12];
	    $result->{chr} = $f[11];
	    $result->{idt} = $f[6];
	    $result->{len1} = $f[4]; #alignment length on ref
	    $result->{len2} = $f[5]; #alignment length on query
	}
	return $result;
}
sub cleanup {
	#remove temp files
	unlink @cleanQ;
}
sub smallerFirst {
	my $start = shift;
	my $end = shift;
	if($start > $end) {
		#make sure end is no smaller than start
		$start = $end + $start;
		$end = $start - $end;
		$start = $start - $end;
	}
	return($start,$end);
}
sub mergeBED {
	my $bed = shift;
	my $gap = shift;
	my $out = shift;
	$gap = $gap || 0;
	$out = $out || "$tmpdir/$$".rand($$).".tmp.bed";
	#keep column 4,5,6, only output distinct columns
	my $nCol = `perl -ne '\@f=split;print scalar \@f and exit;' $bed`;
	!system("$bedtools sort -i $bed | $bedtools merge -d $gap ".($nCol>3? " -c 4,5,6,7,8,9 -o distinct ":"")." > $out") or die "merging $bed fail: $!\n";
	warn("executing: $bedtools sort -i $bed \| $bedtools merge -d $gap ".($nCol>3?" -c 4,5,6,7,8,9 -o distinct ":"")." > $out\n") if $debug >= 1;
	return $out;
}
