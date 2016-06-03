#!/usr/bin/env perl

use strict;
use warnings;
use lib "/home/yunfeiguo/projects/SeqMule_dev/lib";
use SeqMule::Utils;
use threads;
use threads::shared;

die "Usage: $0 <indexed query FASTA> <1.coords 2.coords ...>\n" unless @ARGV >= 2;
#    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
#===============================================================================================================================
#15007271 15007689      32423    32005        419      419      86.64   133797422    51437       0.00     0.81   chr10	m150320_100742_42199_c100794652550000001823158109091525_s1_p0/54078/7931_59368
#0		1	2	3		4	5	6	7		8	  9	  10	  11	12
my $gap = 100; #max gap allowed for two alignments to be stiched together
my $tmpdir = "/tmp";
my $debug = 0;
my $query = shift @ARGV;
my $idx = "$query.fai";
die "no index: $idx\n" unless -e $idx;
my %fa = &SeqMule::Utils::readFastaIdx($idx);
my $total = scalar keys %fa;
my @cleanQ; #files to be removed
for my $i(@ARGV) {
	open IN,'<',$i or die "open($i): $!\n";
	while(<IN>){
		s/\|//g;
		s/^[\t ]+//;
		next if /^[=\/\s\[]|^NUCMER/;
		chomp;
		my @f=split;
		die "ERROR: expected 13 fields: $_\n" unless @f==13;
		warn "DEBUG:@f\n" if $debug;
		#here we only care how much of the query is aligned, regardless of alignment location
		my $query_start = $f[2];
		my $query_end = $f[3];
		my $ref_start = $f[0];
		my $ref_end = $f[1];
		my $id = $f[12];
		my $chr = $f[11];
		($query_start,$query_end) = &smallerFirst($query_start,$query_end);
		($ref_start,$ref_end) = &smallerFirst($ref_start,$ref_end);
		#we need to make sure all regions that are accepted as final mapped regions come from
		#a continuous sequence on the chr
		if(defined $fa{$id}->{queryMapping} ) {
			push @{$fa{$id}->{refMapping}},[$chr,$ref_start,$ref_end];
			push @{$fa{$id}->{queryMapping}},[$chr,$query_start,$query_end];
		} else {
			$fa{$id}->{refMapping} = [[$chr,$ref_start,$ref_end]];
			$fa{$id}->{queryMapping} = [[$chr,$query_start,$query_end]];
		}
	}
	close IN;
}
warn "reading coords done\n";
&convert2BED(\%fa);
warn "conversion to bed done\n";
&mergeBED_extractMaxMapping(\%fa);
warn "merging and extraction done\n";
&output(\%fa);
&cleanup(\%fa);
warn "All done\n";

###################################################################
sub convert2BED {
	my $ref = shift;
	my $count = 0;
	for my $i(keys %$ref) {
		$count++;
		warn "conversion:$count/$total done\n" if $count % 100 == 0;
		next unless defined $ref->{$i}->{refMapping};
		#the mapping may contain overlapping regions
		#we should only use non-overlapped regions
		my $ref_bed = "$tmpdir/$i".rand($$).".ref.bed";
		my $ref_overlapBED = &SeqMule::Utils::genBED($ref->{$i}->{refMapping});
		$ref->{$i}->{refBED} = &mergeBED($ref_overlapBED,0,$ref_bed);
		push @cleanQ,$ref_overlapBED,$ref_bed;
	}
}
sub mergeBED_extractMaxMapping {
	my $ref = shift;
	my $count = 0;
	for my $i(keys %$ref) {
		$count++;
		warn "merge and extractiong: $count/$total done\n" if $count % 100 == 0;
		next unless defined $ref->{$i}->{refBED};
		&findMapping($ref->{$i});
	}
}
sub findMapping {
	my $query_ref = shift;
	my $ref_bed = $query_ref->{refBED};
	my $out = &mergeBED($ref_bed,$gap); #merge with gap allowance
	push @cleanQ,$out;
	#maxMapping = [chr,start,end]
	$query_ref->{maxMapping} = &findMax([&SeqMule::Utils::readBED($out)]);
	$query_ref->{mappedLen} = &getMappedLen($query_ref->{maxMapping},$query_ref->{refMapping},$query_ref->{queryMapping});
}
sub getMappedLen {
	#based on maxMapping found on ref
	#identify mapped regions on query
	#get total length of mapped regions on
	#query
	my $max = shift;
	my $allRef = shift;
	my $allQuery = shift;
	my $chr = $max->[0];
	my $start = $max->[1];
	my $end = $max->[2];
	my $mappedLen = 0;
	my @mappedQuery;
	for my $i(0..$#{$allRef}) {
		#the mapped region must be fully contained
		#as we have merged all mapped regions
		next unless $allRef->[$i]->[0] eq $chr and $allRef->[$i]->[1] >= $start and $allRef->[$i]->[2] <= $end;
		push @mappedQuery,$allQuery->[$i];
	}
	#mapped regions in query should be merged before
	#length is calculated
	my $query_mapped_bed = &SeqMule::Utils::genBED(\@mappedQuery);
	my $query_mapped_nonoverlap_bed = &mergeBED($query_mapped_bed);
	$mappedLen = &SeqMule::Utils::bed2total($query_mapped_nonoverlap_bed);

	push @cleanQ,$query_mapped_bed,$query_mapped_nonoverlap_bed;
	return $mappedLen;
}
sub findMax {
	#given a set of regions in BED format
	#find the longest the region
	my $ref = shift;
	my $max = 0;
	my $pos = [];
	#array of [chr,start,end]
	for my $i(@$ref) {
		my $len = $i->[2]-$i->[1];
		#if there are multiple alignments with equal lengths,
		#we output first one of them
		if ($len > $max) {
			$pos = $i;
			$max = $len;
		}
	}
	return $pos;
}
sub output {
	my $ref = shift;
	my %fa = %$ref;
	print join("\t","FASTA_ID","Mapped_length","Total_length","Mapping_ratio","Mapped_chr","Query_start","Query_end"),"\n";
	for my $i(keys %fa) {
		my $mappedLen = $fa{$i}->{mappedLen} || 0;
		my $totalLen = $fa{$i}->{length};
		my $mapping = ['NA','NA','NA'];
		if(defined $fa{$i}->{maxMapping}) {
			$mapping = $fa{$i}->{maxMapping};
		}
		my $mappingRatio = $totalLen == 0? 0:$mappedLen/$totalLen;
		print join("\t",$i,
			$mappedLen,
			$totalLen,
			sprintf("%.2f",$mappingRatio),
			@$mapping,
		), "\n";
	}
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
	!system("bedtools sort -i $bed | bedtools merge -d $gap > $out") or die "merging $bed fail: $!\n";
	return $out;
}
sub joinThreads {
	for my $i(threads->list(threads::running))
	{
		$i->join(); 
	}
}
