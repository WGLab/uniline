#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

die "Usage: $0 <fasta>\n" unless @ARGV==1;
warn "Scanning for gaps...\n";
my $fa = shift;
my $seqio_obj = Bio::SeqIO->new(-file=>$fa,-format=>'fasta');
my $count = 0;
while(my $i = $seqio_obj->next_seq) {
    my $pos = 0; #position of probe
    my $ingap = 0; #whether we are in a gap
    for my $j(split //,$i->seq) {
	if ($j eq 'N' or $j eq 'n') {
	    print $i->id,"\t$pos\t" if $ingap == 0;
	    $ingap = 1;
	} else {
	    if($ingap) {
		$ingap = 0;
		print $pos,"\tGap$count\n";
		$count++;
	    } else {
		1;
	    }
	}
	$pos++;
    }
    if($ingap) {
	print $pos,"\tGap$count\n";
	$count++;
    }
}
