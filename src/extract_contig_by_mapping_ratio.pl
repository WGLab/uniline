#!/usr/bin/env perl
#
use strict;
use warnings;
use Bio::SeqIO;

die "Usage: $0 <mapping_ratio_file> <contig.fasta> [threshold,default 0.1]\n" unless @ARGV == 2 or @ARGV == 3;
my $mr_file = shift @ARGV;
my $fa = shift @ARGV;
my $threshold = shift @ARGV;
$threshold = $threshold || 0.5;
my %target_contig = &getContigID($mr_file,$threshold);

my $seqio  = Bio::SeqIO->new( '-format' => 'fasta' , -file => $fa);
while (my $seqobj = $seqio->next_seq() ) {
    if($target_contig{$seqobj->id}) {
	print '>',$seqobj->id,"\n";
	print $seqobj->seq,"\n" 
    }
}
warn "All done\n";

###############################SUBROUTINES##################################################
sub getContigID {
    my $file = shift;
    my $threshold = shift;
    my %return;
    my $example = "
    FASTA_ID	Mapped_length	Total_length	Mapping_ratio	Mapped_chr
    006328F	1892518	72416	2613.40%	chr10 chr10
    0		1	2	3		4	5 ...
    ";

    open IN,'<',$file or die "open($file): $!\n";
    my $header = <IN>;
    while(<IN>) {
	my @f = split;
	$return{$f[0]} = 1 if $f[3] <= $threshold;
    }
    close IN;
    return %return;
}
