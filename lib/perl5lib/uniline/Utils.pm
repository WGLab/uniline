package uniline::Utils;

use strict;
use warnings;
use Carp qw/carp croak/;
sub readBED {
    my $bed=shift;
    my @out;
    open IN,'<',$bed or croak("Cannot read $bed: $!\n");
    while(<IN>)
    {
	next if /^browser|^#|^track/;
	next unless /^([^\t]+)\t(\d+)\t(\d+)/;
	if ($3>$2)
	{
	    push @out,[$1,$2,$3];
	} else
	{
	    push @out,[$1,$3,$2];
	}
    }
    close IN;
    return @out;
}

sub readFastaIdx
{
    my $fai=shift;
    my %return;
    open IN,"<",$fai or (carp "Failed to read $fai: $!\n" and return undef);
    while (<IN>)
    {
	my @f=split /\t/;
	croak("5 fields expected at line $. of $fai: $_\n") unless @f==5;
	my ($id,$len,$offset,$nchar_ln,$nbyte_ln)=@f;

	$return{$id}={
	    length=>$len, #length of contig
	    offset=>$offset, #offset where first character in that contig appears
	    nchar_ln=>$nchar_ln, #number of characters per line
	    nbyte_ln=>$nbyte_ln, #number of bytes per line
	};
    }
    close IN;
    return %return;
}

sub genBED
{
    my $content=shift;
    my $out=shift;
    $out=$out||"/tmp/$$".rand($$).".tmp.bed";

    open OUT,'>',$out or croak("ERROR: Failed to write to $out: $!\n");
    for my $i(@$content)
    {
	print OUT join "\t",@$i;
	print OUT "\n";
    }
    close OUT;
    return $out;
}
#convert .fai to .genome file
sub fa2size {
    my $file = shift;
    my $genome = shift;
    my $fai = "$file.fai";
    unless(-e $fai or -l $fai) {
	carp("Call SAMtools to create FASTA index first\n");
	!system("samtools faidx $file") or croak("samtools $file indexing failed: $!\n");
    }
    croak("$fai missing!\n") unless -e $fai or -l $fai;
    return(!system("cut -f 1,2 $fai | sort -k 1,1 -k2,2n > $genome"));
}
1;
