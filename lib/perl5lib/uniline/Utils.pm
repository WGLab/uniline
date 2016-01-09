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
sub extractGap {
    carp("Scanning for gaps...\n");
    my $fa = shift;
    my $out = shift;
    my $count = 0;
    open IN,"<",$fa or die "ERROR: open() failed to read $fa $!\n";
    open OUT,">",$out or die "ERROR: open() failed to write $out $!\n";
    select OUT;
    my ($id,$ingap,$pos,$start,$end,$cur);
    while(my $l = <IN>) {
	if ($l =~ /^>(\S+)/) {
	    if ($ingap) {
		$end = $pos;
		print "$id\t$start\t$end\tGap$count\n";
	    }
	    $id = $1;
	    $pos = 0;
	} else {
	    my $len = length($l);
	    my $i = 0;
	    while ($i < $len) {
		$pos++;
		$cur = substr($l,$i,1);
		if($cur eq 'N' or $cur eq 'n') {
		    if (not $ingap) {
			#get a new gap
			$count++;
			$ingap = 1;
			$start = $pos - 1; #zero start in BED file
		    }
		} elsif ($ingap) {
		    $ingap = 0;
		    $end = $pos - 1;
		    print "$id\t$start\t$end\tGap$count\n";
		} 
		$i++;
	    }
	}
    }
    if ($ingap) {
	$end = $pos;
	print "$id\t$start\t$end\tGap$count\n";
    }
    close IN;
    close OUT;
    carp("NOTICE: all gaps extracted to $out\n");
}
1;
