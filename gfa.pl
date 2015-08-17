#!/usr/bin/env perl
use strict;
use warnings;
use lib "/home/yunfeiguo/projects/SeqMule/lib";
use SeqMule::Utils;

my $FLANKSIZE = 100;
my $nproc = 12;
my $filter = {qcov=>95,idt=>95,e=>1};
my $debug = 1;


my $fa = "GRCh38full.fa"; 
my $gap = "$fa.gap.bed";
my $genome = "$fa.genome";
my $flankBed = "$fa.gap.flank.bed";
my $flankFa = "$fa.gap.flank.fa";
my $db = "/home/yunfeiguo/projects/PacBio_reference_genome/falcon_aln/hx1_20150716/2-asm-falcon/hx1f4full.fa";
my $total = 0;
my $unfilled = 0;

if(0) {
    !system("./extract_gaps.pl $fa > $gap") or die "Gap profiling failed\n" or die "$!\n";
    !system("/home/yunfeiguo/scripts/seq_related/fa2size.pl $fa > $genome") or die "Genome file creation failed\n";
    #prepare blast index
    !system("makeblastdb -in $db -dbtype nucl") or die "$!\n";
    !system("makembindex -input $db -iformat fasta") or die "$!\n";
}
if(1) {
    for my $onegapRecord(&readBED($gap)) {
	my $onegap = &SeqMule::Utils::genBED([$onegapRecord]);
	my $gapid = join("_",@$onegapRecord);
	my $oneflank = "$onegap.flank.bed";
	my $oneflankfa = "$oneflank.fa";
	my $oneresult = "$oneflank.mapping";
	!system("bedtools flank -i $onegap -g $genome -b $FLANKSIZE > $oneflank") or die "Extraction of flanking sequences failed: $!\n";
	!system("bedtools getfasta -fi $fa -bed $oneflank -fo $oneflankfa") or die "BED2FA failed: $!\n";
	warn "$oneflank" if $debug;
	warn "$oneflankfa" if $debug;
	!system("blastn -db $db -query $oneflankfa -num_threads $nproc -task megablast -max_target_seqs 1 -outfmt '10 qseqid qlen sseqid qstart qend sstart send sstrand evalue pident qcovs' > $oneresult") or die "blast failed: $!\n";
	warn "mapping for $gapid: $oneresult\n" if $debug;
	my @filledGap = &fillGap({source=>$fa,id=>$gapid,filter=>$filter,mapping=>$oneresult,db=>$db,gapRecord=>$onegapRecord});
	if(@filledGap) {
	    &outputFilledGap(@filledGap);
	} else {
	    $unfilled ++;
	    warn "$gapid cannot be filled\n";
	}
	$total++;
    }
}
warn "$unfilled/$total remains unfilled\n";

#-----------------------------------------------------------------
sub outputFilledGap {
    #take a lit of filled gaps and their names (gaps that are filled), scores, strands
    #print BED
    for my $i(@_) {
	print join("\t",@$i),"\n";
    }
}
sub fillGap {
    #take mapping result, gap record and database from which gap sequences come
    #return sequences expected to fill the gaps
    my $desc = "
    * filtering by qcov and idt
    * check if only one flanking region has a mapping, if so, go to gap extension mode
   	if yes
	* try to extend flanking sequence into the gaps
	* return if gaps can be fully recovered
	if no
	* for each pair of mappings between two flanking seq, determine if they concordant (distance, strand)
		if yes
		* fill a gap
		if no
		* skip this mapping
    ";
    my $opt = shift;
    my $mapping = &parseMapping($opt->{mapping});
    my $db = $opt->{db};
    my $source = $opt->{source};
    my $source_fai = "$source.fai";
    my $fai = "$db.fai";
    my $gapRecord = $opt->{gapRecord};
    my $filter = $opt->{filter};
    my $id = $opt->{id};
    my $filledGap = "";
    my %sourceProfile = &SeqMule::Utils::readFastaIdx($source_fai);
    $mapping = &filterMapping($filter,$mapping);
    if(&ifTwinMapping($mapping)) {
	warn  "twin mapping\n" if $debug;
	if(1) {

	} else {

	}
    } else {
	warn  "not twin mapping\n" if $debug;
	#if gap occurs at an edge of the target assembly contig
	#then this is possible, otherwise no need to further check
	if($gapRecord->[2] == $sourceProfile{$gapRecord->[0]}->{length} or 
	    $gapRecord->[1] == 0) {
	    return(&extendGap({gapRecord=>$gapRecord,mapping=>$mapping,fai=>$fai,id=>$id}));
	} else {
	    return;
	}
    }
}
sub extendGap {
    my $opt = shift;
    my $gapRecord = $opt->{gapRecord};
    my $mapping = $opt->{mapping};
    my $fai = $opt->{fai};
    my $id = $opt->{id};
    my %targetProfile = &SeqMule::Utils::readFastaIdx($fai);
    die "gap record illegal\n" unless $gapRecord->[2]>$gapRecord->[1]; 
    my $gapLen = $gapRecord->[2]-$gapRecord->[1];
    my @results;
    for my $i(@$mapping) {
	my @oneresult;
	my $status; #full or partial gap closing
	my $strand; #+ or -, used in BED
	modify this part!
	this is only called when gap is on an edge

	if($i->{strand} eq 'plus') {
	    $strand = "+";
	    push @oneresult,$i->{sid},($i->{send}-1);
	    if($targetProfile{$i->{sid}}->{length} >= $i->{send} + $gapLen) {
		#gap can be fully filled
		$status = "full";
	    push @oneresult,($i->{send}+$gaplen);
	    } else {
		#gap can be partially filled
		$status = "partial";
	    push @oneresult,($targetProfile{$i->{sid}}->{length});
	    }
	    push @oneresult,$id,$i->{e},"+",$status;
	} elsif ($i->{strand} eq 'minus') {
	    $strand = "-";
	    push @oneresult,$i->{sid},($i->{sstart}-1),$i->{send},$id,$i->{e};
	    if($i->{send} - $gapLen >= 0) {
		#gap can be fully filled
		$status = "full";
	    } else {
		#gap can be partially filled
		$status = "partial";
	    }
	} else {
	    die "unknown strand".$i->{strand}."\n";
	}
	push @oneresult,$id,$i->{e},$strand,$status;
	push @results,\@oneresult;
    }
    return(@results);
}
sub ifTwinMapping {
    my $mapping = shift;
    my %id;
    for my $i(@$mapping) {
	$id{$i->{qid}} = 1;
    }
    my $idCount =int(keys %id);
    if($idCount == 2) {
	return(1);
    } elsif($idCount == 1 or $idCount == 0) {
	return(0);
    } else {
	die "Only 0,1,2 query IDs are allowed, found $idCount\n";
    }
}
sub filterMapping {
#take a ref to mapping array and return mappings that pass filter
    my $filter = shift;
    my $mapping = shift;
    my @result;

    for my $i(@$mapping) {
	if($filter->{'idt'}) {
	    push @result,$i if $i->{idt} >= $filter->{idt};
	}
	if($filter->{'qcov'}) {
	    push @result,$i if $i->{qcov} >= $filter->{qcov};
	}
	if($filter->{'e'}) {
	    push @result,$i if $i->{e} <= $filter->{e};
	}
    }
    return(\@result);
}
sub parseMapping {
#take raw mapping results, and return hash refs
    my $f = shift;
    my @results;
    open IN,'<',$f or die "open($f): $!\n";
    while(<IN>) {
	my @f = split /,/;
	die "expect 11 fields at line $. of $f\n" unless @f==11;
	#s is for subject, the reference seq
	push @results,{
	    qid	=>	$f[0],
	    qlen	=>	$f[1],
	    sid	=>	$f[2],
	    qstart	=>	$f[3],
	    qend	=>	$f[4],
	    sstart	=>	$f[5],
	    send	=>	$f[6],
	    strand	=>	$f[7],
	    e		=>	$f[8],
	    idt	=>	$f[9],
	    qcov	=>	$f[10],
	};
    }
    close IN;
    return \@results;
}
sub readBED {
    my $bed=shift;
    warn "reading $bed\n";
    my @out;
    open IN,'<',$bed or die "Cannot read $bed: $!\n";
    while(<IN>) {
	chomp;
	next if /^browser|^#|^track/;
	warn "Can't recognize line $.\n" and next unless /^([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)/;
	if ($3>$2) {
	    push @out,[$1,$2,$3,$4];
	} else {
	    push @out,[$1,$3,$2,$4];
	}
    }
    close IN;
    return @out;
}
