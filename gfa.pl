#!/usr/bin/env perl
use strict;
use warnings;
use lib "/home/yunfeiguo/projects/SeqMule/lib";
use SeqMule::Utils;
use Data::Dumper;

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
	#then only one sequence can be mapped
	#when the gap doesn't occur the at edge, it is still possible to have only one sequence mapped
	#if($gapRecord->[2] == $sourceProfile{$gapRecord->[0]}->{length} or 
	#    $gapRecord->[1] == 0) {
	#    return(&extendGap({gapRecord=>$gapRecord,mapping=>$mapping,fai=>$fai,id=>$id}));
	#} else {
	#    return;
	#}
	return(&extendGap({gapRecord=>$gapRecord,mapping=>$mapping,fai=>$fai,id=>$id}));
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
	my ($ctg, $start, $end,
	    $name,$score,
	    $strand,#+ or -, used in BED
	    $status,#full or partial gap closing
	);
	$ctg = $i->{sid};
	$name = $id;
	$score = $i->{e};
	($start,$end,$strand,$status) = &decideGapLocus($targetProfile{$i->{sid}},$i,$gapRecord,$gapLen);
	push @oneresult,$ctg, $start, $end, $name,$score, $strand, $status;
	push @results,\@oneresult;
    }
    return(@results);
}
sub decideGapLocus {
    my $ctgInfo = shift;
    my $mapping = shift;
    my $gapRecord = shift;
    my $gapLen = shift;
    my ($start,$end,$strand,$status);
    &fixMapping($mapping);
    sub fixMapping {
	my $mapping = shift;
	if($mapping->{
    }
    if($mapping->{strand} eq 'plus') {
	#* for N, |----> for mapping
	#case1
	#---*******|--->-------------------------
	#case2
	#----------------------|--->**********---
	$strand = "+";
	if(&isUpstream($mapping,$gapRecord)) {
	    #case2
	    $start = $mapping->{send}; #-1 is not necessary, this is end for flanking seq
	    if($ctgInfo->{length} >= $mapping->{send} + $gapLen) {
		#gap can be fully filled
		$status = "full";
		$end = $start + $gapLen;
	    } else {
		#gap can be partially filled
		$status = "partial";
		$end = $ctgInfo->{length};
	    }
	} else {
	    #case1
	    $end = $mapping->{sstart}-1; #-1 is necessary
	    if($end - $gapLen >= 0) {
		#gap can be fully filled
		$status = "full";
		$start = $end - $gapLen;
	    } else {
		#gap can be partially filled
		$status = "partial";
		$start = 0;
	    }
	}
    } elsif ($mapping->{strand} eq 'minus') {
	#case3
	#---*******<---|-------------------------
	#case4
	#----------------------<---|**********---
	$strand = "-";
	if(&isUpstream($mapping,$gapRecord)) {
	    #case2
	    $end = $mapping->{sstart}; #-1 is not necessary, this is end for flanking seq
	    if($ctgInfo->{length} >= $end + $gapLen) {
		#gap can be fully filled
		$status = "full";
		$start = $end + $gapLen;
	    } else {
		#gap can be partially filled
		$status = "partial";
		$start = $ctgInfo->{length};
	    }
	} else {
	    #case1
	    $start = $mapping->{send}-1; #-1 is necessary
	    if($start - $gapLen >= 0) {
		#gap can be fully filled
		$status = "full";
		$end = $start - $gapLen;
	    } else {
		#gap can be partially filled
		$status = "partial";
		$end = 0;
	    }
	}
    } else {
	die "unknown strand:".$mapping->{strand}."\n";
    }
    if($start == $end) {
	for my $i(keys %{$mapping}) {
	    print "key:$i,value: ",$mapping->{$i},"\n";
	}
	die "equal start and end !\n";
    }
    return($start,$end,$strand,$status);
}
sub isUpstream {
    my $mapping = shift;
    my $gap = shift;
    #chr1:10000-10100
    my ($start,$end) = $mapping->{qid} =~ /^[^:]+:(\d+)-(\d+)$/ or 
    die "unrecognized query ID: ",$mapping->{qid},"\n";
    if($end < $gap->[2]) {
	return(1);
    } else {
	return(0);
    }
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
