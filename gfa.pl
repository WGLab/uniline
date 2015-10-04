#!/usr/bin/env perl
use strict;
use warnings;
use lib "/home/yunfeiguo/projects/SeqMule/lib";
use SeqMule::Utils;
use File::Basename qw/basename/;
use Data::Dumper;
use Getopt::Std;

#param
my $FLANKSIZE = 100; #length of flank sequences for anchoring
my $MIN_GAP_DIST = 2*$FLANKSIZE; #min distance allowed between two gaps
my $GAPLEN_MIS_TOLERANCE = 100; #max difference between filled gap and original gap
my $nproc = 12; #nproc for blast
my $filter = {qcov=>95,idt=>95,e=>1e-8}; #mapping filters
my $debug = 1;
#external exe
my $fa2size = "/home/yunfeiguo/scripts/seq_related/fa2size.pl";
my $extractGap = "/home/yunfeiguo/scripts/seq_related/extract_gaps.pl";


my %opts;
getopts('s',\%opts);
die "Usage: $0 <target assembly> <source assembly>\n".
" -s	skip preprocessing (locating gaps, index db)\n"
unless @ARGV == 2;
my $fa = shift @ARGV;
my $db = shift @ARGV;
my $gap = (basename $fa).".gap.bed";
my $genome = (basename $fa).".genome";
my $flankBed = (basename $fa).".gap.flank.bed";
my $flankFa = (basename $fa).".gap.flank.fa";
#my $db = "/home/yunfeiguo/projects/PacBio_reference_genome/falcon_aln/hx1_20150716/2-asm-falcon/hx1f4full.fa";
my $total = 0;
my $unfilled = 0;

unless($opts{s}) {
    my $bed_for_merge = "/tmp/".rand($$).".tmp.bed.formerge";
    !system("$extractGap $fa > $gap") or die "Gap profiling failed\n" or die "$!\n";
    !system("cat $gap > $bed_for_merge") and !system("bedtools merge -i $bed_for_merge -d $MIN_GAP_DIST > $gap") or die "failed to merge: $!\n";
    #after merging, gap numbers will be gone!
    !system("perl -i -ne '\@f=split;push \@f,\"Gap\$.\";print join(\"\\t\",\@f),\"\\n\"' $gap") or die "failed to add Gap#: $!\n";
    !system("$fa2size $fa > $genome") or die "Genome file creation failed\n";
    #prepare blast index
    !system("makeblastdb -in $db -dbtype nucl 1>&2") or die "makeblastdb: $!\n";
    !system("makembindex -input $db -iformat blastdb 1>&2") or die "makembindex: $!\n";
}
if(1) {
    for my $onegapRecord(&readBED($gap)) {
	#next unless $onegapRecord->[0] =~ /chr17_KI270729v1_random/;
	my $onegap = &SeqMule::Utils::genBED([$onegapRecord]);
	my $gapid = join("_",@$onegapRecord);
	my $oneflank = "$onegap.flank.bed";
	my $oneflankfa = "$oneflank.fa";
	my $oneresult = "$oneflank.mapping";
	!system("bedtools flank -i $onegap -g $genome -b $FLANKSIZE > $oneflank") or die "Extraction of flanking sequences failed: $!\n";
	!system("bedtools getfasta -fi $fa -bed $oneflank -fo $oneflankfa") or die "BED2FA failed: $!\n";
	warn "$oneflank\n" if $debug;
	warn "$oneflankfa\n" if $debug;
	!system("blastn -db $db -query $oneflankfa -num_threads $nproc -task megablast -max_target_seqs 1 -outfmt '10 qseqid qlen sseqid slen qstart qend sstart send sstrand evalue pident qcovs' > $oneresult") or die "blast failed: $!\n";
	warn "mapping for $gapid: $oneresult\n" if $debug;
	my @filledGap = &fillGap({source=>$fa,id=>$gapid,filter=>$filter,mapping=>$oneresult,db=>$db,gapRecord=>$onegapRecord});
	if(@filledGap) {
	    &outputFilledGap(@filledGap);
	} else {
	    $unfilled ++;
	    warn "NOTICE: $gapid cannot be filled\n";
	}
	$total++;
    }
}
warn "NOTICE: $unfilled/$total remains unfilled\n";
#-----------------------------------------------------------------
sub outputFilledGap {
    #take a lit of filled gaps and their names (gaps that are filled), scores, strands, status, upstream/downstream
    #print BED
    for my $i(@_) {
	if($i->[5] eq '-') {
	    #on negative strand, switch start and end
	    my $tmp = $i->[2];
	    $i->[2] = $i->[1];
	    $i->[1] = $tmp;
	}
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
	warn  "DEBUG:twin mapping\n" if $debug;
	return(&twoAnchoredGap({gapRecord=>$gapRecord,mapping=>$mapping,fai=>$fai,id=>$id}));
    } else {
	warn  "DEBUG:not twin mapping\n" if $debug;
	#if gap occurs at an edge of the target assembly contig
	#then only one sequence can be mapped
	#when the gap doesn't occur the at edge, it is still possible to have only one sequence mapped
	return(&extendGap({gapRecord=>$gapRecord,mapping=>$mapping,fai=>$fai,id=>$id}));
    }
}
sub twoAnchoredGap {
    my $opt = shift;
    my $gapRecord = $opt->{gapRecord};
    my $mapping = $opt->{mapping};
    my $fai = $opt->{fai};
    my $id = $opt->{id};
    my %targetProfile = &SeqMule::Utils::readFastaIdx($fai);
    die "gap record illegal\n" unless $gapRecord->[2]>$gapRecord->[1]; 
    my $gapLen = $gapRecord->[2]-$gapRecord->[1];
    my ($anchor1mapping,$anchor2mapping) = &splitMappingbyAnchor($mapping);
    my @results;
    my $count = 0;
    #iterate through all possible pairs to look for concordant pairs
    for my $i(@$anchor1mapping) {
	for my $j(@$anchor2mapping) {
	    next if &isDiscordantAnchor($i,$j,$gapLen);
	    $count++;
	    push @results,();
	    my ($ctg, $start, $end,
		$name,$score,
		$strand,#+ or -, used in BED
		$status,#full or partial gap closing
		$upDownstream, #fill from upstream or downstream?
	    );
	    $ctg = $i->{sid};
	    $name = $id;
	    $score = $i->{e};
	    $strand = ($i->{strand} eq 'plus'? '+':'-');
	    ($start,$end) = &decideDoubleAnchoredGapLocus($i,$j);
	    &checkOrientation($start,$end,$i,$j);
	    $status = "bridge_full";
	    $upDownstream = "NA"; #for double anchors, there is no up down stream
	    if(defined $start && defined $end) {
		push @results, [$ctg, $start, $end, $name,$score, $strand, $status];
	    }
	}
    }
    warn "DEBUG: $count pairs of concordant anchors!\n" if $debug;
    return @results;
}
sub decideDoubleAnchoredGapLocus {
    my $i = shift;
    my $j = shift;
    my $dist1 = abs($i->{sstart} - $j->{send});
    my $dist2 = abs($i->{send} - $j->{sstart});
    my ($start,$end);
    if($dist1 < $dist2) {
	#+ strand
	#           ----->		    ----->
	#-----------------*****************>----------------           
	#- strand
	#           <-----		    <-----
	#-----------------<*****************----------------           
	if ($i->{strand} eq 'plus') {
	    $start = $j->{send};
	    $end = $i->{sstart}-1;
	} else {
	    $start = $j->{send};
	    $end = $i->{sstart} + 1;
	}
    } elsif ($dist1 > $dist2) {
	if ($i->{strand} eq 'plus') {
	    $start = $i->{send};
	    $end = $j->{sstart}-1;
	} else {
	    $start = $i->{send};
	    $end = $j->{sstart} + 1;
	}
    } else {
	warn "ERROR: Inner and outer distances are the same for the following mappings:\n";
	&dumpMapping($i);
	&dumpMapping($j);
	die "\n";
    }
    return($start,$end);
}
sub checkOrientation {
    my $start = shift;
    my $end = shift;
    my $mapping1 = shift;
    my $mapping2 = shift;
    if($mapping1->{strand} eq 'plus' and $start >= $end) {
	warn "ERROR: start($start) must be smaller than end($end) on plus strand\n";
	&dumpMapping($mapping1);
	&dumpMapping($mapping2);
	die;
    } elsif($mapping1->{strand} eq 'minus' and $start <= $end) {
	warn "ERROR: start($start) must be larger than end($end) on minus strand\n";
	&dumpMapping($mapping1);
	&dumpMapping($mapping2);
	die;
    }
}
sub isDiscordantAnchor {
    my $i = shift; #one of anchor mapping
    my $j = shift; #one of the other anchor mapping
    my $gapLen = shift;
    my $discordant = 0;
    my ($chr1,$start1,$end1) = &parseQid($i->{qid});
    my ($chr2,$start2,$end2) = &parseQid($j->{qid});

    if($i->{strand} ne $j->{strand} or  #mapped to different strands
	$i->{sid} ne $j->{sid}) { #mapped to different ctg
	$discordant = 1;
	return($discordant);
    }
    if($i->{strand} eq 'plus') {
	#2nd anchor should not be upstream of 1st anchor
	if( ($start1-$start2) * ($i->{sstart} - $j->{sstart}) > 0) {
	    1;
	} else {
	    $discordant = 1;
	}
    } else {
	#1st anchor should not be upstream of 2nd anchor
	if( ($start1-$start2) * ($i->{sstart} - $j->{sstart}) < 0) {
	    1;
	} else {
	    $discordant = 1;
	}
    }
    #anchors cannot overlap or next to each other
    if( ($i->{strand} eq 'plus' && $i->{sstart} <= $j->{send} && $i->{send} >= $j->{sstart} ) or
	($i->{strand} eq 'minus' && $i->{sstart} >= $j->{send} && $i->{send} <= $j->{sstart}) or
        ($i->{strand} eq 'plus' && ($i->{sstart} - $j->{send} == 1 || $j->{sstart} - $i->{send} == 1)) or
	($i->{strand} eq 'minus' && ($i->{send} - $j->{sstart} == 1 || $j->{send} - $i->{sstart} == 1)) ) {
	$discordant = 1;
    }
    return($discordant) if $discordant;
    #distance between two mapped anchor on s(subject) should be reasonable
    #either dist1 or dist2 is the inner distance (the gap distance)
    my $dist1 = abs($i->{sstart} - $j->{send});
    my $dist2 = abs($i->{send} - $j->{sstart});
    if($dist1 < $dist2) {
	$discordant = 1 if abs($dist1 - $gapLen) > $GAPLEN_MIS_TOLERANCE;
    } elsif ($dist1 > $dist2) {
	$discordant = 1 if abs($dist2 - $gapLen) > $GAPLEN_MIS_TOLERANCE;
    } else {
	warn "ERROR: Inner and outer distances are the same for the following mappings:\n";
	&dumpMapping($i);
	&dumpMapping($j);
	die "\n";
    }
    return($discordant);
}
sub splitMappingbyAnchor {
    my $mapping = shift;
    my %result;
    for my $i(@$mapping) {
	my $qid = $i->{qid};
	if($result{$qid}) {
	    push @{$result{$qid}},$i;
	} else {
	    $result{$qid} = [$i];
	}
    }
    my $nFound = int(keys %result);
    die "Two anchors expected, found $nFound\n" if $nFound != 2;
    return(values %result);
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
	    $upDownstream, #extension from upstream or downstream?
	);
	$ctg = $i->{sid};
	$name = $id;
	$score = $i->{e};
	($start,$end,$strand,$status,$upDownstream) = &decideGapLocus($targetProfile{$i->{sid}},$i,$gapRecord,$gapLen);
	if(defined $start && defined $end) {
	    push @oneresult,$ctg, $start, $end, $name,$score, $strand, $status;
	    push @results,\@oneresult;
	}
    }
    return(@results);
}
sub decideGapLocus {
    my $ctgInfo = shift;
    my $mapping = shift;
    my $gapRecord = shift;
    my $gapLen = shift;
    my ($start,$end,$strand,$status,$upDownstream);
    &fixMapping($mapping);
    if($mapping->{strand} eq 'plus') {
	#* for N, |----> for mapping
	#case1
	#---*******|--->-------------------------
	#case2
	#----------------------|--->**********---
	$strand = "+";
	if(&isUpstream($mapping,$gapRecord)) {
	    $upDownstream = "upstream";
	    #case2
	    #if the following statement is true, then no extension can be done
	    return if $mapping->{send} >= $mapping->{slen};
	    $start = $mapping->{send}; #-1 is not necessary, this is end for flanking seq
	    if($ctgInfo->{length} >= $mapping->{send} + $gapLen) {
		#gap can be fully filled
		$status = "extension_full";
		$end = $start + $gapLen;
	    } else {
		#gap can be partially filled
		$status = "extension_partial";
		$end = $ctgInfo->{length};
	    }
	} else {
	    $upDownstream = "downstream";
	    #case1
	    #no extension can be done if the left end of anchor is on edge of source contig
	    #mappping is 1-based, output is 0-based start, 1-based end
	    return if $mapping->{sstart} <= 1;
	    $end = $mapping->{sstart}-1; #-1 is necessary
	    if($end - $gapLen >= 0) {
		#gap can be fully filled
		$status = "extension_full";
		$start = $end - $gapLen;
	    } else {
		#gap can be partially filled
		$status = "extension_partial";
		$start = 0;
	    }
	}
	if($start > $end) {
	    &dumpMapping($mapping);
	    die "start > end for + strand!\n";
	}
    } elsif ($mapping->{strand} eq 'minus') {
	#case3
	#---*******<---|-------------------------
	#case4
	#----------------------<---|**********---
	$strand = "-";
	if(&isUpstream($mapping,$gapRecord)) {
	    $upDownstream = "upstream";
	    #case4
	    #if the following statement is true, then no extension can be done
	    return if $mapping->{sstart} >= $mapping->{slen}; #out of bound
	    $end = $mapping->{sstart}; #-1 is not necessary, this is end for flanking seq
	    if($ctgInfo->{length} >= $end + $gapLen) {
		#gap can be fully filled
		$status = "extension_full";
		$start = $end + $gapLen;
	    } else {
		#gap can be partially filled
		$status = "partial";
		$start = $ctgInfo->{length};
	    }
	} else {
	    $upDownstream = "downstream";
	    #case3
	    #no extension can be done if the left end of anchor is on edge of source contig
	    #mappping is 1-based, output is 0-based start, 1-based end
	    return if $mapping->{send} <= 1;# out of bound
	    $start = $mapping->{send}-1; #-1 is necessary
	    if($start - $gapLen >= 0) {
		#gap can be fully filled
		$status = "extension_full";
		$end = $start - $gapLen;
	    } else {
		#gap can be partially filled
		$status = "extension_partial";
		$end = 0;
	    }
	}
	if($start < $end) {
	    &dumpMapping($mapping);
	    die "start < end for - strand!\n";
	}
    } else {
	die "unknown strand:".$mapping->{strand}."\n";
    }
    if($start == $end) {
	&dumpMapping($mapping);
	die "equal start and end !\n";
    }
    return($start,$end,$strand,$status,$upDownstream);
}
sub dumpMapping {
    my $mapping = shift;
    for my $i(keys %{$mapping}) {
	print "key:$i,value: ",$mapping->{$i},"\n";
    }
}
sub fixMapping {
    my $mapping = shift;
    if(abs($mapping->{qend} - $mapping->{qstart}) + 1 == $mapping->{qlen}) {
	#full length mapped, nothing to do
	return;
    } else {
	#qstart<qend anytime
	#sstart>send for - strand
	#sstart<send for + strand
	if($mapping->{strand} eq 'plus') {
	    $mapping->{sstart} = $mapping->{sstart} - ($mapping->{qstart}-1);
	    $mapping->{send}   = $mapping->{send} + ($mapping->{qlen} - $mapping->{qend});
	} else {
	    $mapping->{sstart} = $mapping->{sstart} + ($mapping->{qstart}-0);
	    $mapping->{send}   = $mapping->{send} - ($mapping->{qlen} - $mapping->{qend} - 1);
	}
    }
}
sub isUpstream {
    my $mapping = shift;
    my $gap = shift;
    #chr1:10000-10100
    my ($chr,$start,$end) = &parseQid($mapping->{qid});
    if($end < $gap->[2]) {
	return(1);
    } else {
	return(0);
    }
}
sub parseQid {
    my $id = shift;
    $id =~ /^([^:]+):(\d+)-(\d+)$/ or 
    die "unrecognized query ID: $id\n";
    return($1,$2,$3);
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

    warn "DEBUG:before filtering ",int(@$mapping),"\n" if $debug;
    for my $i(@$mapping) {
	my $pass = 1;
	if($filter->{'idt'}) {
	    $pass = 0 if $i->{idt} < $filter->{idt};
	}
	if($filter->{'qcov'}) {
	    $pass = 0 if $i->{qcov} < $filter->{qcov};
	}
	if($filter->{'e'}) {
	    $pass = 0 if $i->{e} > $filter->{e};
	}
	#discard any mapping that doesn't have equal length on target and query
	#we may fix these mappings in the future
	if(abs($i->{sstart}-$i->{send}) != abs($i->{qstart} - $i->{qend})) {
	    $pass = 0;
	}
	push @result,$i if $pass;
    }
    warn "DEBUG:after filtering ",int(@result),"\n" if $debug;
    return(\@result);
}
sub parseMapping {
#take raw mapping results, and return hash refs
    my $f = shift;
    my @results;
    open IN,'<',$f or die "open($f): $!\n";
    while(<IN>) {
	my @f = split /,/;
	die "expect 12 fields at line $. of $f\n" unless @f==12;
	#s is for subject, the reference seq
	push @results,{
	    qid	=>	$f[0],
	    qlen	=>	$f[1],
	    sid	=>	$f[2],
	    slen	=>	$f[3],
	    qstart	=>	$f[4],
	    qend	=>	$f[5],
	    sstart	=>	$f[6],
	    send	=>	$f[7],
	    strand	=>	$f[8],
	    e		=>	$f[9],
	    idt	=>	$f[10],
	    qcov	=>	$f[11],
	};
    }
    close IN;
    return \@results;
}
sub readBED {
    my $bed=shift;
    warn "NOTICE: reading $bed\n";
    my @out;
    open IN,'<',$bed or die "Cannot read $bed: $!\n";
    while(<IN>) {
	chomp;
	next if /^browser|^#|^track/;
	warn "WARNING: Can't recognize line $.\n" and next unless /^([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)/;
	if ($3>$2) {
	    push @out,[$1,$2,$3,$4];
	} else {
	    push @out,[$1,$3,$2,$4];
	}
    }
    close IN;
    return @out;
}
