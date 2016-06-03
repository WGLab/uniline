#!/usr/bin/env perl
use strict;
use warnings;
use lib "/home/yunfeiguo/projects/SeqMule/lib";
use SeqMule::Utils;
use File::Basename qw/basename/;
use Data::Dumper;
use Getopt::Std;

#param
my $GAPLEN_MIS_TOLERANCE = 50_000_000; #max difference between filled gap and original gap, this is a hard filter
my $GAPLEN_TOLERANCE_FACTOR = 2; #variable k used in P_g calculation
my $nproc = 12; #nproc for bwa
my $filter = {qcov=>95,idt=>95,mapq=>-1}; #mapping filters
my $ERROR_RATE = 0.001; #assume pacbio assembly error rate is 0.001
my $BWA_PARAM = "-A1 -B1 -O1 -E1 -L0"; #special parameters for long read assembly to adjust for high error rate
my $debug = 1;
#external exe
my $fa2size = "/home/yunfeiguo/scripts/seq_related/fa2size.pl";
my $extractGap = "/home/yunfeiguo/scripts/seq_related/extract_gaps.pl";
my $fa2fq = "/home/yunfeiguo/scripts/seq_related/fa2fq.pl";


my %opts;
getopts('sf:',\%opts);
die "Usage: $0 <target assembly> <source assembly>\n".
" -s		skip preprocessing (locating gaps, index db)\n".
" -f <INT>	flank size [default: 100bp]\n"
unless @ARGV == 2;
my $fa = shift @ARGV;
my $db = shift @ARGV;
my $FLANKSIZE = $opts{f} || 100; #length of flank sequences for anchoring
my $MIN_GAP_DIST = 2*$FLANKSIZE; #min distance allowed between two gaps
my $gap = (basename $fa).".gap.bed";
my $genome = (basename $fa).".genome";
my $flankBed = (basename $fa).".gap.flank.bed";
my $flankFa = (basename $fa).".gap.flank.fa";
#my $db = "/home/yunfeiguo/projects/PacBio_reference_genome/falcon_aln/hx1_20150716/2-asm-falcon/hx1f4full.fa";
my $total = 0;
my $unfilled = 0;

warn "filtering threshold: ".Dumper($filter)."\n";
unless($opts{s}) {
    my $bed_for_merge = "/tmp/".rand($$).".tmp.bed.formerge";
    !system("$extractGap $fa > $gap") or die "Gap profiling failed\n" or die "$!\n";
    !system("cat $gap > $bed_for_merge") and !system("bedtools merge -i $bed_for_merge -d $MIN_GAP_DIST > $gap") or die "failed to merge: $!\n";
    #after merging, gap numbers will be gone!
    !system("perl -i -ne '\@f=split;push \@f,\"Gap\$.\";print join(\"\\t\",\@f),\"\\n\"' $gap") or die "failed to add Gap#: $!\n";
    !system("$fa2size $fa > $genome") or die "Genome file creation failed\n";
    #prepare bwa index
    !system("bwa index $db") or die "bwa index: $!\n";
}
if(1) {
    for my $onegapRecord(&readBED($gap)) {
	#next unless $onegapRecord->[0] =~ /chr17_KI270729v1_random/;
	my $onegap = &SeqMule::Utils::genBED([$onegapRecord]);
	my $gapid = join("_",@$onegapRecord);
	my $oneflank = "$onegap.flank.bed";
	my $oneflankfa = "$oneflank.fa";
	my $oneflankfq = "$oneflank.fq";
	my $oneresult = "$oneflank.mapping";
	!system("bedtools flank -i $onegap -g $genome -b $FLANKSIZE > $oneflank") or die "Extraction of flanking sequences failed: $!\n";
	!system("bedtools getfasta -fi $fa -bed $oneflank -fo $oneflankfa") or die "BED2FA failed: $!\n";
	warn "$oneflank\n" if $debug;
	warn "$oneflankfa\n" if $debug;
	!system("$fa2fq $ERROR_RATE $oneflankfa > $oneflankfq") or die "fa2fq.pl failed: $!\n";
	!system("bwa mem $db $oneflankfq -t $nproc > $oneresult") or die "bwa failed: $!\n";
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
    unless(-e $source_fai) {
	!system("samtools faidx $source") or die "samtools faidx: $!\n";
    }
    unless(-e $fai) {
	!system("samtools faidx $db") or die "samtools faidx: $!\n";
    }
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
	    $strand = ($i->{strand} eq '+'? '+':'-');
	    ($start,$end) = &decideDoubleAnchoredGapLocus($i,$j);
	    $score = $i->{mapq} + $j->{mapq} + calcPg($gapLen, abs($end-$start));
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
sub calcPg {
    #calculate phred-scaled P_g
    my $L0 = shift;
    my $Lg = shift;
    my $p;

    #L_g - L_0 ~ N(0, k^2*L_0^2)
    my $sd = $GAPLEN_TOLERANCE_FACTOR*$L0;
    if ($L0 == $Lg) {
	$p = pnorm(0.5, 0, $sd) - pnorm(-0.5, 0, $sd);
    } else {
	$p = pnorm(abs($L0-$Lg), 0, $sd) - pnorm(-abs($L0-$Lg), 0, $sd);
    }
    $p = $p == 0 ? 1/$L0 : $p; #set lower limit to prevent overflow
    return -10*log($p);
    #return $p;
}
sub pnorm {
    #x, mean, variance
    my $x = shift;
    my $mean = shift;
    my $sd = shift;
    my $p = `Rscript -e "print(pnorm($x,$mean,$sd))"`;
    chomp $p;
    if ($p =~ /^\[1\] (\S+)$/) {
	$p = $1;
    } else {
	die;
    }
    return $p;
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
	if ($i->{strand} eq '+') {
	    $start = $j->{send};
	    $end = $i->{sstart}-1;
	} else {
	    $start = $j->{send};
	    $end = $i->{sstart} + 1;
	}
    } elsif ($dist1 > $dist2) {
	if ($i->{strand} eq '+') {
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
    if($mapping1->{strand} eq '+' and $start >= $end) {
	warn "ERROR: start($start) must be smaller than end($end) on plus strand\n";
	&dumpMapping($mapping1);
	&dumpMapping($mapping2);
	die;
    } elsif($mapping1->{strand} eq '-' and $start <= $end) {
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
    if($i->{strand} eq '+') {
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
    if( ($i->{strand} eq '+' && $i->{sstart} <= $j->{send} && $i->{send} >= $j->{sstart} ) or
	($i->{strand} eq '-' && $i->{sstart} >= $j->{send} && $i->{send} <= $j->{sstart}) or
        ($i->{strand} eq '+' && ($i->{sstart} - $j->{send} == 1 || $j->{sstart} - $i->{send} == 1)) or
	($i->{strand} eq '-' && ($i->{send} - $j->{sstart} == 1 || $j->{send} - $i->{sstart} == 1)) ) {
	$discordant = 1;
    }
    return($discordant) if $discordant;
    #distance between two mapped anchor on s(subject) should be reasonable
    #either dist1 or dist2 is the inner distance (the gap distance)
    my $dist1 = abs($i->{sstart} - $j->{send});
    my $dist2 = abs($i->{send} - $j->{sstart});
    if($dist1 < $dist2) {
	#disable hard filtering to promote usage of gap closing score
	##$discordant = 1 if abs($dist1 - $gapLen) > $GAPLEN_MIS_TOLERANCE;
	#$discordant = 1 if abs($dist1 - $gapLen) > $GAPLEN_MIS_TOLERANCE * $gapLen;
    } elsif ($dist1 > $dist2) {
	##$discordant = 1 if abs($dist2 - $gapLen) > $GAPLEN_MIS_TOLERANCE;
	#$discordant = 1 if abs($dist2 - $gapLen) > $GAPLEN_MIS_TOLERANCE * $gapLen;
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
	$score = $i->{mapq};
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
    if($mapping->{strand} eq '+') {
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
    } elsif ($mapping->{strand} eq '-') {
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
		$status = "extension_partial";
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
	if($mapping->{strand} eq '+') {
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
	warn "DEBUG: examining ".Dumper($i)."\n" if $debug;
	my $pass = 1;
	if($filter->{'idt'}) {
	    $pass = 0 if $i->{idt} < $filter->{idt};
	}
	if($filter->{'qcov'}) {
	    $pass = 0 if $i->{qcov} < $filter->{qcov};
	}
	if($filter->{'mapq'}) {
	    $pass = 0 if $i->{mapq} < $filter->{mapq};
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
    #bwa mapping is SAM format
#take raw mapping results, and return hash refs
    my $sam = shift;
    my %ref_length;
    my @results;
    open IN,'<',$sam or die "open($sam): $!\n";
    while(<IN>) {
	if (/^@/) {
	    if (/^\@SQ\tSN:(\S+)\tLN:(\d+)/) {
		##@SQ	SN:000000F:1-19000000	LN:19000000
		$ref_length{$1} = $2;
	    }
	    next;
	}
	my @f = split /\t/;
	die "expect >=11 fields at line $. of $sam\n" unless @f >= 11;
	#s is for subject, the reference seq
#E00170:50:H0D6EALXX:8:1101:16001:69819	163	000000F:1-19000000	1	60	97S54M	=	72	220	CAAATGGCCTCCATAAGCAACTCACGCTGGTTTAACTCAGGTTAGTTTACACCATTGGACCTTCCCAAATAGTACATTGTTTTTTTGAGATAAAGCATTGAGAGAGCTCTCCTTAACGTGACACAATGGAAGGACTGGAACACATACCCAC	AAFFFKKKKKKKKKKKKKKKFKKKKFKKKKKKKFKKFAAKKKKKKKKKKKKKKKKKKKKKKKKKFAKFAFKKKKAFFKKKKKKKKKKK77<FF<K7F7,AKFAKAFKKF7AA<AFF7AAFFFKKKKK7AFA7FKF7A<FA<AKKFA<A7<<	NM:i:0	MD:Z:54	AS:i:54	XS:i:0
	#QNAME	FLAG	RNAME	POS	MAPQ	CIGAR	RNEXT	PNEXT	TLEN	SEQ	QUAL			  
	#0	1	2	3	4	5	6	7	8	9	10
	next if $f[5] eq '*';
	my $mapping_detail = &parseCigar($f[5]);
	my $strand = $f[8] & 16? '-':'+'; #flag 16 for read reverse strand
	push @results,{
	    qid		=>	$f[0],
	    qlen	=>	length($f[9]),
	    sid		=>	$f[2],
	    slen	=>	$ref_length{$f[2]},
	    qstart	=>	$mapping_detail->{qstart}, 
	    qend	=>	$mapping_detail->{qend}, 
	    sstart	=>	$f[3] + ($strand eq '+' ? $mapping_detail->{sstart} : $mapping_detail->{send}) - 1, #sstart,send are 1-based
	    send	=>	$f[3] + ($strand eq '+' ? $mapping_detail->{send} : $mapping_detail->{sstart}) - 1,
	    strand	=>	$strand,
	    mapq	=>	$f[4],
	    idt		=>	$mapping_detail->{idt},
	    qcov	=>	$mapping_detail->{qcov},
	};
    }
    close IN;
    return \@results;
}
sub parseCigar {
    #return 1-based coordinate
    my $cigar = shift;
    my $mapping_detail = {};
    warn "DEBUG: parsing $cigar\n" if $debug;
    #cigar characters: M, I, D, S, *, =, H
    #we use the first M as the qstart, last M as qend
    #
    #how to handle I?
    #reference: AT*CG
    #read:      ATTCG
    #qstart, qend will count both M and I
    #start and send will only count M
    #
    #how to handle D?
    #reference: ATCG
    #read:      A*CG
    #qstart, qend will count M
    #sstart and send will count both M and D
    #
    #how to handle H?
    #reference: ATCG
    #read:      -TCG
    #the impact for ref and read is the same
    #
    #how to handle S?
    #count as if it is M unless it occurs at the beginning or end

    #issues with strandness
    #positive strand, nothing special
    #negative strand, the read is converted to reverse complement, and mapped
    #therefore the cigar string is generated relative to the reverse complement
    my @cigar_fields = &cigar2array($cigar);
    $mapping_detail->{qstart} = -1;
    $mapping_detail->{qend} = -1;
    $mapping_detail->{sstart} = -1;
    $mapping_detail->{send} = -1;
    my ($offsetD, $offsetI) = (0, 0); #count D and I for ref and read, respectively
    my $countM = 0;
    my $totalBP = 0;
    my $i = 0; #index for cigar fields
    while ($i < @cigar_fields) {
	if ($cigar_fields[$i] eq 'M' && $mapping_detail->{qstart} == -1) {
	    $mapping_detail->{qstart} = $totalBP - $offsetD + 1;
	    $mapping_detail->{sstart} = $totalBP - $offsetI + 1; #this is relative position to POS field
	    $mapping_detail->{qend} = $totalBP - $offsetD + $cigar_fields[$i + 1];
	    $mapping_detail->{send} = $totalBP - $offsetI + $cigar_fields[$i + 1];
	    $countM += $cigar_fields[$i+1];
	} elsif ($cigar_fields[$i] eq 'M') {
	    #keep updating this until we see the last M
	    $mapping_detail->{qend} = $totalBP - $offsetD + $cigar_fields[$i + 1];
	    #print("total:$totalBP,offsetD:$offsetD,f:".$cigar_fields[$i+1]."\n");
	    $mapping_detail->{send} = $totalBP - $offsetI + $cigar_fields[$i + 1];
	    $countM += $cigar_fields[$i+1];
	} elsif ($cigar_fields[$i] eq 'H') {
	    1;
	} elsif ($cigar_fields[$i] eq 'S') {
	    1;
	} elsif ($cigar_fields[$i] eq 'D') {
	    $offsetD += $cigar_fields[$i + 1];
	} elsif ($cigar_fields[$i] eq 'I') {
	    $offsetI += $cigar_fields[$i + 1];
	}
	$totalBP += $cigar_fields[$i + 1];
	$i += 2;
    }
    $mapping_detail->{idt} = $countM/($totalBP - $offsetD) * 100; #identity percentage
    $mapping_detail->{qcov} = 100 * ($mapping_detail->{qend} - $mapping_detail->{qstart} + 1) / ($totalBP - $offsetD); #percentage of query covered
    return $mapping_detail;
}
sub cigar2array {
    #convert 1M2D3M4I to
    #(M,1,D,2,M,3,I,4)
    my $cigar = shift;
    my @cigar_fields = $cigar =~ /(\d+)(\D+)/g;
    my @return;
    for my $i(0..(@cigar_fields/2 - 1)) {
	push @return,$cigar_fields[$i*2+1],$cigar_fields[$i*2];
    }
    return @return;
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
