#!/usr/bin/env perl
#take a query fasta, and a bed
#the bee specifies expected mapping location
#run lastz using the expected chr to be mapped as ref, against a corresponding query sequence
#we use an entire chromsome as the reference, rather than a particular region on the chr
#such that we can maximize alignment

use strict;
use warnings;
use SeqMule::Utils;
use SeqMule::Parallel;
use File::Basename qw/basename/;
use File::Spec;

my $tmpdir = "/tmp";
my $desc = "
query_mapping.bed has 6 columns, query_ctg_id, start, end, chr, chr_start, chr_end
chr refers to the chromosome the query_ctg_id is mapped to
";
die "Usage: $0 <query.fa> [query_mapping.bed]\n" unless @ARGV == 2 or @ARGV == 1;
my $chr_dir = "/home/rcf-02/yunfeigu/proj_dir/pacbio_reference/data/hg38_chr/unmasked";
my $chrFile = "/home/rcf-02/yunfeigu/proj_dir/pacbio_reference/data/GCA_000001405.15_GRCh38_full_analysis_set.fna";
my $query = shift @ARGV;
my $bed = shift @ARGV;
my $debug = 1;
my %submitted;
my $qsubLaird = "qsub -V -l walltime=12:59:0 -l nodes=1:ppn=1 -A lc_kw -q laird -l mem=2GB -S /bin/bash";
my $qsubMain = "qsub -V -l walltime=12:59:0 -l nodes=1:ppn=1 -A lc_kw -l mem=2GB -S /bin/bash";
my $cwd = $ENV{PWD};
my $result_dir = File::Spec->catfile($cwd,"lastz_result");
mkdir $result_dir or die "mkdir($result_dir): $!\n" unless -d $result_dir;
if($bed) {
	open BED,'<',$bed or die "open($bed):$!\n";
	while (<BED>) {
		chomp;
		my @f = split;
		die "6 fields expected: $_ at $. of $bed\n" unless @f == 6;
		my ($qID,$qStart,$qEnd,$rID,$rStart,$rEnd) = @f;
		&submitLastzJob({qID=>$qID,rID=>$rID,result_dir=>$result_dir,query=>$query,chr_dir=>$chr_dir,qsub=>$qsubLaird});
	}
	close BED;
} else {
	#when no BED file, align each query contig to every chr, one at a time
	my %allContig = &SeqMule::Utils::readFastaIdx("$query.fai");
	my %allChr = &SeqMule::Utils::readFastaIdx("$chrFile.fai");
	my $count = 0;
	for my $i(keys %allContig) {
		for my $j(keys %allChr) {
			if(&getQCount('laird') <= 500) {
				&submitLastzJob({qID=>$i,rID=>$j,result_dir=>$result_dir,query=>$query,chr_dir=>$chr_dir,qsub=>$qsubLaird});
			} elsif (&getQCount('default')+&getQCount('main')+&getQCount('quick')+&getQCount('main_route') <= 300) {
				&submitLastzJob({qID=>$i,rID=>$j,result_dir=>$result_dir,query=>$query,chr_dir=>$chr_dir,qsub=>$qsubMain});
			} else {
				&blockOnQsub(
						{q=>['default','main','quick','main_route'],n=>300},
						{q=>['laird'],n=> 600}
					);
			}
		}
	}
}
warn "All done\n";

sub getQCount {
	#get number of jobs waiting on specified queue
	my $queueName = shift;
	my $count = `qstat -u yunfeigu|grep $queueName|wc -l`;
	chomp $count;
	warn "Jobs on $queueName: $count\n" if $debug;
	return $count;
}
sub blockOnQsub {
	#sleep until the jobs waiting on specified queues
	#go below certain threshold
	while(1) {
		sleep 1;
		for my $i(@_) {
			my $total = 0;
			map {$total += &getQCount($_)} @{$i->{q}};
			if ($total < $i->{n} ) {
				goto END;
			}
		}
	}
	END: 
	return;
}


=head
REF=/home/rcf-02/yunfeigu/proj_dir/pacbio_reference/data/hg38_chr/unmasked/chr1.fa
QUERY=/auto/rcf-proj/kw/yunfeigu/pacbio_reference/lastz/hx1f2a1_one_ctg_ref_hg38/000000F.fasta
lastz $REF $QUERY --notransition --step=20 --ambiguous=iupac --format=maf --gappedthresh=1000000 --identity=90 > one_ctg_on_hg38chr1.maf
=cut
#####################################################################
sub submitLastzJob {
	my $step = "
	1. mkdir invidual dir
	2. extract the query sequence
	3. generate script
	4. submit";
	my $opt = shift;
	#sleep int(rand(7));
	my ($qID,$rID,$result_dir,$query,$chr_dir,$qsub) = ($opt->{qID},$opt->{rID},$opt->{result_dir},$opt->{query},$opt->{chr_dir},$opt->{qsub});
	my $dir = File::Spec->catdir($result_dir,$qID);
	my $fa = File::Spec->catfile($dir,"query_$qID.fa");
	my $ref = File::Spec->catfile($dir,"ref_$rID.fa");
	my $stde = File::Spec->catfile($dir,"stderr");
	my $stdo = File::Spec->catfile($dir,"stdout");
	my $doneFile = File::Spec->catfile($dir,"$qID.$rID.done");
	mkdir $dir or die "mkdir($dir): $!\n" unless -d $dir;
	!system("samtools faidx $query $qID > $fa") or die "samtools faidx: $!\n";
	symlink File::Spec->catfile($chr_dir,"$rID.fa"),$ref or die "link(): $!\n" unless -e $ref or -l $ref;
	my $script = &SeqMule::Parallel::genTempScript(
		"cd \$TMPDIR",
		"lastz $ref $fa --notransition --step=20 --ambiguous=iupac --format=maf --gappedthresh=1000000 --identity=90 > query_${qID}_ref_${rID}.maf",
		"find . -name '*.maf' -size +4k | xargs -I{} cp {} $dir",
		"touch $doneFile");
	!system("$qsub -e $stde -o $stdo $script") or die "$script: $!\n" unless -e $doneFile or $submitted{$qID.$rID};
	#warn("$qsub -e $stde -o $stdo $script\n") or die "$script: $!\n" unless -e $doneFile or $submitted{$qID.$rID};
	$submitted{$qID.$rID} = 1;
}
