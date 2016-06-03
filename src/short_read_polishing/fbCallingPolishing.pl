#!/usr/bin/env perl
die "Usage: <fasta> <bam> <prefix>" unless @ARGV == 3;

my $fa = shift;
my $bam = shift;
my $prefix = shift;
my $rawvcf = "$prefix.vcf";
my $normalizedVCF = "$prefix.normalized.vcf";
my $filteredVCF = "$prefix.filteredD5Q20.vcf";
my $nonoverlapVCF = "$prefix.nonoverlap.vcf";
!system("freebayes --fasta-reference $fa --min-coverage 5 $bam > $rawvcf") or die "freebayes: $!";

!system("vt normalize -r $fa $rawvcf > $normalizedVCF") or die "vt: $!";

!system("/home/yunfeiguo/projects/SeqMule/exe/samtools/bcftools/vcfutils.pl varFilter -d 5 $normalizedVCF | perl -ane 'print if /^#/ || \$F[5] ne \".\" && \$F[5] > 20' > $filteredVCF") or die "filter: $!";

!system("/home/yunfeiguo/projects/PacBio_reference_genome/short_read_polishing/rmOverlap.py $filteredVCF $nonoverlapVCF") or die "rmOverlap:$!";

!system("/home/yunfeiguo/projects/PacBio_reference_genome/short_read_polishing/polishContigByIndelSNV.py $fa $nonoverlapVCF") or die "polish:$!";
