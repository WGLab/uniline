#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/basename/;
use File::Spec;
use Getopt::Std;

my %options;
getopts("nl",\%options);
die "Usage: $0 <abs path to reference chr dir> <abs path to query fasta>\n".
" -n		no sleep\n".
" -l		run in local\n"
unless @ARGV>=2;
#my $chr_dir = "/home/yunfeiguo/database/hg_index/GRCh38_patch/GRCh38p4_chr";
my $chr_dir = shift @ARGV;
my $tmpdir = "/tmp";
my $pwd = $ENV{PWD};
#my $query = "/home/yunfeiguo/projects/PacBio_reference_genome/falcon_aln/hx1_20150716/2-asm-falcon/p_ctg.fa";
my $query = shift @ARGV;
my $run_dir = "run";
my $queryPrefix= "000000F";
mkdir $run_dir unless -d $run_dir;
chdir $run_dir;
$pwd = File::Spec->catfile($pwd,$run_dir);
my $qsub_option = <<SCRIPT;
#!/bin/bash
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
#\$ -l h_vmem=8G
set -e
SCRIPT
#\$ -m ea
#\$ -M guoyunfei1989\@gmail.com
#
warn "ref chr dir: $chr_dir\n";
warn "query fasta: $query\n";

for my $i(glob File::Spec->catfile($chr_dir,"*.fa")) {
    #my ($prefix) = (basename $i)=~/^(chr[\dXYM]{1,2})\.fa$/ or warn "ERROR: failed to match $i\n" and next;
    my ($prefix) = (basename $i)=~/^(.*?)\.fa$/ or warn "ERROR: failed to match $i\n" and next;
    my $runtimetmpdir = rand($$)."mummer";
    $prefix = "${queryPrefix}_$prefix";
    $prefix =~ s/[\:\*]/_/g;
    #after prefix for script and subfolder is done, process file name to make it unambiguous
    $i =~ s/([\:\*])/\\$1/g;
    my $script = File::Spec->catfile($tmpdir,"runmummer.$prefix.".rand($$).".sh");
    mkdir($prefix) or die "mkdir($prefix): $!\n" unless -d $prefix;
    open OUT,">",$script or die "ERROR: failed to write to $script: $!\n";
    print OUT $qsub_option,"\n";
    #redirect stderr and stdout
    print OUT "#\$ -e ".File::Spec->catfile($prefix,"stderr")."\n";
    print OUT "#\$ -o ".File::Spec->catfile($prefix,"stdout")."\n";
    #use file staging to reduce IO
    print OUT "TMP=$runtimetmpdir\n";
    print OUT "mkdir \$TMP\n";
    print OUT "cd \$TMP\n";
    print OUT "nucmer -c 400 -l 150 --prefix=$prefix $i $query\n";
    print OUT "show-coords -r -c -l -k $prefix.delta > $prefix.coords\n";
    print OUT "show-tiling $prefix.delta > $prefix.tiling\n";
    print OUT "if [ -s $prefix.tiling ]; then mummerplot --postscript $prefix.tiling -p $prefix;fi\n";
    #print OUT "delta-filter -m $prefix.delta > $prefix.filtered.delta\n";
    #print OUT "show-coords -r -c -l $prefix.filtered.delta > $prefix.filtered.coords\n";
    #print OUT "show-tiling $prefix.filtered.delta > $prefix.filtered.tiling\n";
    #print OUT "if [ -s $prefix.filtered.tiling ]; then mummerplot --postscript $prefix.filtered.tiling -p $prefix.filtered;fi\n";
    print OUT "cp -r * ".File::Spec->catfile($pwd,$prefix)." \n";
    print OUT "touch ".File::Spec->catfile($pwd,$prefix,"$prefix.done")."\n";
    print OUT "rm -rf \$TMP\n";
    close OUT;
    unless (-e File::Spec->catfile($pwd,$prefix,"$prefix.done")) {
	!system("chmod +x $script\n") or die "chmod +x $script: $!\n" ;
	if ($options{'l'}) {
	    !system("$script &") or die "$script: $!\n" ;
	} else {
	    !system("qsub $script\n") or die "$script: $!\n" ;
	}
	sleep 5 unless $options{'n'} ;
	warn("qsub $script\n");
    }
}
warn "Results written to $run_dir\n";
warn "All done\n";

#process coords file
#perl -ne 'next if /^[=\/\s]|^NUCMER/;s/\|//g;print' hx1_50kb_on_hg38.filtered.coords > hx1_50kb_on_hg38.filtered.whitespace.coords
