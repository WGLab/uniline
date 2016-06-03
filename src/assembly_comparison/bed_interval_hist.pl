#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/basename/;
use Getopt::Std;

my %opts;
getopts('f:',\%opts);
die "Usage: $0 <bed>\n".
" -f <INT>	min interval length. Default: no filtering.\n" unless @ARGV == 1;
my $bed = shift @ARGV;
my @intervalLen = &getBedIntervalLen($bed);
my $prefix = basename $bed;
my $pdf = "$prefix.pdf";
my $interval = "$prefix.interval";
if(defined $opts{f}) {
    warn "Filtering intervals smaller than $opts{f}\n";
    @intervalLen = grep {$_>=$opts{f}} @intervalLen;
}
die "failed to generate intervals: $!\n" unless &writeInterval($interval,@intervalLen);
die "Failed to plot: $!\n" unless &drawInterval($pdf,$interval);
warn "Plot $pdf done.\n";


#############################################
sub getBedIntervalLen{
    warn "Reading intervals\n";
    my $bed = shift;
    my @return;
    open IN,"<",$bed or die "open($bed): $!\n";
    while(<IN>) {
	my @f=split;
	die "ERROR: 3 fields expected at $. of $bed\n" unless @f>=3;
	push @return,($f[2]-$f[1]);
    }
    close IN;
    return(@return);
}
sub writeInterval {
    warn "Writting intervals\n";
    my $file = shift;
    my @intervalLen = @_;
    open OUT,">",$file or die "open($file): $!\n";
    print OUT join("\n",@intervalLen);
    close OUT;
    return(1);
}
sub drawInterval {
    warn "Plotting histogram...\n";
    my $rScript = "/tmp/".rand($$).".R";
    my $rSrc = "
require('ggplot2')
args <- commandArgs(TRUE)
input <- as.character(args[1])
output <- as.character(args[2])
data = read.table(input)
data\$V1 = as.numeric(data\$V1)
pdf(output)
print(ggplot(data,aes(x=data\$V1))+geom_histogram()+xlab('Interval length (bp)'))
######
#the following function superimposes a poisson distribution to the histogram
#print(ggplot(data.exp2,aes(x=data\$V1)) + geom_histogram() + xlab('filtered read length')+stat_function(geom='line',fun = function(...,lambda,total){dpois(...,lambda)*total}, args = list(lambda = mean(data.exp2\$data.exp), total = dim(data.exp2)[1]),colour = 'red',fill=NA)
dev.off()
";
open OUT,">",$rScript or die "open($rScript): $!\n";
print OUT $rSrc;
close OUT;
!system("Rscript --vanilla $rScript $interval $pdf") or die "Rscript($rScript): $!\n";
return(1);
}
