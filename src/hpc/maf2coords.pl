#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
#take a few MAF files
#convert them into coords format
#used by MUMmer
#assume that only two sequences
#are compared in each alignment

my $IDT_DEFAULT = 100;
my %opts;
getopts('f',\%opts);
die "Usage: $0 <1.maf 2.maf ...>\n-f	treat argument as list of files\n" unless @ARGV >= 1;
if($opts{f}) {
    my @allfiles = `cat $ARGV[0]`;
    chomp @allfiles;
    @ARGV = @allfiles;
}
for my $maf(@ARGV) {
    &maf2coords($maf);
}
warn "All done\n";
#####################SUBROUTINE######################
sub maf2coords {
    #take a maf
    #parse it and print in coords
    my $maf = shift;
    my $aln = [];
    open IN,'<',$maf or die "open($maf): $!\n";
    while(<IN>) {
        #example of MAF
        #comments are marked by #, alignments are separated by blank lines
###maf version=1 scoring=lastz.v1.03.73
## lastz.v1.03.73 --notransition --gap=1000,1 --step=20 --ambiguous=iupac --format=maf --gappedthresh=1000000 --identity=80 --progress=10 --maxwordcount=1 --masking=0 
#a score=1100772
#s chr1           146285887 12312 + 248956422 GATCTCGGCTCACTGCAAGCTC...
#s 001072F-004-01       466 12252 +     12718 GATCTCGGCTCACTGCAAGCTC...
#type ID          start     size  strand totalLen sequence
        next if /^#/;
        chomp;
        if(/^\s*$/) {
            if (@$aln) {
                print join("\t",&parseMAF($aln)),"\n";
                $aln = [];
            } else {
                1;
            }
        } else {
            my @f = split;
            if ($f[0] =~ /^s/) {
                die "ERROR: 7 fields expected at line $. of $maf\n"
                unless @f == 7;
                push @$aln,\@f;
            } else {
                #not actual alignment
                #skip for now
                1;
            }
        }
    }
    if (@$aln) {
        print join("\t",&parseMAF($aln)),"\n";
        $aln = [];
    }

}
sub parseMAF {
    #assume coords is 1-based, maf is 0-based
    #take 2 elements, [type,ID,start,size,strand,totalLen,sequence]
    #                   0   1   2     3       4   5        6
    #1st is ref ID in alignment
    #2nd is query ID in alignment
    #output in coords format
    #coords example
    #refs     refend          qs     qend        reflen   qlen     idt       refTotalLen      qTotalLen          refcov   qcov    refID     qID
    #16627    17364  |        1      738  |      738      738  |   100.00  | 248956422   561012  |     0.00     0.13  | chr1      KE141541.1
    #16865    17372  |        1      508  |      508      508  |   100.00  | 248956422      508  |     0.00   100.00  | chr1      ADDF02338563.1

    my $aln = shift;
    die "Expect 2 elements in \$aln\n" unless @$aln == 2;
    my ($ref,$query) = @$aln;
    my ($refs,$refend,$qs,$qend,$reflen,$qlen,$idt,$refTotalLen,$qTotalLen,$refcov,$qcov,$refID,$qID);
    #first calculate numbers unrelated to strandness
    $reflen = $ref->[3];
    $qlen = $query->[3];
    $refTotalLen = $ref->[5];
    $qTotalLen = $query->[5];
    $qcov = $qTotalLen==0? 0:$qlen/$qTotalLen;
    $refcov = $refTotalLen==0? 0:$reflen/$refTotalLen;
    $refID = $ref->[1];
    $qID = $query->[1];
    $idt = $IDT_DEFAULT;

    #now calculate numbers related to strand
    if($ref->[4] eq '+') {
        $refs = $ref->[2] + 1; #0 to 1 based
        $qs = $query->[2] + 1;
        $refend = $refs + $reflen - 1;
        $qend = $qs + $qlen - 1;
    } elsif ($ref->[4] eq '-') {
        $refs = $ref->[2] - 1; #0 to 1 based
        $qs = $query->[2] - 1;
        $refend = $refs - $reflen + 1;
        $qend = $qs - $qlen + 1;
    } else {
        die "+ or - for strand only: ",$ref->[4],"\n";
    }

    return($refs,$refend,$qs,$qend,$reflen,$qlen,$idt,$refTotalLen,$qTotalLen,$refcov,$qcov,$refID,$qID);
}
