package uniline::Stats;
use strict;
use warnings;

my $K = 2; #set tolerance factor to 2 for now, more statistical properties should be explored
#calculate phred-scaled p value
sub calcGapScore {
    my $ret;
    if (@_ == 4) {
	#double anchors
	return(&calcDoubleAnchorScore(@_));
    } elsif (@_ == 3) {
	#single anchor
	return(&calcSingleAnchorScore(@_));
    } else {
	croak("ERROR: 3 or 4 arguments expected\n");
    }
}
sub calcDoubleAnchorScore {
    my ($e1, $e2, $l0, $l1) = @_;
#we need to calculate two probabilities
    #prob of randomly alignment
    #prob of predicting wrong length
    my ($p, $pRandHit1, $pRandHit2, $pWrongL); 
    $pRandHit1 = 1 - exp(-$e1);
    $pRandHit2 = 1 - exp(-$e2);
    if($l0 > 20) {
	#use N(l0,k^2l0^2) for approximation
	if ($l0 == $l1) {
	    $pWrongL = &runRCMD(["2*pnorm(0.5/$K*$l0, lower.tail = TRUE) - 1"]);
	    chomp $pWrongL;
	    $pWrongL =~ s/^\[\d+\]\s+(.*)$/$1/;
	} else {
	    1;
	}
    } else {
	#use Poisson(l0)
	if ($l0 == $l1) {
	    1;
	} else {
	    1;
	}
    }
    $p = -(log($pRandHit1) + log($pRandHit2) + log($pWrongL))/log(10);
    return($p)
}
sub calcSingleAnchorScore {
    my ($e, $l0, $l1) = @_;
    return;
}
#run commands in array ref $cmd
#return output
sub runRCMD {
    my $cmd = shift;
    my $tmp = "/tmp/".rand($$).$$.".R";
    open OUT,">",$tmp or croak("NOTICE: open($tmp) failed, $!\n");
    print OUT join("\n",@$cmd),"\n";
    close OUT;
    return(`Rscript --vanilla $tmp`);
}
1;
