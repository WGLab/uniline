#!/usr/bin/env python
import sys, argparse, os
sys.path.append('/home/yunfeiguo/projects/PacBio_reference_genome/short_read_polishing')
from polisher import polisher

def main():
    parser = argparse.ArgumentParser(description='Polish draft assembly with variant calling results')
    #by default, argparse automatically adds -h, --help 
    parser.add_argument('fasta', metavar = 'FASTA', type=str, 
	    		help='FASTA file to be polished')
    parser.add_argument('vcf', metavar = 'VCF', type=str, 
	    		help='Sorted, non-overlapping, quality-filtered VCF')
    parser.add_argument('--het', metavar = 'threshold', nargs = '?', type=float, default=0.1,	    
	    help='Heterozygosity threshold, below this is considered as homozygous. Default: 0.1')
    parser.add_argument('--no-validate', dest='validate',action = 'store_false',
	    		help='skip VCF validation')
    parser.add_argument('--version', action = 'version', version='alpha')
    args = parser.parse_args()

    p = polisher(args.fasta, args.vcf, args.het)
    p.outputPolished(os.path.basename(args.fasta) + '.fixed', args.validate)
if __name__ == '__main__':
    main()
