import vcf, logging, os
logging.basicConfig(level=logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')

def generateNonOverlapVCF(vcfFile, outFile = ""):
    '''remove overlaps in sorted vcf
    only keep first record among of a
    group of overlapping records
    '''
    previousRecord = None    
    out = os.path.basename(vcfFile)
    out = out.replace(".vcf",".nonoverlap.vcf")
    if outFile != "":
	out = outFile
    vcfReader = vcf.Reader(open(vcfFile,'r'))
    vcfWriter = vcf.Writer(open(out,'w'),vcfReader)
    overlapCount = 0
    total = 0
    for i in vcfReader:
	total += 1
        if previousRecord != None \
	   and previousRecord.CHROM == i.CHROM \
	   and i.POS < previousRecord.POS + len(previousRecord.REF):
	    overlapCount += 1
            continue
	previousRecord = i		
	vcfWriter.write_record(i)
    logging.info("Total {0} records, {1} overlapping records removed.".format(total, overlapCount))
    logging.info("output writtent to {0}".format(out))
    return
