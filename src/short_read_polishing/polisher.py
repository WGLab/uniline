import vcf, logging
from Bio import SeqIO
from sets import Set
logging.basicConfig(level=logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')
class polisher:
    def __init__(self, fasta, vcf, het = 0.1):
	self.vcf = vcf
	self.fasta = fasta
	#variants below this threshold are considered homozygous
	self.heterozygosityThreshold = het
	self.varCount = 0
	self.varUsedCount = 0
	self.baseCount = 0
	self.baseCorrectedCount = 0
	self.baseToCorrectCount = 0
	self.SNVcount = 0
	self.lt5IndelCount = 0
	self.gt5IndelCount = 0
    def isVCFValid(self):
	'''make sure VCF records are sorted
	in ascending order, and have no
	overlap	within each contig
	'''
	previousRecord = None
	vcfReader = vcf.Reader(open(self.vcf,'r'))
	count = 0
	for i in vcfReader:
	    count += 1
	    if count % 100000 == 0:
		logging.info("isVCFValid(): {0} variants scanned...\n".format(count))
            if previousRecord != None:
	        if previousRecord.CHROM == i.CHROM:
		    if i.POS < previousRecord.POS + len(previousRecord.REF):
			print(i.CHROM)
			print(i.POS)
			print(previousRecord.POS)
			return False
	    previousRecord = i		  
	return True

    def outputPolished(self, file, validation = True):
	'''output polished seq to file
	using buffered writer to
	reduce number of write requests
	'''
	if validation:
	    logging.info("Starting validating VCF")
	    assert(self.isVCFValid())
	logging.info("Starting output polished FASTA")
	fh = open(file, 'w', 10*1024*1024) #buffer size: 10MB
	#WARNING: not completely compatible with VCF4.2
	vcfReader = vcf.Reader(open(self.vcf,'r'))
	fastaDict = SeqIO.index(self.fasta, 'fasta')
	visitedChrom = Set()
	try:
	    while True:
	        currentVCFRecord = next(vcfReader)
	        self.varCount += 1
	        if currentVCFRecord.heterozygosity <= self.heterozygosityThreshold:
		    self.varUsedCount += 1
	            break
	except StopIteration:
	    fh.close()
	    logging.info("No VCF record.")
	    return
	#make sure to handle the following corner cases
	#output all contigs regardless whether they have variants or not
	#finish the current contig after the last variant
	currentFastaPos = 0
	#the actual sequence will be loaded into this variable in the form of SeqRecord obj
	currentChrom = fastaDict[currentVCFRecord.CHROM]
        fh.write(">{0} {1}\n".format(currentChrom.id, currentChrom.description))
	while True:
	    assert(currentFastaPos < len(currentChrom) and currentVCFRecord.POS <= len(currentChrom))
	    fh.write(str(currentChrom.seq[currentFastaPos:currentVCFRecord.POS - 1]))
	    fh.write(str(currentVCFRecord.ALT[0]))
	    self.baseCorrectedCount += len(str(currentVCFRecord.ALT[0]))
	    self.baseToCorrectCount += len(str(currentVCFRecord.REF))
	    currentFastaPos = currentVCFRecord.POS + len(currentVCFRecord.REF) - 1
	    if len(currentVCFRecord.REF) == 1 and len(currentVCFRecord.ALT[0]) == 1:
		self.SNVcount += 1
	    elif abs(len(currentVCFRecord.REF) - len(currentVCFRecord.ALT[0])) < 5:
		self.lt5IndelCount += 1
	    else:
		self.gt5IndelCount += 1
	    try:
		while True:
	            currentVCFRecord = next(vcfReader)
	            self.varCount += 1
		    if self.varCount % 100000 == 0:
			logging.info("outputPolished(): {0} variants scanned...\n".format(self.varCount))
		    if currentVCFRecord.heterozygosity <= self.heterozygosityThreshold:
		        self.varUsedCount += 1
			break
	    except StopIteration:
	        break
	    if currentChrom.id != currentVCFRecord.CHROM:
		#we are about to switch contig
		#output the remaining sequences
		fh.write(str(currentChrom.seq[currentFastaPos:len(currentChrom)]) + '\n')
	        self.baseCount += len(currentChrom)
		visitedChrom.add(currentChrom.id)
		#switch contig
		currentChrom = fastaDict[currentVCFRecord.CHROM]
		currentFastaPos = 0
		fh.write(">{0} {1}\n".format(currentChrom.id, currentChrom.description))
	if currentFastaPos < len(currentChrom):
	    #this is the last contig that we haven't finished
	    fh.write(str(currentChrom.seq[currentFastaPos:len(currentChrom)]) + '\n')
	    self.baseCount += len(currentChrom)
	    visitedChrom.add(currentChrom.id)
	#output unvisited chromosomes
	for i in fastaDict:
	    if not i in visitedChrom:
		currentChrom = fastaDict[i]
		fh.write(">{0} {1}\n".format(currentChrom.id, currentChrom.description))
		fh.write(str(currentChrom.seq) + '\n')
	        self.baseCount += len(currentChrom)
	fh.close()
	logging.info("FASTA:{0}\nVCF:{1}".format(self.fasta, self.vcf))
	logging.info("{0} variants scanned, {1} passing heterozygosity threshold ({2}) are considered homozygous.".format(self.varCount, self.varUsedCount, self.heterozygosityThreshold))
	logging.info("{0} base pairs in source FASTA ({1}).".format(self.baseCount, self.fasta))
	logging.info("{0} base pairs were corrected (replaced/deleted/inserted) to {1} base pairs.".format(self.baseToCorrectCount, self.baseCorrectedCount))
	logging.info("{0} SNVs were corrected, {1} indels (<5bp) were corrected, {2} indels (>=5bp) were corrected".format(self.SNVcount, self.lt5IndelCount, self.gt5IndelCount))
	logging.info("{0} base pairs in polished sequences.".format(self.baseCount + self.baseCorrectedCount - self.baseToCorrectCount))
	logging.info("Polished FASTA written to {0}".format(file))
	return
def main():
    #unit test
    #test = polisher('tiny_example/original_seq.fa', 'tiny_example/indel_in_middle.bam.vcf')
    #assert(test.isVCFValid())
    #test2 = polisher('tiny_example/original_seq.fa', 'tiny_example/indel_in_middle.wrongorder.vcf')
    #assert(test2.isVCFValid() == False)
    #test3 = polisher('tiny_example/original_seq.fa', 'tiny_example/indel_in_middle.overlap.vcf')
    #assert(test3.isVCFValid() == False)
    #test3 = polisher('tiny_example/indel_in_middle.fa', 'tiny_example/indel_in_middle.bam.vcf')
    #test3.outputPolished('tiny_example/indel_in_middle.fixed.fa')
    test4 = polisher('example/000000F_left.fixed.fa','example/one_ilmn_fq_on_000000FleftFixed.filteredD5Q20.vcf')
    test4.outputPolished('/tmp/y.fa')
    #test5 = polisher('/home/yunfeiguo/projects/PacBio_reference_genome/falcon_aln/hx1_20150716/2-asm-falcon/hx1f4full.fa','/home/yunfeiguo/projects/PacBio_reference_genome/short_read_polishing/firstRoundPolish/samtoolsFirstRound.filteredD5Q20.valid.NonOverlap.vcf.gz')
    #test5.outputPolished('/home/yunfeiguo/projects/PacBio_reference_genome/short_read_polishing/firstRoundPolish/hx1f4full.fixed.fa')
    print("all tests passed")
if __name__ == '__main__':
    main()
