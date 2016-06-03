#!/usr/bin/env python
import pysam
import sys
import os
import logging

#logging.basicConfig(level=logging.DEBUG, format = '%(asctime)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')
#min and max distance allowed for a pair
def isSamePair(read1, read2):
    '''determine whether or not two reads come from the same pair'''
    if read1[:-1] != read2[:-1]:
	return False
    return (read1[-1] == '1' and read2[-1] == '2') or (read1[-1] == '2' and read2[-1] == '1')
def isConcordantPair(aln1, aln2):
    '''given two lists of alignments
    determine if any two alignments are concordant
    takes O(nm), n length of aln1, m length of aln2
    '''
    assert(len(aln1) > 0 and len(aln2) > 0)
    for i in aln1:
	for j in aln2:
	    if isConcordantAln(i, j):
		return True
    return False
def isConcordantAln(i, j):
    '''determine if two alignments are concordant
    criteria for concordant alignment:
    same chr, same orientation, within minConcordantDist and maxConcordantDist
    '''
    if i.reference_id != j.reference_id:
	return False
    if i.is_reverse != j.is_reverse:
	return False
    distance = abs(i.get_reference_positions()[0] - j.get_reference_positions()[0])
    if distance > maxConcordantDist or distance < minConcordantDist:
	return False
    return True

def main():
	if len(sys.argv) < 4:
	    print('Usage: {0} <min distance> <max distance> <1.bam 2.bam ...>'.format(os.path.basename(__file__)))
	    raise SyntaxError()
	assert pysam.__version__ == '0.9.0'
	global minConcordantDist, maxConcordantDist
	mapQFilter = 30
	minConcordantDist = int(sys.argv[1])
	maxConcordantDist = int(sys.argv[2])
	for i in xrange(3,len(sys.argv)):
	    bam = sys.argv[i]
	    print(bam)
	    bamHandle = pysam.AlignmentFile(bam, 'rb')
	    totalProperlyAlignedPairs = 0
	    totalMisJoinedPairs = 0
	    prevReadOneName = None #first read in pair
	    prevReadName = None
	    read1Aln = []
	    read2Aln = []
	    totalAlignment = 0
	    unalignedCount = 0
	    lowMapQAlignment = 0
	    totalReads = 0
	    totalPairs = 0
	    readsWithMultipleHits = 0
	    pairsWithMultipleHits = 0
	    currentReadName = None
	    while True:
	#example input	
	##000004F-028-01window1read1	0	chr14	89918901	60	100M	*	0	0	TCTCCTGGATAGCAATAGGGGTATTGTCTGACCTGAGTGACAGGCATACCAGATCTCTGCATGAAATGGAAGGCATTCTGGAATGACACAGGACCAGGGT	*	NM:i:0	MD:Z:100	AS:i:100	XS:i:0
	#000004F-028-01window1read2	0	chr14	90019710	60	100M	*	0	0	GATTTATCTCGAGGTGGTAGATGACAAGTCACCAAAGGAATGGCCCAATCAAGCACAAGTTGTTTTAGGTAGGGGTGACTGGAGGATGAAGTTGGGAAGC	*	NM:i:0	MD:Z:100	AS:i:100	XS:i:0
	#cases to handle
	#assume reads are sorted by name
	#window1read1
	#window1read2
	#window2read1
	#window2read1
	#window2read1
	#window2read2
	#window2read2
	#window3read1
	#window4read1
		try:
		    currentLine = bamHandle.next()
		    totalAlignment += 1
		except StopIteration:
		    break		
		#assume input is sorted by read name
		if currentLine.is_unmapped:
		    unalignedCount += 1
		    continue
		if currentLine.mapping_quality < mapQFilter:
		    lowMapQAlignment += 1
		    continue
		currentReadName = currentLine.query_name
		logging.debug('examining ' + currentReadName)
		if prevReadName != None and prevReadName != currentReadName:
		    totalReads += 1
		prevReadName = currentReadName
		if currentReadName == prevReadOneName:
		    read1Aln.append(currentLine)
		    logging.debug('read1 aln ' + currentReadName)
		    continue
		elif prevReadOneName != None and isSamePair(currentReadName, prevReadOneName):
		    logging.debug('read2 aln ' + currentReadName)
		    read2Aln.append(currentLine)
		else:
		    logging.debug('new pair' + currentReadName)
		    prevReadOneName = currentReadName
		    if len(read1Aln) > 1:
		        readsWithMultipleHits += 1
		    if len(read2Aln) > 1:
		        readsWithMultipleHits += 1
		    if len(read1Aln) > 0 and len(read2Aln) > 0:
			totalPairs += 1
			if len(read1Aln) > 1 and len(read2Aln) > 1:
			    pairsWithMultipleHits += 1
			else:
		            totalProperlyAlignedPairs += 1
			    logging.debug(read1Aln[0].qname + read2Aln[0].qname + ' properly aligned')		    
		            if (not isConcordantPair(read1Aln, read2Aln)):
			        totalMisJoinedPairs += 1
			        logging.debug(read1Aln[0].qname + read2Aln[0].qname + ' not concordant')		    
			    else:
			        logging.debug(read1Aln[0].qname + read2Aln[0].qname + ' concordant')		    
		    else:
			logging.debug('no enough alignments')
			if len(read1Aln) > 0 :
			    logging.debug(read1Aln[0].qname + ' not aligned')
		    read1Aln = []
		    read2Aln = []
		    read1Aln.append(currentLine)

	    if prevReadName != None and prevReadName != currentReadName:
	        totalReads += 1
	    if len(read1Aln) > 0 and len(read2Aln) > 0:
		totalPairs += 1
		if len(read1Aln) > 1 and len(read2Aln) > 1:
		    pairsWithMultipleHits += 1
		else:
	            totalProperlyAlignedPairs += 1
	            logging.debug(read1Aln[0].qname + read2Aln[0].qname + ' aligned')		    
	            if (not isConcordantPair(read1Aln, read2Aln)):
	    	         totalMisJoinedPairs += 1
	            else:
	                 logging.debug(read1Aln[0].qname + read2Aln[0].qname + ' concordant')		    
	    else:
                logging.debug('no enough alignments')
	        if len(read1Aln) > 0 :
	             logging.debug(read1Aln[0].qname + ' not aligned')
	    logging.info('input: {0}'.format(bam))
	    logging.info('{0} ({1}%) out of {2} properly aligned pairs were mis-joined.'.format(totalMisJoinedPairs, 100*totalMisJoinedPairs*1.0/totalProperlyAlignedPairs, totalProperlyAlignedPairs))
	    logging.info('total alignments examined: {0}'.format(totalAlignment))
	    logging.info('unmapped alignments: {0}'.format(unalignedCount))
	    logging.info('alignments failing mapping quality filter (<{0}): {1}'.format(mapQFilter, lowMapQAlignment))
	    logging.info('total reads examined: {0}'.format(totalReads))
	    logging.info('total reads with multiple alignments: {0}'.format(readsWithMultipleHits))
	    logging.info('total read pairs (both reads must map and pass filtering) examined: {0}'.format(totalPairs))
	    logging.info('total read pairs with multiple alignments: {0}'.format(pairsWithMultipleHits))
	
if __name__ == '__main__':
    main()
