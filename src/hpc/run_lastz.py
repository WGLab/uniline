#!/usr/bin/env python
'''take a query fasta, and optionally a bed
run lastz with hg38 as reference
if a bed is given, then each sequence in query fasta should be mapped a specific reference sequence
'''
import os
import sys
import logging
import tempfile
import time
logging.basicConfig(level=logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')
#global variables
cleanQ = []
submitInterval = 1 #number of seconds between qstat
checkInterval = 1 #sec
debug = False
maxConcurrentLaird = 100
maxConcurrentBcl = 500
maxConcurrentOther = 100
jobReg = {} #record job queue, status, sentinel file path
#chr_dir = "/auto/rcf-proj/kw/yunfeigu/pacbio_reference/lastz/ref_hg38+decoy_query_hx1f4full_no_hg38+decoy.mummer/GRCh38_full_analysis_set_plus_decoy_hla.noMetaChar_chr"
#chrFile = "/auto/rcf-proj/kw/yunfeigu/pacbio_reference/lastz/ref_hg38+decoy_query_hx1f4full_no_hg38+decoy.mummer/GRCh38_full_analysis_set_plus_decoy_hla.noMetaChar.fa"
chr_dir = "/home/yunfeiguo/projects/PacBio_reference_genome/lastz/ref_novel_seq_absent_in_yh/chroms"
chrFile = "/home/yunfeiguo/projects/PacBio_reference_genome/lastz/ref_novel_seq_absent_in_yh/novel_by_ml_absent_in_yh.nometa.fasta"
#chrFile = "/home/yunfeiguo/database/hg_index/GRCh38_analysis/GRCh38_full_analysis_set_plus_decoy_hla.noMetaChar.fa"
#chr_dir = "/home/yunfeiguo/database/hg_index/GRCh38_analysis/GRCh38_full_analysis_set_plus_decoy_hla.noMetaChar_chr"
qsubLaird = "qsub -V -l walltime=199:59:0 -l nodes=1:ppn=1 -A lc_kw -q laird -l mem=2GB -S /bin/bash"
qsubMain = "qsub -V -l walltime=23:59:0 -l nodes=1:ppn=1 -A lc_kw -l mem=2GB -S /bin/bash"
qsubBcl = "qsub -V -cwd -l h_vmem=7G -S /bin/bash"
#prototype of lastz options
lastzOpt = "--notransition --gap=1000,1 --step=20 --ambiguous=iupac --format=maf --gappedthresh=10000 --identity=80 --progress=10 --maxwordcount=1 --masking=0 ";
cwd = os.getcwd()
#noSleep = False
noSleep = True

def safeMkdir(dir):
    if not os.path.isdir(dir):
        try:
            os.mkdir(dir)
        except OSError:
            logging.critical("%s exists" % dir)
def readFastaIdx(idx):
    fh = None
    try:
       fh = open(idx,'r')
    except IOError as e:
       logging.critical("I/O error: {0}".format(idx))
       raise
    ctg = {}
    for line in fh:
        f = line.split()
        assert len(f) == 5
        (id,length,offset,nCharLine,nByteLine) = f
        ctg[id] = {
                'length':length,
                'offset':offset,
                'nCharLine':nCharLine,
                'nByteLine':nByteLine,
                }
    fh.close()
    return(ctg)
def genTempScript(prefix,suffix,cmd):
    assert type(cmd) is list
    fd, tmpFile = tempfile.mkstemp(prefix = 'ref'+prefix, suffix='pairwiseLastz'+str(suffix))
    os.write(fd,'#!/bin/bash\n')
    os.write(fd,'set -e\n') #let shell run the script, exit at first error
    os.write(fd,'set -o pipefail\n') #let shell run the script, exit at first error
    for c in cmd:
        os.write(fd,c+'\n')
    os.close(fd)
    cleanQ.append(tmpFile) #clean at the very end
    return(tmpFile)
def clean():
    for f in cleanQ:
        os.remove(f)
def getJobCount(qName):
    '''
    check all jobs not finished
    remove any if finished
    return count for a particular queue
    '''
    #when all jobs are submitted, wait before check
    if not noSleep:
        if len(jobReg) == len([v for k,v in jobReg.items() if v['submitted']]):
            time.sleep(checkInterval)
            #pass
        else:
            time.sleep(submitInterval)
            #pass
    #make sure we mark a job as finished once it is done
    for i in jobReg:
        if jobReg[i]['finished']:
            pass
        else:
            if debug:
                print("checking %s" % jobReg[i]['doneFile'])
            if os.path.exists(jobReg[i]['doneFile']):
                if debug:
                    print("%s exists"%jobReg[i]['doneFile'])
                    print("%s is done"%i)
                jobReg[i]['finished'] = True
    total = len([v for k,v in jobReg.items() if v['submitted'] and v['q']==qName])
    finished = len([v for k,v in jobReg.items() if v['submitted'] and v['finished'] and v['q']==qName])
    if debug:
        print("total:%d,finished:%d"%(total,finished))
    return(total-finished)

def submitLastzJob(result_dir, fa, qCount, rID, q, qsub):
    '''
    example lastz command
    REF=/home/rcf-02/yunfeigu/proj_dir/pacbio_reference/data/hg38_chr/unmasked/chr1.fa
    QUERY=/auto/rcf-proj/kw/yunfeigu/pacbio_reference/lastz/hx1f2a1_one_ctg_ref_hg38/000000F.fasta
    lastz $REF $QUERY --notransition --step=20 --ambiguous=iupac --format=maf --gappedthresh=10000 --identity=90 > one_ctg_on_hg38chr1.maf
    
    steps to submit
    	1. mkdir invidual dir
    	2. extract the query sequence
    	3. generate script
    	4. submit"
    '''
    if not noSleep:
        time.sleep(submitInterval)
    dir = os.path.join(result_dir,rID)
    ref = os.path.join(dir,"ref_%s.fa" % rID)
    stde = os.path.join(dir,"stderr"+str(qCount))
    stdo = os.path.join(dir,"stdout"+str(qCount))
    jobID = rID + str(qCount)
    if debug:
        print(jobID)
    individualDoneFile = os.path.join(dir,"%s_query%d.done" % (rID,qCount))
    safeMkdir(dir)
    jobReg[jobID] = {
            'doneFile':individualDoneFile,
            'q':q,
            'submitted':False,
            'finished':False,
            }
    script = None
    if os.path.lexists(ref):
        os.remove(ref)
    #create the symlink for reference regardless it existence
    os.symlink(os.path.join(chr_dir,"%s.fa"%rID),ref)
    script = genTempScript(rID, qCount, [
        "cd $TMPDIR",
		"lastz %s %s %s > ref_%s_vs_query%d.maf" % (ref,fa,lastzOpt,rID,qCount),
		"find . -name '*.maf' -size +4k | xargs -I{} cp {} %s" % (dir),
		"touch %s" % (individualDoneFile),
        ])
    if script is None:
        logging.critical("Nonetype for script")
    if debug:
        print("done file exists?")
        print(os.path.exists(individualDoneFile))
        print("submitted?")
        print(jobReg[jobID]['submitted'])
    if os.path.exists(individualDoneFile) or jobReg[jobID]['submitted']: 
        pass
    else:
        os.system("%s -e %s -o %s %s" % (qsub, stde, stdo, script))
        #print("%s -e %s -o %s %s" % (qsub, stde, stdo, script))
        jobReg[jobID]['submitted'] = True
#####################################################################
desc = '''
query_mapping.bed has 6 columns, query_ctg_id, start, end, chr, chr_start, chr_end
chr refers to the chromosome the query_ctg_id is mapped to'''
if len(sys.argv) == 0 or len(sys.argv) == 1:
    sys.exit("Usage: %s <query1.fa query2.fa ...>" % sys.argv[0])
query = sys.argv[1:len(sys.argv)]
bed = None
result_dir = os.path.join(cwd,"lastz_result")
safeMkdir(result_dir)
	#align each query contig to every chr, one at a time
#allChr = ["alt","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"]
allChr = readFastaIdx(chrFile+".fai")
for j in range(0,len(query)):
    oneQuery = os.path.abspath(query[j])
    print(oneQuery)
    for i in allChr:
        #if getJobCount('laird') < maxConcurrentLaird:
        #    submitLastzJob(result_dir = result_dir, fa = oneQuery, qCount = j, rID = i, q = 'laird', qsub = qsubLaird)
        #elif getJobCount('other') < maxConcurrentOther:
        #    submitLastzJob(result_dir = result_dir, fa = oneQuery, qCount = j, rID = i, q = 'other', qsub = qsubMain)
        if getJobCount('all.q') < maxConcurrentBcl:
            submitLastzJob(result_dir = result_dir, fa = oneQuery, qCount = j, rID = i, q = 'all.q', qsub = qsubBcl)
        else:
            if debug:
                print("queuing")
            while True:
		if getJobCount('all.q') < maxConcurrentBcl:
                #if getJobCount('laird') < maxConcurrentLaird or getJobCount('other') < maxConcurrentOther:
                    break
#clean()
logging.info("All done")
