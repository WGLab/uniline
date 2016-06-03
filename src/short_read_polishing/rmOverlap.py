#!/usr/bin/env python
import sys
sys.path.append('/home/yunfeiguo/projects/PacBio_reference_genome/short_read_polishing')
import utils
if len(sys.argv) >= 3:
    utils.generateNonOverlapVCF(sys.argv[1], sys.argv[2])
else:    
    utils.generateNonOverlapVCF(sys.argv[1])
