#!/usr/bin/env python

from distutils.spawn import find_executable
import warnings
import os
import subprocess
import re
version = {'samtools':'0.1.19-44428cd', 'blastn':'2.2.31+', 'perl':'v5.16.0', 'bedtools':'v2.25.0'}

def getVersion(name):
    ret = ''
    if name == 'samtools':
	o = subprocess.check_output('samtools 2>&1; exit 0',shell=True)
	#Version: 0.1.19-44428cd
	ret = re.search('Version: (.*?)\n',o).groups()[0]
    elif name == 'blastn':
	o = subprocess.check_output('blastn -version',shell=True)
	#Version: 0.1.19-44428cd
	ret = re.search('blastn: (.*?)\n',o).groups()[0]
    elif name == 'perl':
	o = subprocess.check_output('perl -v', shell = True)
	ret = re.search('This is perl 5, version \d+, subversion \d+ \((.*?)\)',o).groups()[0]
    elif name == 'bedtools':
	o = subprocess.check_output('bedtools --version', shell = True)
	ret = re.search('bedtools (.*?)\n', o).groups()[0]
    else:
	raise Exception('NOT IMPLEMENTED for {0}'.format(name))
    return(ret)

for i in version:
    exe = find_executable(i)
    if not exe:
	raise Exception('{0} not found. Please install it first.'.format(i))
    v = getVersion(i)
    if v != version[i]:
	warnings.warn('''{0} version is {1} (not {2}).
	You can proceed, but it may not work as expected.'''.format(i,v, version[i]))
