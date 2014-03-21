#!/usr/bin/env python
'''
struct_chains2R.py <chain_prefix>
'''

import sys,os,re

prefix = sys.argv[1]

def find_files(prefix):
	'''
	will return a list of files in the current directory that contain the prefix
	'''
	return [f for f in os.listdir('.') if re.search(prefix,f)]

def parse_files(filelist,prefix):
	'''
	output a file prefix_sum.txt with the summary chains
	'''
	fo = open(prefix+'_sum.txt','w')
	for f in filelist:
		fi = open(f,'r')
		rep = f.split('_')[-1]
		tmp1 = [re.sub(':','',r.strip().lstrip()) for r in fi if re.search('[0-9]{3,}:',r)]
		tmp2 = [re.sub('[ ]{1,}',' ',r) for r in tmp1]
		tmp3 = [re.sub('--','- -',r) for r in tmp2]
		tmp4 = [fo.write(rep+' '+r+'\n') for r in tmp3]
		