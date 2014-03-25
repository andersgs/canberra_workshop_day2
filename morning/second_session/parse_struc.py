#!/usr/bin/env python
'''
parse_struc.py <chain_prefix>

e.g.,

./parse_struc.py chain_K2
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
	file_count = 0
	out = []
	for f in filelist:
		fi = open(f,'r')
		rep = f.split('_')[-1]
		tmp1 = []
		count = 0
		for line in fi:
			tmp = line.strip().lstrip()
			if re.search('[0-9]{3,}:',tmp):
				tmp1.append(re.sub(':','',tmp))
			else:
				if re.search('BURNIN completed',tmp):
					count = 1
				else:
					if count==1 and re.search('Rep#:',tmp) and file_count==0:
						header=tmp
						header = re.sub('[ ]{2,}',';',header)
						header = re.sub('[ ,]{1}','_',header)
						header = re.sub(';',' ',header)
						header = re.sub('Rep#:','Step',header)
						header = 'Rep '+ header
						count +=1
						file_count+=1
					else:
						continue
		#tmp1 = [re.sub(':','',r.strip().lstrip()) for r in fi if re.search('[0-9]{3,}:',r)]
		tmp2 = [re.sub('[ ]{1,}',' ',r) for r in tmp1]
		tmp3 = [re.sub('--','- -',r) for r in tmp2]
		tmp4 = [rep+' '+r+'\n' for r in tmp3]
		out.extend(tmp4)
	
	fo.write(header+'\n')
	for r in out:
		fo.write(r)
	fo.close()
	return out
	
if __name__=="__main__":
	flist = find_files(prefix)
	parse_files(flist,prefix)
