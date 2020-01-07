#! /usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################################################
#    CNVar: accurate knowledge-based copy number variation prediction using max likelihood method       #
#                                                                                                       # 
#    Copyright (C) 2019 Enming He(emhe@wegene.com)                                                      #
#                                                                                                       #
#    This program is free software: you can redistribute it and/or modify                               #
#    it under the terms of the GNU General Public License as published by                               #
#    the Free Software Foundation, either version 3 of the License, or                                  #
#    (at your option) any later version.                                                                #
#                                                                                                       #
#    This program is distributed in the hope that it will be useful,                                    #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                                     #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                      #
#    GNU General Public License for more details.                                                       #
#                                                                                                       #
#    You should have received a copy of the GNU General Public License                                  #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.                             #
#########################################################################################################

import os
import sys
import subprocess
from config import PYTHONPATH, BWAPATH, SAMTOOLSPATH
import time 


def run(cmd):
    print('CMD:%s' % cmd)
    signal = subprocess.check_call(cmd, shell=True).decode()
    return signal


def sed(infile):
	with open(infile + '.tmp','w') as w:
		with open(infile) as f:
			line = f.readline()
			line = line.replace(':', ' ')
			w.write(line)

			for line in f:
				w.write(line)
	
	os.rename(infile + '.tmp', infile)

	return 0


def samtoolsFaidx(refsimlist, dirname, region):
	localtime = time.asctime(time.localtime(time.time()))
	print('%s extract target region in Fasta' % localtime)
	with open(refsimlist) as f:
		for line in f:
			name = line.strip()
			print(name)
			# almost equal to "sed -i 's1/:/ /g' $dir/$i"
			sed('%s/%s' % (dirname, name))
			cmd = '%s/samtools faidx %s/%s %s > %s/tmp' % (SAMTOOLSPATH,dirname, name, region, dirname)
			subprocess.check_call(cmd, shell=True)
			os.rename('%s/tmp' % dirname,  '%s/%s' % (dirname, name))
	
	return 0


if __name__ == '__main__':
	if len(sys.argv):
		print("usage: python %s <refsimlist> <dir> <region>" % sys.argv[0])
		exit(1)

	refsimlist = sys.argv[1]
	dirname = sys.argv[2]
	region = sys.argv[3]

	samtoolsFaidx(refsimlist, dirname, region)
