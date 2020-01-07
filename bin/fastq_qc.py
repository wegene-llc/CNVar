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
from pathlib import Path
from config import PYTHONPATH, BWAPATH, SAMTOOLSPATH
import subprocess

def fastqBwaSort(fq1, fq2, name, refgenome, MQ_cutoff, threads):
	MQ_cutoff = int(MQ_cutoff)
	threads = int(threads)

	if Path(name).is_dir():
		pass
	else:
		os.mkdir(name)
	os.chdir(name)

	if MQ_cutoff == 60:
		cmd = '%s/bwa mem -M -t %s -R \"@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA\"  %s %s %s |%s/samtools sort -@ %s > %s.bam' % (BWAPATH, 
		threads, name, name, name, refgenome, fq1, fq2, SAMTOOLSPATH, threads, name)
	else:
		cmd = '%s/bwa mem -M -t %s -R \"@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA\"  %s %s %s |%s/samtools view  -bh -q %s |%s/samtools sort -@ %s > %s.bam' \
			% (BWAPATH, threads, name, name, name, refgenome, fq1, fq2, SAMTOOLSPATH ,MQ_cutoff ,SAMTOOLSPATH, threads, name)
	
	print(cmd)
	subprocess.check_call(cmd, shell=True)

	subprocess.check_call('%s/samtools index %s.bam' % (SAMTOOLSPATH, name), shell=True)

	return 0


if __name__ == '__main__':
	if len(sys.argv) < 5:
		print('usage: python %s <abs_F1> <abs_F2> <NAME> <refgenome> <MQ cutoff> <threads>' % sys.argv[0])
		exit(1)
	
	fq1 = sys.argv[1]
	fq2 = sys.argv[2]
	name = sys.argv[3]
	refgenome = sys.argv[4]
	MQ_cutoff = sys.argv[5]
	threads = sys.argv[6]

	fastqBwaSort(fq1, fq2, name, refgenome, MQ_cutoff, threads)
    
