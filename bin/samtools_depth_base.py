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
import anno_base_in_depth
import subprocess
import time
from config import PYTHONPATH, BWAPATH, SAMTOOLSPATH

def samtoolsDepthBase(bam, region, refGenome, outfilename):
    localtime = time.asctime(time.localtime(time.time()))
    depthfile = '%s.depth.txt' % outfilename
    cmd =  "%s/samtools depth -a %s -r %s > %s" % (SAMTOOLSPATH, bam, region, depthfile)
    print('%s CMD: %s' % (localtime, cmd) )
    output = subprocess.check_output(cmd, shell=True).decode()
    depthBaseFile = '%s.base.depth.txt' % outfilename
    anno_base_in_depth.annoBaseInDepth(depthfile, depthBaseFile, refGenome)

    return 0


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('usage: python %s <bam> <region> <refGenome> <outfilename>' % sys.argv[0])
        exit(1)
    
    bam = sys.argv[1]
    region = sys.argv[2]
    refGenome = sys.argv[3]
    outfilename = sys.argv[4]

    samtoolsDepthBase(bam, region, refGenome, outfilename)