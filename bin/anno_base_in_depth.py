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
import gzip
from config import  BWAPATH, SAMTOOLSPATH

def getRegionFromDepthfile(depthfile):
    with open(depthfile) as f:
        firstline = f.readline()
        chrom, start = firstline.strip().split()[0:2]
        lastline = firstline
        for line in f:
            lastline = line
        end = lastline.strip().split()[1]

    return chrom, int(start), int(end)


def annoBaseInDepth(depthfile, outfile, refGenome):

    chrom, start, end = getRegionFromDepthfile(depthfile)
    fasta = subprocess.check_output('%s/samtools faidx %s %s:%s-%s' % (SAMTOOLSPATH, refGenome, chrom, start, end), shell=True).decode()
    string = ''
    for line in fasta.split()[1:]:
        string += line

    i = 0
    with gzip.open(outfile + '.gz','wt') as o:
        with open(depthfile) as f:
            for line in f:
                o.write(line.strip() + '\t' + string[i] + '\n')
                i += 1

    return 0


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('usage: python %s <depth.txt> <outfile> <refGenome>' % ( sys.argv[0]))
        exit(1)
    
    depthfile = sys.argv[1]
    outfile = sys.argv[2]
    refGenome = sys.argv[3]
    
    annoBaseInDepth(depthfile, outfile, refGenome)
        