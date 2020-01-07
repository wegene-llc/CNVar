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
from itertools import combinations
import subprocess
from multiprocessing import Pool
from config import PYTHONPATH, BWAPATH, SAMTOOLSPATH
import shutil
import time
import anno_base_in_depth

def run(cmd):
    print('CMD:%s' % cmd)
    signal = subprocess.check_call(cmd, shell=True)
    return signal


def samtoolsMergeIndexDepthAnno( outdir, pairs, source_dir, region, refGenome ):
    localtime = time.asctime(time.localtime(time.time()))
    print('%s mix %s with %s' % (localtime, pairs[0], pairs[1]))
    dirname = outdir + "/" + pairs[0] + "_" + pairs[1]
    if os.path.exists(dirname):
        shutil.rmtree(dirname, True)
        os.mkdir(dirname)
    else:
        os.mkdir(dirname)

    cmd = "%s/samtools merge -@ 1 -f %s/%s_%s.bam %s/%s_homo/%s_homo.bam %s/%s_homo/%s_homo.bam" % (SAMTOOLSPATH, 
    dirname,pairs[0],pairs[1],source_dir,pairs[0],pairs[0],source_dir,pairs[1],pairs[1])
    run(cmd)

    cmd = "%s/samtools index %s/%s_%s.bam" % (SAMTOOLSPATH, dirname, pairs[0],pairs[1]) 
    run(cmd)

    cmd = "%s/samtools depth -a -r %s %s/%s_%s.bam > %s/%s_%s.bam.depth.txt" % (SAMTOOLSPATH, region, dirname,
    pairs[0],pairs[1],dirname, pairs[0],pairs[1]) 
    run(cmd)

    anno_base_in_depth.annoBaseInDepth('%s/%s_%s.bam.depth.txt' % (dirname, pairs[0],pairs[1]), 
    '%s/%s_%s.bam.base.depth.txt' % (dirname, pairs[0],pairs[1]), refGenome)
    
    return 0


def mixAlpha(refsimlist, outdir, mix_reads_num, source_dir, region, threads, refGenome):
    threads = int(threads)
    List = []
    with open(refsimlist) as f:
        for i in f:
            List.append(i.strip())

    mix = list(combinations(List, 2)) 

    #run Multiprocesing
    p = Pool(int(threads))
    for pairs in mix:
        result = p.apply_async(samtoolsMergeIndexDepthAnno, args=(outdir, pairs, source_dir, region, refGenome))
    result.get()  
    p.close()
    p.join()


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("usage:python mix_alpha.py <RefSimList.txt> <outdir> <mix reads number> <source dir>  <region> <threads>, <refGenome>")
        exit(1)

    refsimlist = sys.argv[1]
    outdir = sys.argv[2]
    mix_reads_num = sys.argv[3]
    source_dir = sys.argv[4]
    region = sys.argv[5]
    threads = sys.argv[6]
    refGenome = sys.argv[7]

    mixAlpha(refsimlist, outdir, mix_reads_num, source_dir, region, threads, refGenome)


