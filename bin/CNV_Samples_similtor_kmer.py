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
import k_mers
from multiprocessing import Pool
import fastq_qc
import samtools_depth_base
import shutil
import time

def bwaSamtoolsSortDepth(fq1, fq2, name, refgenome, MQ_cutoff, threads,bam, region, outfilename):
    fastq_qc.fastqBwaSort(fq1, fq2, name, refgenome, MQ_cutoff, threads)
    samtools_depth_base.samtoolsDepthBase(bam, region, refgenome, outfilename)

    return 0

def cnvSampleSimulatorKmer(List, outdir, refdir, k, insert_size, region, refgenome, MQ_cutoff ,threads):
    localtime = time.asctime(time.localtime(time.time()))
    print('%s simulate kmer reads' % localtime)
    # kmer reads
    p = Pool(int(threads))
    with open(List) as f:
        for line in f:
            name = line.strip()
            print(name)
            fa = '%s/%s' % (refdir, name)
            #step1
            genotype = 'homo'
            if os.path.exists('%s/%s_%s' % (outdir, name, genotype)):
                shutil.rmtree('%s/%s_%s' % (outdir, name, genotype), True)
                os.makedirs( '%s/%s_%s' % (outdir, name, genotype) )
            else:
                os.makedirs( '%s/%s_%s' % (outdir, name, genotype) )
            outfile1 = '%s/%s_%s/%s_%s_1.fq' % (outdir, name, genotype, name, genotype)
            outfile2 = '%s/%s_%s/%s_%s_2.fq' % (outdir, name, genotype, name, genotype)
            result = p.apply_async(k_mers.kMers, args=(fa, outfile1, outfile2, k, insert_size))
        result.get()

    p.close()
    p.join()

    # bwa AND samtools sort AND samtools depth base
    p2 = Pool(int(threads))
    with open(List) as f:
        for line in f:
            name = line.strip()
            genotype = 'homo'
            fq1 = '%s/%s_%s/%s_%s_1.fq' % (outdir, name, genotype, name, genotype)
            fq2 = '%s/%s_%s/%s_%s_2.fq' % (outdir, name, genotype, name, genotype)
            bam = '%s/%s_%s/%s_%s.bam' % (outdir, name, genotype, name, genotype)
            outfilename =  '%s' % bam
            result = p2.apply_async(bwaSamtoolsSortDepth, args=(fq1, fq2, '%s/%s_%s/%s_%s' % (outdir, name, genotype, name, genotype), 
            refgenome, MQ_cutoff, 1 ,bam, region, outfilename)) # threads 1
        result.get()
    p2.close()
    p2.join()


if __name__ == '__main__':
    if len(sys.argv) < 8:
        print('usage: python %s <List> <outdir> <ref dir> <k-value> <insert size> <region> <threads>' % sys.argv[0])
        exit(1)

    List = sys.argv[1]
    outdir = sys.argv[2]
    refdir = sys.argv[3]
    k = sys.argv[4]
    insert_size = sys.argv[5]
    region = sys.argv[6]
    refgenome = sys.argv[7]
    threads = sys.argv[8]

    cnvSampleSimulatorKmer(List, outdir, refdir, k, insert_size, region, refgenome, MQ_cutoff ,threads)