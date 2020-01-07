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
import re
import argparse
import cnv_similator
import time

def parse_region(input):
    m1 = re.search('(\d+):(\d+)-(\d+)',input)
    CHROM = int(m1.group(1))
    START = int(m1.group(2))
    END = int(m1.group(3))
    return CHROM, START, END


def cnvFastaSimulator(mututionInfo, refFasta, outdir):
    localtime = time.asctime( time.localtime(time.time()) )
    print('%s simluate mutation Fasta' % localtime)
    with open (mututionInfo) as f:
            for i in f:
                if i[0] == '#':
                    continue
                iterms = i.split()
                Type = iterms[0][0:3]
                print(iterms[3])
                chrom, start, end = parse_region(iterms[1])
                if iterms[2] != '0':
                    insert_chrom, insert_start, insert_end = parse_region(iterms[2])
                else:
                    insert_chrom, insert_start, insert_end = 0,0,0
                
                out_fq = '%s/%s' % (outdir, iterms[3])

                cnv_similator.cnvSimulator( refFasta, 
                out_fq, 
                chrom, 
                start, 
                end, 
                Type, 
                insert_start,
                insert_end,
                outdir)

    return 0


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-m",required=True, help='CNV information')
    parser.add_argument("-f",required=True, help='reference FASTA')
    parser.add_argument("-o",required=True, help='outdir')
    args = parser.parse_args()

    mututionInfo = args.m
    refFasta =  args.f
    outdir = args.o

    cnvFastaSimulator(mututionInfo, refFasta, outdir)

