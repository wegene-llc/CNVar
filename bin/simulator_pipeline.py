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
import CNV_Samples_similtor_kmer
import mix_alpha

def simulatorPipeline(RefSimList, outDir, ref, refDir, insertSize, region, kmer, MQ_cutoff ,threads):
    CNV_Samples_similtor_kmer.cnvSampleSimulatorKmer(RefSimList, outDir, refDir, kmer, insertSize, region, ref, MQ_cutoff ,threads)
    mix_reads_num = 0
    mix_alpha.mixAlpha(RefSimList, outDir, mix_reads_num, outDir, region, threads, ref)

    return 0

if __name__ == '__main__':

    if len(sys.argv) < 8:
        print('usage: python %s <RefSimList> <outdir> <ref.fa> <ref dir> <insert size> <region> <k-mers> <MQ_cutoff> <threads>' % sys.argv[0])
        exit(1)

    RefSimList = sys.argv[1]
    outDir = sys.argv[2]
    ref = sys.argv[3]
    refDir = sys.argv[4]
    insertSize = sys.argv[5]
    region = sys.argv[6]
    kmer = sys.argv[7]
    MQ_cutoff = sys.argv[8]
    threads = sys.argv[9]

    simulatorPipeline(RefSimList, outDir, ref, refDir, insertSize, region, kmer, MQ_cutoff ,threads)