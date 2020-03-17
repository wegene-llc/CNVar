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
import subprocess
from pathlib import Path
import samtools_faidx
import simulator_pipeline
import multiprocessing
import calculate_samtools_depth
import samtools_depth_base
import read_depth_method_matrix_mem_optimize_GCmodify_GCprepare
import possion_FixPossionModel_siteBYsite_optimal_GCmodify
import format_vcf
import cnv_fasta_simulator
from config import PYTHONPATH, BWAPATH, SAMTOOLSPATH
import time

def get_region(infile):

    def parse_region(string):
        m = re.search('(\w+):(\d+)-(\d+)', string)
        chrom = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        length = end - start

        return (chrom, start, end, length)

    # collect region info in Mutation info
    Dict = {}
    Dict_detect = {}
    with open(infile) as f:
        for line in f:
            if line[0] == '#':
                continue

            cnv_type, mutation_region, background_region, cnv_name = line.strip().split()

            if re.search('_', cnv_name):
                print('Error: symbol \'_\' should not be used in variant name %s' % cnv_name)
                exit(1)
            if re.search('homo', cnv_name):
                print('Error: keyword \'homo\' should not be used in variant name %s' % cnv_name)
                exit(1)

            parsed_mutation_region = parse_region(mutation_region)

            if background_region != '0':
                parsed_background_region = parse_region(background_region)
                if parsed_mutation_region[0] != parsed_background_region[0]:
                    print('Erorr: In this version, we only accept duplication insert in the same chromesome.')
                    exit(1)
                else:
                    if parsed_mutation_region[0] not in Dict:
                        Dict[parsed_mutation_region[0]] = []
                    Dict[parsed_mutation_region[0]].append(parsed_mutation_region[1])
                    Dict[parsed_mutation_region[0]].append(parsed_mutation_region[2])
                    Dict[parsed_mutation_region[0]].append(parsed_background_region[1])
                    Dict[parsed_mutation_region[0]].append(parsed_background_region[2]+parsed_mutation_region[3])
            else:
                if parsed_mutation_region[0] not in Dict:
                    Dict[parsed_mutation_region[0]] = []
                Dict[parsed_mutation_region[0]].append(parsed_mutation_region[1])
                Dict[parsed_mutation_region[0]].append(parsed_mutation_region[2])

            if parsed_mutation_region[0] not in Dict_detect:
                Dict_detect[parsed_mutation_region[0]] = []
            Dict_detect[parsed_mutation_region[0]].append(parsed_mutation_region[1])
            Dict_detect[parsed_mutation_region[0]].append(parsed_mutation_region[2])

    # summary the regions may be used
    def extract_region(Dict):
        Dict_region = {}
        for i in Dict:
            Dict_region[i] = (min(Dict[i]) , max(Dict[i]))

        if len(Dict_region) > 1:
            print('Error: In this version, we only comparing CNVs in the one chromesome at once.')
            exit(1)

        for i in Dict_region:
            chrom = i
            start = Dict_region[i][0] - 1000 # extend 1kbp backward
            if start < 0:
                start = 0
            end = Dict_region[i][1] + 1000 # extend 1kbp forward

        region = (chrom, start, end)
        return region

    total_region = extract_region(Dict)
    detect_region = extract_region(Dict_detect)

    return detect_region, total_region


def run(cmd):
    localtime = time.asctime( time.localtime(time.time()) )
    print('%s CMD:%s' % (localtime, cmd))
    signal = subprocess.check_call(cmd, shell=True)
    return signal


def simulate(args):
    MQ_cutoff = int(args.MQ_cutoff)
    if Path(args.outdir).is_dir():
        myoutdir = os.path.abspath(args.outdir)
    else:
        os.mkdir(args.outdir)
        myoutdir = os.path.abspath(args.outdir)
    dirname, filename = os.path.split(os.path.abspath(sys.argv[0]))

    if re.search('homo',myoutdir):
        print('Error: keyword \'homo\' should not use in dirname %s ' % myoutdir)
        exit(1)

    #get template
    detect_region,region = get_region(args.mutation_info)
    chrom = region[0]

    cmd = '%s/samtools faidx %s %s > %s/wildtype' % (SAMTOOLSPATH, args.reference_genome, chrom, myoutdir)  
    run(cmd)

    # simulate fasta
    cnv_fasta_simulator.cnvFastaSimulator(args.mutation_info, '%s/wildtype' % myoutdir, myoutdir)

    with open('%s/refsimlist' % myoutdir, 'w') as w:
        with open(args.mutation_info) as f:
            for line in f:
                if line[0] == '#':
                    continue
                items = line.strip().split()
                w.write(items[3] + '\n')
        w.write('wildtype\n')

    # extend simulation region 1Kbp forward and backward
    if int(region[1]) - 1000 < 0:
        extend_start = 0
    else:
        extend_start =  int(region[1]) - 1000
    
    extend_end = int(region[2]) + 1000

    samtools_faidx.samtoolsFaidx('%s/refsimlist' % myoutdir, myoutdir, '%s:%s-%s' % (region[0],extend_start,extend_end) )

    #simulate 'k-mer' samples containing mutations
    simulator_pipeline.simulatorPipeline('%s/refsimlist' % myoutdir, myoutdir, args.reference_genome, myoutdir, args.insert_size, 
         '%s:%s-%s' % (detect_region[0],detect_region[1],detect_region[2]), args.read_length, MQ_cutoff , args.threads)

    # equal to "find %s -name '*gz' > %s/references.list" % (myoutdir, myoutdir)
    with open('%s/references.list' % myoutdir, 'w') as w:
        for root,dirs,files in os.walk(myoutdir):
            for file in files: 
                filename = os.path.join(root,file)
                if filename.endswith('.gz'):
                    w.write(filename + '\n')

    w = open('%s/tmp' % myoutdir,'w')
    with open('%s/references.list' % myoutdir) as f:
        for line in f:
            if re.search('homo',line):
                d = int(args.read_length) * 2
            else:
                d = int(args.read_length) * 4

            w.write('%s\t%s\n' % (line.strip(), d))
    w.close()

    os.rename('%s/tmp' % myoutdir, '%s/references.list' % myoutdir)

    return 0


def call(args):
    MQ_cutoff = int(args.MQ_cutoff)

    if Path(args.outdir).is_dir():
        myoutdir = os.path.abspath(args.outdir)
    else:
        os.mkdir(args.outdir)
        myoutdir = os.path.abspath(args.outdir)

    dirname, filename = os.path.split(os.path.abspath(sys.argv[0]))

    #if bai exists:
    baiName1 = os.path.splitext(args.bam)[0] + '.bai'
    baiName2 = args.bam + '.bai'

    if os.path.exists(baiName1) or os.path.exists(baiName2):
        pass
    else:
        cmd = '%s/samtools index %s' % (SAMTOOLSPATH, args.bam)
        run(cmd)

    #get region
    detect_region,region = get_region(args.mutation_info)
    chrom = region[0]

    outfilename = '%s/%s' % (myoutdir,os.path.splitext(args.bam)[0].split('/')[-1])


    if args.gc_config and args.average_depth:
        start_pos = detect_region[1] - 1000
        end_pos = detect_region[2] + 1000
        # if MQ_cutoff != 60, filt out reads MQ lower than cutoff
        if MQ_cutoff ==  60: 
            samtools_depth_base.samtoolsDepthBase(args.bam, '%s:%s-%s' % (detect_region[0],detect_region[1],detect_region[2]), args.reference_genome, outfilename)
            with open('%s.sample.list' % outfilename,'w') as w:
                w.write('%s.base.depth.txt.gz\n' % outfilename) 
              
            samplelist = '%s.sample.list' % outfilename       
        else:
            chrBam = '%s.%s.bam' % (outfilename, detect_region[0])
            cmd = '%s/samtools view -bh -q %s %s %s:%s-%s > %s' % (SAMTOOLSPATH, MQ_cutoff, args.bam, detect_region[0],start_pos, end_pos, chrBam)
            run(cmd)
            cmd = '%s/samtools index %s' % (SAMTOOLSPATH, chrBam)
            run(cmd)
            samtools_depth_base.samtoolsDepthBase(chrBam, '%s:%s-%s' % (detect_region[0],detect_region[1],detect_region[2]), args.reference_genome, chrBam)
            with open('%s.sample.list' % chrBam,'w') as w:
                w.write('%s.base.depth.txt.gz\n' % chrBam) 

            samplelist = '%s.sample.list' % chrBam
    else:
        # chrLength
        bamHeader = subprocess.check_output('%s/samtools view -H %s|grep SN:%s' % (SAMTOOLSPATH, args.bam, chrom), shell=True).decode()
        m = re.search('LN:(\d+)', bamHeader)
        chrLength = m.group(1)

        chrBam = '%s.%s.bam' % (outfilename, detect_region[0])

        if MQ_cutoff == 60:
            cmd = '%s/samtools view -bh %s %s > %s' % (SAMTOOLSPATH, args.bam, detect_region[0], chrBam)
        else:
            cmd = '%s/samtools view -bh -q %s %s %s > %s' % (SAMTOOLSPATH, MQ_cutoff, args.bam, detect_region[0], chrBam)
        run(cmd)

        cmd = '%s/samtools index %s' % (SAMTOOLSPATH, chrBam)
        run(cmd)

        samtools_depth_base.samtoolsDepthBase(chrBam, detect_region[0], args.reference_genome, chrBam)

        with open('%s.sample.list' % chrBam,'w') as w:
            w.write('%s.base.depth.txt.gz\n' % chrBam) 

        samplelist = '%s.sample.list' % chrBam

    if args.gc_config:
        gc_config = args.gc_config
    else:
        read_depth_method_matrix_mem_optimize_GCmodify_GCprepare.readDepthMethodMatrixMemOptimizeGCmodifyGCprepare('%s.base.depth.txt.gz' % chrBam, 
        args.binsize, '%s.GCconfig.json' % chrBam, args.slide_window, '%s:%s-%s' % (chrom, 1000, int(chrLength)-1000))

        gc_config = '%s.GCconfig.json' % (chrBam)

    if args.average_depth:
        average_depth = args.average_depth
    else:
        average_depth = calculate_samtools_depth.calculateSamtoolsDepth('%s.base.depth.txt.gz' % chrBam)
        print('average depth: %s ' % average_depth)

    start_pos = detect_region[1] + 1000
    end_pos = detect_region[2] - 1000

    start_pos2 = start_pos - start_pos % 50
    end_pos2 = end_pos - end_pos % 50

    possion_threads = 1
    possion_FixPossionModel_siteBYsite_optimal_GCmodify.possionFixPossionModelSiteBYSiteOptionGCmodify(args.reference_list, samplelist, '%s.result.txt' % outfilename, 
    possion_threads, '%s:%s-%s' % (detect_region[0], start_pos2, end_pos2), average_depth, gc_config, args.binsize, args.slide_window, 
    args.logmin, args.mismatch )

    format_vcf.cnvar2vcf( '%s.result.txt' % outfilename, '%s.vcf' % outfilename, args.mutation_info )

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='\
    CNVar: accurate knowledge-based copy number variation prediction using max likelihood method')

    subparsers = parser.add_subparsers()

    #simulate references
    parser_list = subparsers.add_parser('simulate', help='simulate references')
    parser_list.add_argument('-mt','--mutation_info',help="Mutation info <FILE>",required=True)
    parser_list.add_argument('-rg','--reference_genome', help='Reference genome (Fasta) <FILE>',required=True)
    parser_list.add_argument('-rl','--read_length', help="read length (bp) <INT>",required=True)
    parser_list.add_argument('-is','--insert_size', help="insert size (bp) <INT>",required=True)
    parser_list.add_argument('-o', '--outdir', help="output directory ",required=True)
    parser_list.add_argument('-t', '--threads', help="threads <INT>, default use all CPUs ", default=multiprocessing.cpu_count())
    parser_list.add_argument('-mq', '--MQ_cutoff', help="Mapping Quality cutoff <INT>, default: 10", default=10)
    parser_list.set_defaults(func=simulate)

    #calling
    parser_list = subparsers.add_parser('call', help='genotype calling')
    parser_list.add_argument('-mt','--mutation_info',help="Mutation info <FILE>",required=True)
    parser_list.add_argument('-r','--reference_list', help='Reference list <FILE>',required=True)
    parser_list.add_argument('-b','--bam', help="bam <FILE>",required=True)
    parser_list.add_argument('-c', '--gc_config', help="GC config (json). If absent, GC config will be calculate with the bam file. <FILE>")
    parser_list.add_argument('-d', '--average_depth', help="Average Depth. If absent. Average Depth will be calculate with the bam file. <FLOAT>")
    parser_list.add_argument('-l', '--logmin', help="minimal logP, default:1e-7 <FLOAT>",default=1e-7)
    parser_list.add_argument('-ms', '--mismatch', help='assume mismatch rate, default:1e-4 <FLOAT>', default=1e-4)
    parser_list.add_argument('-rg','--reference_genome', help='Reference genome (Fasta) <FILE>',required=True)
    parser_list.add_argument('-bz', '--binsize', help="binsize（bp). Develop only. default:150 <INT>",default=150)
    parser_list.add_argument('-sw', '--slide_window', help="slide window (bp). Develop only. default:50 <INT>",default=50)
    parser_list.add_argument('-mq', '--MQ_cutoff', help="Mapping Quality cutoff <INT>, default: 10", default=10)
    parser_list.add_argument('-o', '--outdir', help="output directory",required=True)

    parser_list.set_defaults(func=call)

    arglist = sys.argv[1:]
    if len(arglist) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args(arglist)
        args.func(args)
