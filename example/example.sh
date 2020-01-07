#! /bin/bash
# author: Enming He
# date: 2019.11.30
set -e
dir=`pwd`

#simulate
python bin/CNVar.py simulate \
-mt example/example.cnv.info.txt \
-rg $dir/example/example.fasta \
-rl 150 -is 300 \
-o example/example_simulate

#call
python bin/CNVar.py call \
-mt example/example.cnv.info.txt \
-rg $dir/example/example.fasta \
-r example/example_simulate/references.list \
-b example/example.bam \
-c example/WGS_workflows_testData.GCprofile.json \
-d 39 \
-o example/example_call


# result

echo "

DONE
##########################################################
result is show in:  $dir/example/example_call/example.vcf
##########################################################
"