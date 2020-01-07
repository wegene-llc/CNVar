#### 0. Quick test
```
sh example/example.sh
```


#### 1. Prerequisites:

 python==3.6.5

 bwa==0.7.10-r789

 samtools==1.9

Python packages are mentioned in requirements.txt 

Please specify paths in bin/config.py

#### 2. Prepare your CNV information.

Present your known CNV information in a 'mutation info file', use the format below:

```
#MutationType Position InsertRegion Name
deletion 16:215400-234700 0 SEA
duplication 16:203016-213826 16:300000-300000 anti11
```

The first column is the mutation type, which is deletion or duplication.

The second column is the mutation region.

The third column is the start position where the duplication sequence inserts. if a region is inputted, the start position will be randomly chosen in the input region. 0 for deletions. Generally, insert position close to mutation region is recommended (for example: within 10kbp), makes the simulate region not being too large. 

The fourth column is the name of the mutation.


For now, we only analyze CNVs in one chromosome at once.
And please notice as the number of CNVs in the mutation format file increases, the runtime of each sample will raise exponentially. To make the computing effective, please try to put unrelated CNVs into different mutation info files (In the next version of CNVar, we will try to make this automatic).


#### 3. Simulate references

For each mutation info file, generate a set of references.

Example:
```
python CNVar.py simulate \
--mutation_info example/example.cnv.info.txt \
--reference_genome < reference genome: hg19 Fasta > \
--read_length 150 \
--insert_size 300 \
--outdir example/example_simulate
--MQ_cutoff <Mapping Quality cutoff> # optional
```

Paths of references are recorded in <outdir>/references.list

Notice:

(1) Only hg19 reference genome has been tested. The same below.

(2) Mapping Quality cutoff with default value 10, to filter out reads with Mapping Quality lower than 10, such as the multiple mapping reads. The same below.


#### 4. Genotype Calling

If you not sure the average depth and the GC bias of the bam. CNVar will calculate these by sampling the chromosome in the CNV information file. ( please don't test this with the example/example.bam, for it does not contain a full chromosome )

Example:
```
python CNVar.py call \
--mutation_info example/example.cnv.info.txt \
--reference_list example/example_simulate/references.list \
--reference_genome < reference genome: hg19 Fasta > \
--outdir <ourdir> \
--bam <bam> \ 
--gc_config <GC config file> \ # optional
--average_depth <average_depth> # optional
--MQ_cutoff <Mapping Quality cutoff> # optional
```

Notice:

(1) This command will generate a GC config file in JSON format, represents the GC bias of the chromosome in the mutation info file. You can input your GC bias info by using --gc_config <GC config file> to skip the GC bias sampling step in CNVar.

(2) If you already know the average depth of the input bam, input it with --average_depth <average_depth> or CNVar will calculate it by sampling the whole chromosome.


#### 5. Result

The result is presented in VCF4.2 format.

Example:
```
##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=LEN,Number=1,Type=String,Description="Length of SV">
##INFO=<ID=VARIANT>,Number=1,Type=String,Description="Name of Variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample
16 215400 1 - <DEL> . PASS SVTYPE=DEL;SVMETHOD=CNVar;END=234700;LEN=19300;VARIANT=SEA GT 1/1 example.bam
```
