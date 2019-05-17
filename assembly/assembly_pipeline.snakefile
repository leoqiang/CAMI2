#
# Assembly pipeline: metaspades + megahit + links + long reads correction
#
# Version 1.0
# - SPAdes 3.13.0
# - megahit 1.2.1 (beta)
# - LINKS 1.8.6
# - fmlrc 1.0.0
#
# ----------------------------------------
import os
import logging
import logging.config
import yaml
import pprint

# settings
DIR_WORK='asmb_pipeline'
THREADS=20
# input (config)
SAMPLE_NAME=config['sample']
DIR_ROOT=os.path.join(config['run_dir'], SAMPLE_NAME)
FASTQ_READSHORT=config['sr']
FASTQ_READSLONG=config['lr']



# ----------------------------
#  setup logging under SAMPLE_RESULTS_DIR
# ----------------------------
# using snakemake utility func 'srcdir' to get the directory of snakefile
logging_cfg_file = srcdir("logging.yml")
logging_cfg = yaml.load(open(logging_cfg_file, 'r'))
logging_cfg['handlers']['sample_debug'] = {
    'class': 'logging.FileHandler',
    'level': 'DEBUG',
    'formatter': 'basic',
    'filename': os.path.join(config['run_dir'], SAMPLE_NAME, SAMPLE_NAME+'_analysis.log'),
    'delay': 'True'
}
logging_cfg['root']['handlers'].append('sample_debug')
logging.config.dictConfig(logging_cfg)
logging.info(" Running metagenomics assembly pipeline sample {!r} ...".format(config['sample']))
logging.info(pprint.pformat(config))



# ----------------------------
# File & DIR definition
# ----------------------------
#
# S0: inputs handling
FASTQ_READSHORT_SL=os.path.join(DIR_ROOT, SAMPLE_NAME+"_sr_interleaved.fq.gz")
FASTA_READSLONG=os.path.join(DIR_ROOT, SAMPLE_NAME+"_lr.fa")
#
# S4: LINKS scaffolded assembly
FASTA_LINKS_MEGATHIT=os.path.join(DIR_ROOT, 'links', SAMPLE_NAME+'.scaffolds.fa')
DIR_TOOL_LINKS='/data/testdir/qiafu/tool_test/CAMI2-201903/tools/links_v1.8.6'
DIR_LINKS=os.path.join(DIR_ROOT, 'links')
#
# S3: FMLRC corrected long reads
FASTA_READSLONG_COR=os.path.join(DIR_ROOT, 'fmlrc', SAMPLE_NAME+'_fmlrc_lr_cor.fa')
DIR_TOOL_ROPEBWT='/data/testdir/qiafu/tool_test/CAMI2-201903/tools/fmlrc-longread.correction/ropebwt2'
DIR_TOOL_FMLRC='/data/testdir/qiafu/tool_test/CAMI2-201903/tools/fmlrc-longread.correction/fmlrc-1.0.0'
DIR_FMLRC=os.path.join(DIR_ROOT, 'fmlrc')
#
# S2: MEGAHIT sr.cor reads assembly
FASTA_MEGAHIT=os.path.join(DIR_ROOT, 'megahit', 'sense', 'megahit-sense.contigs.fa')
DIR_MEGAHIT=os.path.join(DIR_ROOT, 'megahit')
#
# S1: SPAdes corrected reads
FASTQ_READSHORT_base= os.path.splitext(os.path.basename(FASTQ_READSHORT))[0]
FASTQ_READSHORT_COR_P1=os.path.join(DIR_ROOT, 'spades_sr_cor', 'corrected', FASTQ_READSHORT_base+'_1.00.0_0.cor.fastq.gz')
FASTQ_READSHORT_COR_P2=os.path.join(DIR_ROOT, 'spades_sr_cor', 'corrected', FASTQ_READSHORT_base+'_2.00.0_0.cor.fastq.gz')
FASTQ_READSHORT_COR_SE=os.path.join(DIR_ROOT, 'spades_sr_cor', 'corrected', FASTQ_READSHORT_base+'__unpaired.00.0_0.cor.fastq.gz')
DIR_SPADES=os.path.join(DIR_ROOT, 'spades_sr_cor')


# ----------------------------
# Rules
# ----------------------------

rule all:
    input:
        FASTA_LINKS_MEGATHIT


rule inputs_handling:
    input:
        FASTQ_READSHORT,
        FASTQ_READSLONG
    output:
        FASTQ_READSHORT_SL,
        FASTA_READSLONG
    shell:
        """
        cd {DIR_ROOT}
        ln -s {input[0]} {output[0]}
        module load seqtk/1.2
        seqtk seq -A {input[1]} > {output[1]}
        module purge
        """


rule links_scarfolding:
    input:
        fasta_asmb=FASTA_MEGAHIT,
        fasta_lr=FASTA_READSLONG_COR
    output:
        FASTA_LINKS_MEGATHIT
    params:
        rundir=DIR_LINKS
    shell:
       """
       cd {params.rundir}
       echo {input.fasta_lr} > tmp_lr.fof
       {DIR_TOOL_LINKS}/LINKS -k 21 -t 20 -f {input.fasta_asmb} -s tmp_lr.fof -b {SAMPLE_NAME} -d 4000,5000,6000,8000,10000,15000,20000
       """


rule fmlrc_lr_cor:
    input:
        FASTQ_READSHORT,
        FASTA_READSLONG
    output:
        FASTA_READSLONG_COR
    params:
        rundir=DIR_FMLRC
    threads: config['threads']
    log:
        ropebwt=os.path.join(DIR_FMLRC, 'log/ropebwt.log'),
        fmlrc=os.path.join(DIR_FMLRC, 'log/fmlrc.log')
    shell:
        """
        cd {params.rundir}
        gunzip -c {input[0]} | awk 'NR % 4 == 2' | sort | tr NT TN | {DIR_TOOL_ROPEBWT}/ropebwt2 -LR | tr NT TN | {DIR_TOOL_FMLRC}/fmlrc-convert short_reads_msbwt.npy > {log.ropebwt} 2>&1 
        {DIR_TOOL_FMLRC}/fmlrc -p {threads} short_reads_msbwt.npy {input[1]} {output[0]} > {log.fmlrc} 2>&1 
        """


rule megahit:
    input:
        FASTQ_READSHORT_COR_P1,
        FASTQ_READSHORT_COR_P2,
        FASTQ_READSHORT_COR_SE
    output:
        FASTA_MEGAHIT
    params:
        rundir=DIR_MEGAHIT
    threads: config['threads']        
    shell:
        """
        cd {params.rundir}
        /data/testdir/qiafu/tool_test/CAMI2-201903/tools/MEGAHIT-1.2.1-beta-Linux-static/bin/megahit -t {threads} -1 {input[0]} -2 {input[1]} -r {input[2]} --presets meta-sensitive --out-prefix megahit-sense -o sense
        """


rule spades_sr_cor:
    input:
        FASTQ_READSHORT_SL
    output:
        FASTQ_READSHORT_COR_P1,
        FASTQ_READSHORT_COR_P2,
        FASTQ_READSHORT_COR_SE
    params:
        rundir=DIR_SPADES
    threads: config['threads']
    shell:
        """
        module load SPAdes/3.13.0
        cd {params.rundir}
        spades.py -t {threads} --only-error-correction --12 {input} -o .
        module purge
        """
