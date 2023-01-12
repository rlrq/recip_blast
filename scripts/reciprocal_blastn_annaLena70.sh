#!/bin/bash

## functions for running blastn to Anna-Lena's RenSeq sequences, and then reciprocal blast back to Col-0 NLR

source /mnt/chaelab/rachelle/src/run_blast.sh

## get sequences and cat to single file
## arg1: directory
## arg2: file name prefix
## arg3: genes
## arg4: acc
## arg5: feature
## arg6: output file name (optional)
get_seq () {
    local dir=${1}
    local prefix=${2}
    local fa=${dir}/${prefix}.fa
    local genes=${3}
    local acc=${4}
    local feature=${5}
    local complete=${6}
    local adj_dir=${7}
    local fasta=${8:-${fa}sta}
    local start_t=$(date '+%Y-%m-%d %T')
    if [ "${complete}" == 'True' ] && [ "${adj_dir}" == 'True' ]; then
        /mnt/chaelab/rachelle/scripts/get_seqs/get_seqs.sh -g ${genes} -a ${acc} -d ${dir} -f ${feature} --complete --adjust-dir --out ${fasta}
    elif [ "${complete}" == 'True' ]; then
        /mnt/chaelab/rachelle/scripts/get_seqs/get_seqs.sh -g ${genes} -a ${acc} -d ${dir} -f ${feature} --complete --out ${fasta}
    elif [ "${adj_dir}" == 'True' ]; then
        /mnt/chaelab/rachelle/scripts/get_seqs/get_seqs.sh -g ${genes} -a ${acc} -d ${dir} -f ${feature} --adjust-dir --out ${fasta}
    else
        /mnt/chaelab/rachelle/scripts/get_seqs/get_seqs.sh -g ${genes} -a ${acc} -d ${dir} -f ${feature} --out ${fasta}
    fi
}

## get sequences of contigs from output of blastn_to_al70 for contigs that have alignments with min len and id
## arg1: output directory
## arg2: .tsv file produced by blastn_to_al70 (blast --> asn --> tsv)
## arg3: output prefix
## arg4: minimum alignment length (bases)
## arg5: minimum alignment identity (%)
## arg6: output file name (optional)
filter_contig() {
    local dir=${1}
    local prefix=${2}
    local ftsv=${3}
    local minlen=${4}
    local minid=${5}
    local all_contigs="/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta"
    local out_root="${prefix}.minid${minid}_minlen${minlen}.contigs"
    local contig_names="${dir}/${out_root}.names.txt"
    awk -v id=${minid} -v len=${minlen} '{ if( $3 >= id && $4 >= len ) print }' $ftsv | cut -f2 | sort | uniq > $contig_names
    ## get contig sequences
    local contig_hits=${6:-${out_root}.fasta}
    python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/TE_SV/src'); import get_contigs; get_contigs.get_contigs('${all_contigs}','${contig_names}','${contig_hits}')"
}

###############
##  EXAMPLE  ##
###############

# ## getting sequences
# testdir=/mnt/chaelab/rachelle/hap_tags/test_src
# testprefix=testdm2
# get_seq $testdir $testprefix AT3G44400,AT3G44480,AT3G44630,AT3G44670 ref gene

# ## first blast
# fasta=/mnt/chaelab/rachelle/hap_tags/test_src/testdm2.fasta
# testdir=${testdir}/results
# testprefix=${testprefix}_al70
# blastn_to_al70 $testdir $testprefix $fasta

# ## filter + reciprocal blast
# minlen=1000
# minid=90
# blast1_tsv=${testdir}/${testprefix}.blastn_summary.tsv ## output of first blast
# filter_contig ${testdir} ${testprefix}.blastn ${blast1_tsv} ${minlen} ${minid}
# contigs_fasta=${testdir}/${testprefix}.blastn.minid${minid}_minlen${minlen}.contigs.fasta ## output of filter_cotig
# blastn_to_col0nlr ${testdir} ${testprefix} ${contigs_fasta}
