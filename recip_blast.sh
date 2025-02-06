#!/bin/bash

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

source /mnt/chaelab/rachelle/src/run_blast6.sh
source "${SCRIPT_DIR}/scripts/reciprocal_blastn_annaLena70.sh"

## if no arguments, show manual ## TODO
if [[ $# -eq 0 ]]; then
    man -l ${SCRIPT_DIR}/MANUAL_recip_blast.1
    exit 1
fi

DB_ANNALENA='/mnt/chaelab/shared/blastdb/anna_lena70.contigs/anna_lena70.contigs'
DB_COL0NLR='/mnt/chaelab/shared/blastdb/nlr164/all_NLR.164.db'
DB1_DEFAULT="${DB_ANNALENA}"
DB2_DEFAULT="${DB_COL0NLR}"
TASK_DEFAULT='default'
MINLEN_DEFAULT=0
MINID_DEFAULT=0
FEATURE_DEFAULT='gene'

params=${@}

while (( "$#" )); do
    case "$1" in
        -g|--gene) GENE="${2}";;
        -f|--fasta) FASTA="${2}";;
        -p|--prefix) PREFIX="${2}";;
        -d|--dir) DIR="${2}";;
        --task1) TASK1="${2}";; ## valid tasks:'blastn', 'default' (megablast), 'dc' (discontiguous megablast)
        --task2) TASK2="${2}";;
        --db1) DB1="${2}";;
        --db2) DB2="${2}";;
        --minlen) MINLEN="${2}";;
        --minid) MINID="${2}";;
        --complete) COMPLETE='True';;
        --feature) FEATURE="${2}";;
        --adjust-dir) ADJ_DIR='True';;
        -h|--help) man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1; exit 0;;
        --readme) cat ${SCRIPT_DIR}/README; exit 0;;
    esac
    shift
done

DB1="${DB1:-${DB1_DEFAULT}}"
DB2="${DB2:-${DB2_DEFAULT}}"
TASK1="${TASK1:-${TASK_DEFAULT}}"
TASK2="${TASK2:-${TASK_DEFAULT}}"
MINLEN="${MINLEN:-${MINLEN_DEFAULT}}"
MINID="${MINID:-${MINID_DEFAULT}}"
COMPLETE="${COMPLETE:-False}"
ADJ_DIR="${ADJ_DIR:-False}"
FEATURE="${FEATURE:-${FEATURE_DEFAULT}}"

## throw error if directory not provided
if [ -z "${DIR}" ]; then
    echo "Directory required. Please provide directory using '-d <path to directory>'"
    exit 1
## throw error if GENE or FASTA are not provided
elif [ -z "${GENE}" ] && [ -z "${FASTA}" ]; then
    echo "Gene ID(s) or an input fasta file required. Please provide gene IDs using '-g <geneIDs>' where <genes> is a comma-delimited sequence of (uppercase) gene IDs (e.g. '-g AT1G01010' for a single gene, or '-g AT1G01010,AT1G01020' for two genes, etc.). Alternatively, please provide a fasta file containing sequences to search against using '-f <path to fasta file>'"
    exit 1
## else if any combination of 2 of ACCS_F, ACC, and INPUT_DIR are provided for some reason
elif ! [ -z "${GENE}" ] && ! [ -z "${FASTA}" ]; then
    echo "Please only use either '-g <geneIDs>' or '-f <path to fasta file>'. These parameters are mutually exclusive."
    exit 1
fi

## create output directory if it doesn't exist
mkdir -p ${DIR}
DIR="$(realpath ${DIR})"
echo "Output files will be generated in ${DIR}"

## write log
script_path="$SCRIPT_DIR/${SCRIPT_NAME}"
printf -- "${params}\n\n${script_path}\n\n-g|--gene:\t${GENE}\n-f|--fasta:\t${FASTA}\n-p|--prefix:\t${PREFIX}\n-d|--dir:\t${DIR}\n--feature:\t${FEATURE}\n--complete:\t${COMPLETE}\n--task1:\t${TASK1}\n--task2:\t${TASK2}\n--db1:\t${DB1}\n--db2:\t${DB2}\n--minlen:\t${MINLEN}\n--minid:\t${MINID}\n" > ${DIR}/${PREFIX}_rblast.log

## getting sequences
if ! [ -z "${GENE}" ]; then
    fasta=${DIR}/${PREFIX}_${FEATURE}.fasta
    start_t=$(date '+%Y-%m-%d %T')
    get_seq ${DIR} ${PREFIX} ${GENE} ref ${FEATURE} ${COMPLETE} ${ADJ_DIR} ${fasta}
else
    fasta=${FASTA}
fi

## first blast (Anna Lena)
blast1_dir=${DIR}/blast1

blast1_al70_prefix=${PREFIX}.annaLena70
blast1_tsv=${blast1_dir}/${blast1_al70_prefix}.blastn_summary.tsv ## blast1 output
mkdir -p ${blast1_dir}
echo "Executing first blast"
blastn_${TASK1} ${DB1} ${blast1_dir} ${blast1_al70_prefix} ${fasta} ${blast1_tsv}
# blastn_to_al70 ${blast1_dir} ${PREFIX} ${fasta} ${blast1_tsv} ## uncomment if it's decided that this is planned to only work for Anna Lena's dataset

# blast1_prefix=${PREFIX}.blast1
# blast1_tsv=${blast1_dir}/${blast1_prefix}.blastn-${TASK1}.tsv ## blast1 output
# mkdir -p ${blast1_dir}
# echo "Executing first blast"
# run_blast6 blastn -query ${blast1_dir} ${blast1_prefix} ${fasta} ${blast1_tsv} \
#            -task ${TASK1} ${DB1}

## filter + reciprocal blast
contigs_fasta=${blast1_dir}/${blast1_al70_prefix}.minid${MINID}_minlen${MINLEN}.contigs.fasta
filter_contig ${blast1_dir} ${blast1_al70_prefix} ${blast1_tsv} ${MINLEN} ${MINID} ${contigs_fasta}

# ## filter + merge + extract sequences

# blast1_fa=${blast1_dir}/${blast1_prefix}.
# python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import merge_hits; "

## reciprocal blast
blast2_dir=${DIR}/blast2
blast2_col0_prefix=${blast1_al70_prefix}.col0nlr
blast2_tsv=${blast2_dir}/${blast2_col0_prefix}.blastn_summary.tsv ## blast2 (reciprocal blast) output
mkdir -p ${blast2_dir}
echo "Executing second blast"
blastn_${TASK2} ${DB2} ${blast2_dir} ${blast2_col0_prefix} ${contigs_fasta} ${blast2_tsv}
# blastn_to_col0nlr ${blast2_dir} ${blast1_al70_prefix} ${contigs_fasta} ## uncomment if it's decided that this is planned to only work for dataset of known Col-0 NLRs

# ## reciprocal blast
# blast2_dir=${DIR}/blast2
# blast2_prefix=${blast1_prefix}.blast2
# blast2_tsv=${blast2_dir}/${blast2_prefix}.blastn_summary.tsv ## blast2 (reciprocal blast) output
# mkdir -p ${blast2_dir}
# echo "Executing second blast"
# blastn_${TASK2} ${DB2} ${blast2_dir} ${blast2_prefix} ${contigs_fasta} ${blast2_tsv}
