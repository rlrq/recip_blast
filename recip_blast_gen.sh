#!/bin/bash

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
SCRIPT_NAME="$(basename "$0")"
DIR_SRC=/mnt/chaelab/rachelle/src

OUTFMT6_COLS=( qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore )
OUTFMT6_NCOL=${#OUTFMT6_COLS[@]}
COL_QSEQID=1
COL_SSEQID=2
COL_START=9
COL_END=10
COL_PIDENT=3
COL_BITSCORE=14
OUTFMT="6 ${OUTFMT6_COLS[@]}"

source /mnt/chaelab/rachelle/src/run_blast6.sh
source "${SCRIPT_DIR}/scripts/reciprocal_blastn_annaLena70.sh"

## if no arguments, show manual ## TODO
if [[ $# -eq 0 ]]; then
    man -l ${SCRIPT_DIR}/MANUAL_recip_blast.1
    exit 1
fi

PREFIX_DEFAULT=recipblast
BG_FASTA_DEFAULT=/mnt/chaelab/shared/genomes/TAIR10/fasta/TAIR10_pep_20101214.fasta
FA_TAIR10=/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta
GFF_TAIR10=/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.gff
FA_ARAPORT11=/mnt/chaelab/shared/genomes/AraPort11/TAIR10_Chr.all.fasta
GFF_ARAPORT11=/mnt/chaelab/shared/genomes/AraPort11/Araport11_GFF3_genes_transposons.201606.gff
REFERENCE_FA_DEFAULT="${FA_TAIR10}"
REFERENCE_GFF_DEFAULT="${GFF_TAIR10}"
DB_ANNALENA='/mnt/chaelab/shared/blastdb/anna_lena70.contigs/anna_lena70.contigs'
DB_COL0NLR='/mnt/chaelab/shared/blastdb/nlr164/all_NLR.164.db'
# DB1_DEFAULT="${DB_ANNALENA}"
# DB2_DEFAULT="${DB_COL0NLR}"
TASK_DEFAULT='default'
TASK_BLASTN_DEFAULT=dc-megablast
TASK_BLASTP_DEFAULT=blastp
MINLEN_DEFAULT=0
MINID_DEFAULT=0
# FEATURE_DEFAULT='gene'
SEQ_TYPE_DEFAULT='n'
RELAX_DEFAULT=1

## executables
get_seq=/mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh
extract_gff_features=${DIR_SRC}/extract_gff_features.py

params=${@}

while (( "$#" )); do
    case "$1" in
        -g|--gene|-s|--seqid) GENE="${2}";; ## seqid of query sequences (within --fasta input file)
        -f|--fasta) FASTA="${2}";; ## should contain all sequences for query
        -b|--bg-fasta) BG_FASTA="${2}";; ## should contain all sequences for background
        -p|--prefix) PREFIX="${2}";;
        -d|--dir) DIR="${2}";;
        -m|--mask) MASK_ID="${2}";; ## IDs of sequences in BG_FASTA to ignore (will be included in report but hits to it will be treated as 0 bitscore when it comes to the 'mask accept' and 'strict-mask accept' column [new column!] of a potential homologue)
        --seq-type) SEQ_TYPE="${2}";; ## 'n' or 'p' (where n = nucleotide and p = peptide)
        --reference-fasta) REFERENCE_FA="${2}";;
        --reference-gff) REFERENCE_GFF="${2}";; ## GFF file containing information on genes in -g/--gene
        --query-fasta) QUERY_FA="${2}";; ## fasta file in which to find the homologues
        --task1) TASK1="${2}";; ## valid tasks:'blastn', 'default' (megablast), 'dc' (discontiguous megablast)
        --task2) TASK2="${2}";;
        --topn) TOPN="${2}";; ## top N (by bitscore) number of blast1 hits per query to consider for blast2
        # --db1) DB1="${2}";;
        # --db2) DB2="${2}";;
        --minlen) MINLEN="${2}";;
        --minid) MINID="${2}";;
        # --complete) COMPLETE=1;;
        # --feature) FEATURE="${2}";;
        --relax) RELAX=1;; ## should be raised if query sequences are present in --bg-fasta but with different IDs
        --strict) RELAX=0;; ## whichever comes last (--relax or --strict) in the command will be used
        # --adjust-dir) ADJ_DIR='True';;
        # -h|--help) man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1; exit 0;;
        --readme) cat ${SCRIPT_DIR}/README; exit 0;;
    esac
    shift
done

SEQ_TYPE="${SEQ_TYPE:-${SEQ_TYPE_DEFAULT}}"
if [ "${SEQ_TYPE}" == 'p' ]; then
    echo "Sequence type is peptide (--seq-type p). Setting '--feature CDS' and ignoring '--complete' flag."
    # FEATURE=CDS
    # COMPLETE=0
    TASK1="${TASK1:-${TASK_BLASTP_DEFAULT}}"
    TASK2="${TASK2:-${TASK_BLASTP_DEFAULT}}"
    BLAST_PROGRAMME=blastp
else
    # FEATURE="${FEATURE:-${FEATURE_DEFAULT}}"
    # COMPLETE="${COMPLETE:-0}"
    TASK1="${TASK1:-${TASK_BLASTN_DEFAULT}}"
    TASK2="${TASK2:-${TASK_BLASTN_DEFAULT}}"
    BLAST_PROGRAMME=blastp
fi

BG_FASTA="${BG_FASTA:-${BG_FASTA_DEFAULT}}"
PREFIX="${PREFIX:-${PREFIX_DEFAULT}}"
REFERENCE_FA="${REFERENCE_FA:-${REFERENCE_FA_DEFAULT}}"
REFERENCE_GFF="${REFERENCE_GFF:-${REFERENCE_GFF_DEFAULT}}"
# DB1="${DB1:-${DB1_DEFAULT}}"
# DB2="${DB2:-${DB2_DEFAULT}}"
MINLEN="${MINLEN:-${MINLEN_DEFAULT}}"
MINID="${MINID:-${MINID_DEFAULT}}"
# ADJ_DIR="${ADJ_DIR:-False}"


# if [ "${SEQ_TYPE}" == 'p' ]; then
#     echo "Sequence type is peptide (--seq-type p). Setting '--feature CDS' and ignoring '--complete' flag"
#     FEATURE=CDS
#     COMPLETE=0
# fi

## throw error if directory not provided
if [ -z "${DIR}" ]; then
    echo "Directory required. Please provide directory using '-d <path to directory>'"
    exit 1
## throw error if FASTA not provided
elif [ -z "${BG_FASTA}" ]; then
    echo "Input fasta file of background sequences required. Use '--bg-fasta <path to fasta file>'."
    exit 1
# ## throw error if GENE not provided
# elif [ -z "${GENE}" ]; then
#     echo "Sequence ID(s) required. Please provide sequence IDs using '-g <seqids>' or '-s <seqids>', where <seqids> is a comma-delimited sequene of sequence IDs (as used in --fasta file)."
#     exit 1
## throw error if GENE or FASTA are not provided
elif [ -z "${GENE}" ] && [ -z "'${FASTA}" ]; then
    echo "Gene ID(s) or an input fasta file required. Please provide gene IDs using '-g <geneIDs>' where <genes> is a comma-delimited sequence of (case-sensitive) gene IDs (e.g. '-g AT1G01010' for a single gene, or '-g AT1G01010,AT1G01020' for two genes, etc.). Alternatively, please provide a fasta file containing sequences to search against using '-f <path to fasta file>'"
    exit 1
## else if any combination of 2 of ACCS_F, ACC, and INPUT_DIR are provided for some reason
elif ! [ -z "${GENE}" ] && ! [ -z "${FASTA}" ]; then
    echo "Please only use either '-g <geneIDs>' or '-f <path to fasta file>'. These parameters are mutually exclusive."
    exit 1
## throw error if QUERY_FA is not provided
elif [ -z "${QUERY_FA}" ]; then
    echo "Query FASTA is required. Use '--query-fasta <path to fasta file>'."
    exit 1
fi


## create output directory if it doesn't exist
mkdir -p ${DIR}
DIR="$(realpath ${DIR})"
echo "Output files will be generated in ${DIR}"

## write log
script_path="$SCRIPT_DIR/${SCRIPT_NAME}"
printf -- "${params}\n\n${script_path}\n\n-g|--gene|-s|--seqid:\t${GENE}\n--query-fasta:\t${QUERY_FA}\n-f|--fasta:\t${FASTA}\n-b|--bg-fasta:\t${BG_FASTA}\n-m|--mask:\t${MASK_ID}\n--seq-type:\t${SEQ_TYPE}\n-p|--prefix:\t${PREFIX}\n-d|--dir:\t${DIR}\n--feature:\t${FEATURE}\n--complete:\t${COMPLETE}\n--task1:\t${TASK1}\n--task2:\t${TASK2}\n--db1:\t${DB1}\n--db2:\t${DB2}\n--minlen:\t${MINLEN}\n--minid:\t${MINID}\n" > ${DIR}/${PREFIX}_rblast.log

# ## getting sequences
# if ! [ -z "${GENE}" ]; then
#     fasta=${DIR}/${PREFIX}_${FEATURE}.fasta
#     start_t=$(date '+%Y-%m-%d %T')
#     get_seq ${DIR} ${PREFIX} ${GENE} ref ${FEATURE} ${COMPLETE} ${ADJ_DIR} ${fasta}
# else
#     fasta=${FASTA}
# fi

## getting annotations and sequences
dir_qgenes=${DIR}/query_genes
genes_gid=${dir_qgenes}/${PREFIX}_query.gid
mkdir -p ${dir_qgenes}
if ! [ -z "${GENE}" ]; then
    genes_gff=${dir_qgenes}/${PREFIX}_query_gene.gff
    if [ "${SEQ_TYPE}" == 'p' ]; then
        genes_fasta=${dir_qgenes}/${PREFIX}_query_protein.fasta
    else
        genes_fasta=${dir_qgenes}/${PREFIX}_query_${FEATURE}.fasta
    fi
    python3 -c "seqids = set('${GENE}'.split(',')); from Bio import SeqIO; to_keep = [r for r in SeqIO.parse('${BG_FASTA}', 'fasta') if r.id in seqids]; SeqIO.write(to_keep, '${genes_fasta}', 'fasta')"
    # ## get gene IDs into file
    # echo ${GENE} | tr ',' '\n' > ${genes_gid}
    # ## subset GFF annotations for these genes
    # echo "Filtering reference annotations for query genes"
    # ${extract_gff_features} ${REFERENCE_GFF} ${genes_gid} ${genes_gff} GFF
    # ## if no annotations, abort
    # if [ $(grep -vP '^$' ${genes_gff} | wc -l) -eq 0 ]; then
    #     echo "No annotations found (${GENE}). Are you sure the gene IDs and/or GFF file are correct?"
    #     exit 2
    # fi
    # ## get gene sequences
    # echo "Extracting query gene sequences"
    # start_t=$(date '+%Y-%m-%d %T')
    # if [ "${COMPLETE}" == '1' ]; then
    #     ${get_seq} --gene ${GENE} --acc ref --feature ${FEATURE} \
    #                --complete --adjust-dir \
    #                --reference ${REFERENCE_FA} --gff ${genes_gff} \
    #                --dir ${dir_qgenes} --out ${genes_fasta} --no-bed
    # else
    #     ${get_seq} --gene ${GENE} --acc ref --feature ${FEATURE} \
    #                --adjust-dir \
    #                --reference ${REFERENCE_FA} --gff ${genes_gff} \
    #                --dir ${dir_qgenes} --out ${genes_fasta} --no-bed
    # fi
    # ## translate (if necessary)
    # if [ "${SEQ_TYPE}" == 'p' ]; then
    #     /mnt/chaelab/rachelle/src/translate.py ${genes_fasta} ${genes_fasta}
    # fi
    # ## if no sequence entries, abort
    # if [ $(grep '>' ${genes_fasta} | wc -l) -eq 0 ]; then
    #     echo "Unable to extract gene sequences (${GENE}). Are you sure the gene IDs and/or GFF file are correct?"
    #     exit 2
    # fi
else
    genes_fasta=${dir_qgenes}/${PREFIX}_query.fasta
    grep '>' ${FASTA} | sed 's/^>//' | tr ' ' '\t' | cut -f1 > ${genes_gid}
    GENE=$( cat ${genes_gid} | tr '\n' ',' | sed 's/,$//' )
    cp ${FASTA} ${genes_fasta}
fi


## first blast (query genes to query fasta)
dir_blast1=${DIR}/blast1
mkdir -p ${dir_blast1}
tsv_blast1=${dir_blast1}/${PREFIX}.blast1.${BLAST_PROGRAMME}-${TASK1}.tsv
# bed_blast1=${tsv_blast1%*.tsv}.merged.bed
gid_blast1=${tsv_blast1%.tsv}.query.gid
fa_blast1=${tsv_blast1%.tsv}.query.fasta
## run blast
echo "Executing first BLAST of query genes against query sequences"
run_blast6 ${BLAST_PROGRAMME} -query ${genes_fasta} -subject ${QUERY_FA} \
           -outfmt "${OUTFMT}" \
           -out ${tsv_blast1} -task ${TASK1}
## if no hits, abort (there is a header row so an empty file has 1 non-empty line)
if [ $(grep -vP '^$' ${tsv_blast1} | wc -l) -eq 1 ]; then
    echo "No BLAST hits found in query genome. Exiting."
    exit 3
fi
## append newline because the output doesn't end with newline and this might throw problems in some programmes
printf "\n" >> ${tsv_blast1}
## filter (in-place) by bitscore if TOPN is provided
if [[ ! -z ${TOPN} ]]; then
    python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); import data_manip; data = [x.split('\t') for x in data_manip.splitlines('${tsv_blast1}')]; header = data[0]; data = data[1:]; i_qacc = header.index('qseqid'); i_bitscore = header.index('bitscore'); qaccs = set(x[i_qacc] for x in data); bitscore_cutoff_by_qacc = {qacc: sorted([float(x[i_bitscore]) for x in data if x[i_qacc] == qacc], reverse = True)[${TOPN}-1] for qacc in qaccs}; data_filtered = [x for x in data if float(x[i_bitscore]) >= bitscore_cutoff_by_qacc[x[i_qacc]]]; data_manip.write_delim('${tsv_blast1}', [header] + data_filtered, delim = '\t', coerce = True)"
fi
## get gid of hits
awk -v OFS='\t' -v sseqid=${COL_SSEQID} 'NR>1 { print $(sseqid) }' ${tsv_blast1} | sort | uniq > ${gid_blast1}
## get seq of hits
python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from file_manip import splitlines; gids = set(splitlines('${gid_blast1}')); from Bio import SeqIO; seqs_to_write = [record for record in SeqIO.parse('${QUERY_FA}', 'fasta') if record.id in gids]; SeqIO.write(seqs_to_write, '${fa_blast1}', 'fasta')"


## second blast (query genes to input fasta)
dir_blast2=${DIR}/blast2
mkdir -p ${dir_blast2}
tsv_blast2=${dir_blast2}/${PREFIX}.blast2.${BLAST_PROGRAMME}-${TASK2}.tsv
sum_recip=${tsv_blast2%.tsv}.intersect_ref.summary.txt
## run blast on background sequences
echo "Executing reciprocal BLAST"
run_blast6 ${BLAST_PROGRAMME} -query ${fa_blast1} -subject ${BG_FASTA} \
           -outfmt "${OUTFMT}" \
           -out ${tsv_blast2} -task ${TASK2}
## run blast on query sequences
tmp=${DIR}/tmp.tsv
${BLAST_PROGRAMME} -query ${fa_blast1} -subject ${genes_fasta} \
                   -outfmt "${OUTFMT}" \
                   -out ${tmp} -task ${TASK2}
## append newline because the output doesn't end with newline and this might throw problems in some programmes
printf "\n" >> ${tsv_blast2}
printf "\n" >> ${tmp}
cat ${tmp} >> ${tsv_blast2}
rm ${tmp}
## summarise (well I say summarise but it's really just blast columns + repeat the sseqid column at end, and remove header
python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from data_manip import splitlines, write_delim; data = [x.split('\t') for x in splitlines('${tsv_blast2}')[1:]]; data = [x + [x[${COL_SSEQID}-1]] for x in data]; write_delim('${sum_recip}', data, delim = '\t', coerce = True)"


## summarise recip blast run
dir_final=${DIR}/final
mkdir -p ${dir_final}
txt_report=${dir_final}/${PREFIX}.report.txt
final_genes_gid=${dir_final}/${PREFIX}.final.gid
fa_final=${dir_final}/${PREFIX}.final.fasta
## writing report
echo "Writing report"
python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from generate_final_report import generate_report; generate_report('${sum_recip}', '${GENE}'.split(','), '${txt_report}', col_qseqid=${COL_QSEQID}-1, col_bitscore=${COL_BITSCORE}-1, annotated=None, write_mask = True, mask = ([] if '${MASK_ID}' == '' else '${MASK_ID}'.split(',')), generic = True, feature_type = '${FEATURE}')"
## filter query GFF for annotations of best candidates (i.e. candidates which top bitscore is to a query gene)
echo "Extracting IDs of best candidate homologues"
## if not relax (i.e. reject if highest bitscore to target is equal to non-target)
if [ "${RELAX}" == '0'  ]; then
    echo "--strict mode--"
    awk 'NR>1 {if($7=="T") {print $1}}' ${txt_report} > ${final_genes_gid}
else
    echo "--relax mode--"
    awk 'NR>1 {if($4==1) {print $1}}' ${txt_report} > ${final_genes_gid}
fi
## if no homologues, abort
if [ $(grep -vP '^$' ${final_genes_gid} | wc -l) -eq 0 ]; then
    if [ "${RELAX}" == '0' ]; then
        echo "No homologues found. Try using --relax."
    else
        echo "No homologues found."
    fi
    exit 3
fi
## get sequences
echo "Extracting sequence for best candidates homologues"
python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from file_manip import splitlines; gids = set(splitlines('${final_genes_gid}')); from Bio import SeqIO; seqs_to_write = [record for record in SeqIO.parse('${QUERY_FA}', 'fasta') if record.id in gids]; SeqIO.write(seqs_to_write, '${fa_final}', 'fasta')"

exit 0
