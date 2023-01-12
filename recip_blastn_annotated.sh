#!/bin/bash

# ## test command (executed in /mnt/chaelab/rachelle/scripts/recip_blast/tmp)
# ../recip_blastn_annotated.sh --dir . --gene AT1G01010 --query-fasta /mnt/chaelab/shared/genomes/JiaoSchneeberger/C24.chr.all.v2.0.fasta --query-gff /mnt/chaelab/shared/genomes/JiaoSchneeberger/C24.protein-coding.genes.v2.5.2019-10-09.gff3

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
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

## executables
get_seq=/mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh
extract_gff_features=${DIR_SRC}/extract_gff_features.py

source /mnt/chaelab/rachelle/src/run_blast6.sh

params=${@}

# ## if no arguments, show manual ## TODO
# if [[ $# -eq 0 ]]; then
#     man -l ${SCRIPT_DIR}/MANUAL_recip_blast.1
#     exit 1
# fi

PREFIX_DEFAULT=recipblast
DB_ANNALENA=/mnt/chaelab/shared/blastdb/anna_lena70.contigs/anna_lena70.contigs
DB_COL0NLR=/mnt/chaelab/shared/blastdb/nlr164/all_NLR.164.db
FA_TAIR10=/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta
GFF_TAIR10=/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.gff
FA_ARAPORT11=/mnt/chaelab/shared/genomes/AraPort11/TAIR10_Chr.all.fasta
GFF_ARAPORT11=/mnt/chaelab/shared/genomes/AraPort11/Araport11_GFF3_genes_transposons.201606.gff
# DB1_DEFAULT="${DB_ANNALENA}"
# DB2_DEFAULT="${DB_COL0NLR}"
REFERENCE_FA_DEFAULT="${FA_TAIR10}"
REFERENCE_GFF_DEFAULT="${GFF_TAIR10}"
TASK_DEFAULT=dc-megablast
MINLEN_DEFAULT=0
MINID_DEFAULT=0
FEATURE_DEFAULT='gene'
FLANK_DEFAULT=200

while (( "$#" )); do
    case "$1" in
        -p|--prefix) PREFIX="${2}";;
        -d|--dir) DIR="${2}";;
        -g|--gene) GENE="${2}";;
        --reference-fasta) REFERENCE_FA="${2}";;
        --reference-gff) REFENCE_GFF="${2}";; ## GFF file containing information on genes in -g/--gene
        --query-fasta) QUERY_FA="${2}";; ## genome in which to find the homologues
        --query-gff) QUERY_GFF="${2}";;
        --task1) TASK1="${2}";; ## valid tasks: blastn, blastn-short, dc-megablast, megablast, rmblastn
        --task2) TASK2="${2}";;
        --minlen) MINLEN="${2}";;
        --minid) MINID="${2}";;
        --flank) FLANK="${2}";;
        # --complete) COMPLETE='True';;
        # --feature) FEATURE="${2}";;
        # --adjust-dir) ADJ_DIR='True';;
        # -h|--help) man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1; exit 0;;
        # --readme) cat ${SCRIPT_DIR}/README; exit 0;;
    esac
    shift
done

# DB1="${DB1:-${DB1_DEFAULT}}"
# DB2="${DB2:-${DB2_DEFAULT}}"
PREFIX="${PREFIX:=${PREFIX_DEFAULT}}"
REFERENCE_FA="${REFERENCE_FA:-${REFERENCE_FA_DEFAULT}}"
REFERENCE_GFF="${REFERENCE_GFF:-${REFERENCE_GFF_DEFAULT}}"
TASK1="${TASK1:-${TASK_DEFAULT}}"
TASK2="${TASK2:-${TASK_DEFAULT}}"
MINLEN="${MINLEN:-${MINLEN_DEFAULT}}"
MINID="${MINID:-${MINID_DEFAULT}}"
# COMPLETE="${COMPLETE:-False}"
# ADJ_DIR="${ADJ_DIR:-False}"
# FEATURE="${FEATURE:-${FEATURE_DEFAULT}}"
FLANK=${FLANK:-${FLANK_DEFAULT}}


## throw error if directory not provided
if [ -z "${DIR}" ]; then
    echo "Directory required. Please provide directory using '-d <path to directory>'"
    exit 1
## throw error if GENE not provided
elif [ -z "${GENE}" ]; then
    echo "Gene ID(s) required. Please provide gene ID(s) using '-g <geneIDs>' where <genesIDs> is a comma-delimited sequence of case-sensitive gene IDs (e.g. '-g AT1G01010' for a single gene, or '-g AT1G01010,AT1G01020' for two genes, etc.)."
## throw error if QUERY_FA or QUERY_GFF is not provided
elif [ -z "${QUERY_FA}" ] || [ -z "${QUERY_GFF}" ]; then
    echo "Query FASTA and GFF are required. Use '--query-fasta <path to fasta file>' and '--query-gff <path to gff file>'."
    exit 1
fi

## create output directory if it doesn't exist
mkdir -p ${DIR}
DIR="$(realpath ${DIR})"
echo "Output files will be generated in ${DIR}"

## write log
printf -- "${params}\n\n-p|--prefix:\t${PREFIX}\n-d|--dir:\t${DIR}\n-g|--gene:\t${GENE}\n--reference-fasta:\t${REFERENCE_FA}\n--reference-gff:\t${REFERENCE_GFF}\n--query-fasta:\t${QUERY_FA}\n--query-gff:\t${QUERY_GFF}\n--task1:\t${TASK1}\n--task2:\t${TASK2}\n--minlen:\t${MINLEN}\n--minid:\t${MINID}\n" > ${DIR}/${PREFIX}_rblastn.log

## getting annotations and sequences
dir_qgenes=${DIR}/query_genes
mkdir -p ${dir_qgenes}
genes_gid=${dir_qgenes}/${PREFIX}_query_gene.gid
genes_gff=${dir_qgenes}/${PREFIX}_query_gene.gff
genes_fasta=${dir_qgenes}/${PREFIX}_query_gene.fasta
## get gene IDs into file
echo ${GENE} | tr ',' '\n' > ${genes_gid}
## subset GFF annotations for these genes
echo "Filtering reference annotations for query genes"
${extract_gff_features} ${REFERENCE_GFF} ${genes_gid} ${genes_gff} GFF
## get gene sequences
echo "Extracting query gene sequences"
start_t=$(date '+%Y-%m-%d %T')
${get_seq} --gene ${GENE} --acc ref --feature gene \
           --complete --adjust-dir \
           --reference ${REFERENCE_FA} --gff ${genes_gff} \
           --dir ${dir_qgenes} --out ${genes_fasta} --no-bed

## first blast (query genes to query fasta)
dir_blast1=${DIR}/blast1
mkdir -p ${dir_blast1}
tsv_blast1=${dir_blast1}/${PREFIX}.blast1.blastn-${TASK1}.tsv
bed_blast1=${tsv_blast1%*.tsv}.merged.bed
gff_intersect=${tsv_blast1%.tsv}.intersect_query.gene.gff
gid_intersect=${tsv_blast1%.tsv}.intersect_query.gid
fa_intersect=${tsv_blast1%.tsv}.intersect_query.gene.fasta
## run blast
echo "Executing first BLASTN of query genes against query genome"
run_blast6 blastn -query ${genes_fasta} -subject ${QUERY_FA} \
           -outfmt "${OUTFMT}" \
           -out ${tsv_blast1} -task ${TASK1}
## append newline because the output doesn't end with newline and this might throw problems in some programmes
printf "\n" >> ${tsv_blast1}
## get hit ranges
echo "Extract BLAST hit ranges"
awk -v OFS='\t' -v minid=${MINID} -v minlen=${MINLEN} -v sseqid=${COL_SSEQID} -v start=${COL_START} -v end=${COL_END} -v pident=${COL_PIDENT} '{if (NR>1 && $(pident)>=minid) { if ($(start)<$(end) && $(end)-$(start)+1>=minlen) { print $(sseqid),$(start)-1,$(end) } else if ($(start)-$(end)+1>=minlen) { print $(sseqid),$(end)-1,$(start) }}}' ${tsv_blast1} | bedtools sort -i stdin | bedtools merge -i stdin > ${bed_blast1}
## get gff of intersecting genes (feature type = gene)
echo "Finding genes in query genome that intersect with hit ranges"
bedtools intersect -wb -a ${bed_blast1} -b ${QUERY_GFF} | cut -f4- | awk '$3=="gene"' > ${gff_intersect}
## get gid of intersecting genes
cut -f9 ${gff_intersect} | grep -Po '(^ID=|;ID=)\K[^;]+' > ${gid_intersect}
## get sequences of intersecting genes
echo "Extracting sequences of intersecting genes"
${get_seq} --gene $(tr '\n' ',' < ${gid_intersect} | sed 's/,\+/,/g' | sed 's/,$//') \
           --acc ref --feature gene \
           --complete --adjust-dir \
           --reference ${QUERY_FA} --gff ${gff_intersect} \
           --dir ${dir_blast1} --out ${fa_intersect} --no-bed

## second blast (intersected genes to reference fasta)
dir_blast2=${DIR}/blast2
mkdir -p ${dir_blast2}
tsv_blast2=${dir_blast2}/${PREFIX}.blast2.blastn-${TASK2}.tsv
bed_intersect_recip=${tsv_blast2%.tsv}.intersect_ref.bed
sum_recip=${tsv_blast2%.tsv}.intersect_ref.summary.txt
## run blast on reference genome
echo "Executing reciprocal BLAST"
run_blast6 blastn -query ${fa_intersect} -subject ${REFERENCE_FA} \
           -outfmt "${OUTFMT}" \
           -out ${tsv_blast2} -task ${TASK2}
## append newline because the output doesn't end with newline and this might throw problems in some programmes
printf "\n" >> ${tsv_blast2}
## intersect with reference annotations
echo "Getting reciprocal annotations"
awk -v OFS='\t' -v sseqid=${COL_SSEQID} -v start=${COL_START} -v end=${COL_END} '{if (NR>1) { if ($(start)<$(end)) { print $(sseqid),$(start)-1,$(end),$0 } else { print $(sseqid),$(end)-1,$(start),$0 }}}' ${tsv_blast2} | bedtools intersect -loj -a stdin -b ${REFERENCE_GFF} > ${bed_intersect_recip}
## retain only intersections with genes or intersections with nothing
## write file with columns from BLAST output, and append column with intersected reference gene ID
echo "Getting reference genes from reciprocal annotations"
## extracting blast columns and GFF attributes column
awk -v OFS='\t' -v blastcols=${OUTFMT6_NCOL} \
    '{ if ($(blastcols+3+3)=="gene" || $(blastcols+3+3)==".") {print $0} }' ${bed_intersect_recip} |
    cut -f4-$((${OUTFMT6_NCOL}+3)),$((${OUTFMT6_NCOL}+3+9)) > ${sum_recip}
## extract gene ID from GFF attributes column, and replace GFF attributes column with gene IDs
python3 -c "import re; import sys; sys.path.append('${DIR_SRC}'); from data_manip import splitlines, write_delim; dat = [x.split('\t') for x in splitlines('${sum_recip}')]; dat = [x[:-1] + ['.' if x[-1] == '.' else re.search('(?<=ID=)[^;]+', re.search('(^|;)ID=[^;]+', x[-1]).group(0)).group(0)] for x in dat]; write_delim('${sum_recip}', dat, delim = '\t', coerce = True)"

## summarise recip blast run
dir_final=${DIR}/final
dir_gb=${dir_final}/genbank
mkdir -p ${dir_final} ${dir_gb}
txt_report=${dir_final}/${PREFIX}.report.txt
final_genes_gid=${dir_final}/${PREFIX}.final.gid
gff_final=${dir_final}/${PREFIX}.final.gff
fa_final_gene=${dir_final}/${PREFIX}.final.gene.fasta
fa_final_cds=${dir_final}/${PREFIX}.final.CDS.fasta
## writing report
echo "Writing report"
python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from generate_final_report import generate_report; generate_report('${sum_recip}', '${GENE}'.split(','), '${txt_report}', col_qseqid=${COL_QSEQID}-1, col_bitscore=${COL_BITSCORE}-1)"
## filter query GFF for annotations of best candidates (i.e. candidates which top bitscore is to a query gene)
echo "Extracting IDs of best candidate homologues"
awk 'NR>1 {if($4==1) {print $1}}' ${txt_report} > ${final_genes_gid}
echo "Making genbank files for best candidate homologues"
python3 -c "import sys; sys.path.append('${DIR_SRC}'); from extract_seq_and_annotate import get_flanked_seq; from reference import AnnotatedFasta; from data_manip import splitlines; ref = AnnotatedFasta('${QUERY_FA}', '${QUERY_GFF}', name = 'query'); get_flanked_seq(ref, splitlines('${final_genes_gid}'), flank = ${FLANK}, feature_type = 'gene', prefix = '${PREFIX}.final', gff_subset = '${gff_final}', gb_dir = '${dir_gb}')"
## get sequences
echo "Extracting sequence for best candidates homologues"
${get_seq} --gene $(tr '\n' ',' < ${final_genes_gid} | sed 's/,\+/,/g' | sed 's/,$//') \
           --acc ref --feature gene \
           --complete --adjust-dir \
           --reference ${QUERY_FA} --gff ${gff_final} \
           --dir ${dir_final} --out ${fa_final_gene} --no-bed
${get_seq} --gene $(tr '\n' ',' < ${final_genes_gid} | sed 's/,\+/,/g' | sed 's/,$//') \
           --acc ref --feature CDS \
           --adjust-dir \
           --reference ${QUERY_FA} --gff ${gff_final} \
           --dir ${dir_final} --out ${fa_final_cds} --no-bed

