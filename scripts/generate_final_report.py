import re
import sys
import itertools

sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines, write_delim

## f_recip_summary: str of path to input file (first X columns are BLAST tsv output, last column is the GID of a gene that intersects with the BLAST hit, where ideally at least one of them should be in query_genes)
## fout: str of path to output file
## query_genes: list of str (gene ID)
def generate_report(f_recip_summary, query_genes, fout, col_qseqid = 0, col_bitscore=13,
                    annotated = True, generic = False, mask = [], write_mask = False, feature_type = "gene"):
    seq_type = "seq" if generic else "gene"
    query_genes = set(query_genes)
    mask = set(mask).union(query_genes)
    dat = [x.split('\t') for x in splitlines(f_recip_summary)]
    candidates = set([x[col_qseqid] for x in dat])
    dat_to_write = {}
    for candidate in candidates:
        tmp = {}
        grouped = {}
        grouped_query = {}
        for hit in [x for x in dat if x[col_qseqid]==candidate]:
            btsc = float(hit[col_bitscore])
            grouped[btsc] = grouped.get(btsc, []) + [hit]
            if hit[-1] in query_genes:
                grouped_query[btsc] = grouped_query.get(btsc, []) + [hit]
        ## get stats for best bitscore
        hits_hi = grouped[max(grouped)]
        tmp["best " + seq_type] = ','.join(sorted(set([x[-1] for x in hits_hi])))
        tmp["best bitscore"] = hits_hi[0][col_bitscore]
        if not grouped_query:
            tmp["best query rank"] = '.'
            tmp["best query " + seq_type] = '.'
            tmp["best query bitscore"] = '.'
            tmp["strict accept"] = 'F'
            tmp["mask accept"] = 'F'
            tmp["strict-mask accept"] = 'F'
        else:
            btsc_hi_query = max(grouped_query)
            hits_hi_query = grouped_query[btsc_hi_query]
            hits_higher = tuple(itertools.chain(*[hits for btsc, hits in grouped.items() if btsc > btsc_hi_query]))
            rank_hi_query = len(hits_higher) + 1
            tmp["best query rank"] = rank_hi_query
            tmp["best query " + seq_type] = ','.join(sorted(set([x[-1] for x in hits_hi_query])))
            tmp["best query bitscore"] = hits_hi_query[0][col_bitscore]
            tmp["strict accept"] = 'T' if (rank_hi_query == 1 and
                                           len(hits_hi) == len(hits_hi_query)) else 'F'
            tmp["mask accept"] = 'T' if (all(hit[-1] in mask for hit in hits_higher)) else 'F'
            tmp["strict-mask accept"] = 'T' if (tmp["mask accept"] == 'T' and
                                                len(hits_hi) == len(hits_hi_query)) else 'F'
        dat_to_write[candidate] = tmp
    if write_mask:
        stats_cols = ["best " + seq_type, "best bitscore", "best query rank", "best query " + seq_type, "best query bitscore", "strict accept", "mask accept", "strict-mask accept"]
    else:
        stats_cols = ["best " + seq_type, "best bitscore", "best query rank", "best query " + seq_type, "best query bitscore", "strict accept"]
    if annotated is None:
        pattern = "(?<=Reference\|).+?(?=\|" + feature_type + "\|)"
        get_gid = lambda candidate: (candidate if not re.search(pattern, candidate) else
                                     re.search(pattern, candidate).group(0))
    elif annotated:
        pattern = "(?<=Reference\|).+?(?=\|" + feature_type + "\|)"
        get_gid = lambda candidate: re.search(pattern, candidate).group(0)
    else:
        get_gid = lambda candidate: candidate
    to_write = [[get_gid(candidate)] +
                [dat_to_write[candidate][stats_col] for stats_col in stats_cols]
                for candidate in candidates]
    to_write = [["candidate"] + stats_cols] + to_write
    write_delim(fout, to_write, delim = '\t', coerce = True)
    return
