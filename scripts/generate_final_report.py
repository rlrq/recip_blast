import re
import sys

sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import splitlines, write_delim

## f_recip_summary: str of path to input file (first X columns are BLAST tsv output, last column is the GID of a gene that intersects with the BLAST hit, where ideally at least one of them should be in query_genes)
## fout: str of path to output file
## query_genes: list of str (gene ID)
def generate_report(f_recip_summary, query_genes, fout, col_qseqid = 0, col_bitscore=13):
    query_genes = set(query_genes)
    dat = [x.split('\t') for x in splitlines(f_recip_summary)]
    candidates = set([x[col_qseqid] for x in dat])
    dat_to_write = {}
    for candidate in candidates:
        tmp = {}
        ranked = sorted([x for x in dat if x[col_qseqid]==candidate],
                             key = lambda y:float(y[col_bitscore]),
                             reverse = True)
        tmp["best gene"] = ranked[0][-1]
        tmp["best bitscore"] = ranked[0][col_bitscore]
        query_hits_i = [i for i,x in enumerate(ranked) if x[-1] in query_genes]
        if not query_hits_i:
            tmp["best query rank"] = '.'
            tmp["best query gene"] = '.'
            tmp["best query bitscore"] = '.'
        else:
            best_query_hit = ranked[query_hits_i[0]]
            tmp["best query rank"] = query_hits_i[0]+1
            tmp["best query gene"] = best_query_hit[-1]
            tmp["best query bitscore"] = best_query_hit[col_bitscore]
        dat_to_write[candidate] = tmp
    stats_cols = ["best gene", "best bitscore", "best query rank", "best query gene", "best query bitscore"]
    to_write = [[re.search("(?<=Reference\|).+?(?=\|gene\|)", candidate).group(0)] +
                [dat_to_write[candidate][stats_col] for stats_col in stats_cols]
                for candidate in candidates]
    to_write = [["candidate"] + stats_cols] + to_write
    write_delim(fout, to_write, delim = '\t', coerce = True)
    return
