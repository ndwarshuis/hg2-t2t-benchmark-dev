import Bio.SeqIO.FastaIO as F  # type: ignore
import pylcs  # type: ignore
import pandas as pd  # type: ignore
from multiprocessing import Pool


with open(snakemake.input[0], "r") as f:
    d = {r.name: str(r.seq) for r in F.FastaIterator(f)}


def edit_distance(s1, s2):
    return len(s1) + len(s2) - 2 * pylcs.lcs_string_length(s1, s2)


def mk_column(s1):
    return [edit_distance(s1, s2) for s2 in d.values()]


with Pool(snakemake.threads) as p:
    df = pd.DataFrame(
        list(p.map(mk_column, d.values())),
        columns=list(d.keys()),
        index=list(d.keys()),
    )

df.reset_index(names=["gene"]).to_csv(snakemake.output[0], sep="\t", index=False)
