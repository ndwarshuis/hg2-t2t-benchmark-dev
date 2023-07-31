from tempfile import NamedTemporaryFile
import subprocess as sp
from pathlib import Path
from multiprocessing import Pool
from collections import namedtuple
from functools import partial

Paths = namedtuple(
    "Paths",
    [
        "alleles_in",
        "bam_in",
        "index_log",
        "align_log",
        "fa_out",
        "sam_out",
        "merged_out",
    ],
)

# since we have tons of bam files to merge, need to break the merge into batches
# so we don't hit the 1024 (usually) open file limit per process in linux.
BAMCHUNK = 800


def fasta_name(ps: Paths, i: int) -> str:
    return ps.fa_out / f"{i}.fasta"


def bam_name(ps: Paths, i: int) -> str:
    return ps.sam_out / f"{i}.bam"


def index_align(ps: Paths, i: int) -> None:
    fasta_path = fasta_name(ps, i)
    sam_path = bam_name(ps, i)
    with open(ps.index_log, "w") as f:
        sp.run(
            ["bwa", "index", fasta_path],
            stderr=f,
            check=True,
        )
    with open(sam_path, "wb") as fo, open(ps.align_log, "w") as fa:
        mem = sp.Popen(
            ["bwa", "mem", "-a", fasta_path, ps.alleles_in],
            stderr=fa,
            stdout=sp.PIPE,
        )
        sp.run(
            ["samtools", "view", "-b"],
            stdin=mem.stdout,
            stdout=fo,
            check=True,
        )
        mem.wait()
        if mem.returncode != 0:
            exit(1)


def chunky_list(n, xs):
    return [xs[i : i + n] for i in range(0, len(xs), n)]


def merge_bam(bams: list[Path], out):
    sp.run(["samtools", "merge", "-o", "-", *bams], check=True, stdout=out)


# def merge_sub_bam(ps: Paths, bams: list[Path]):
#     fp = tempfile.NamedTemporaryFile()
#     merge_bam(bams, fp)
#     return fp


def main(smk):
    ps = Paths(
        alleles_in=Path(smk.input["alleles"]),
        bam_in=Path(smk.input["bam"]),
        index_log=Path(smk.log["index"]),
        align_log=Path(smk.log["align"]),
        fa_out=Path(smk.output["fasta"]),
        sam_out=Path(smk.output["sam"]),
        merged_out=Path(smk.output["merged"]),
    )

    n = 0
    for p in [ps.fa_out, ps.sam_out]:
        p.mkdir(exist_ok=True, parents=True)

    with open(ps.bam_in, "r") as fi:
        while header := next(fi, None):
            with open(fasta_name(ps, n), "w") as fo:
                seq = next(fi)
                fo.write(f">{header.strip()}\n")
                fo.write(seq.strip() + "\n")
            next(fi)
            next(fi)
            n += 1

    with Pool(processes=smk.threads) as pool:
        pool.map(partial(index_align, ps), range(n))

    bams_to_merge = [bam_name(ps, i) for i in range(n)]
    with open(ps.merged_out, "wb") as f:
        if n < BAMCHUNK:
            merge_bam(bams_to_merge, f)
        else:
            chunked = [
                (c, NamedTemporaryFile()) for c in chunky_list(BAMCHUNK, bams_to_merge)
            ]
            for cs, t in chunked:
                merge_bam(cs, t)
            merge_bam([t.name for _, t in chunked], f)
            for _, t in chunked:
                t.close()


main(snakemake)  # noqa: F821
