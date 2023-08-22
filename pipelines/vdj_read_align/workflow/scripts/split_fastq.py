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
        # "motifs_in",
        "bam_in",
        "allele_index_log",
        # "motif_index_log",
        "allele_align_log",
        # "motif_align_log",
        "fa_out",
        "sam_out",
    ],
)

# since we have tons of bam files to merge, need to break the merge into batches
# so we don't hit the 1024 (usually) open file limit per process in linux.
BAMCHUNK = 800


def fasta_name(ps: Paths, i: int) -> str:
    return ps.fa_out / f"{i}.fasta"


def bam_name(ps: Paths, i: int, rss: bool) -> str:
    return ps.sam_out / f"{i}_{'rss' if rss else 'allele'}.bam"


def index_align(ps: Paths, i: int) -> None:
    fasta_path = fasta_name(ps, i)
    allele_path = bam_name(ps, i, False)
    rss_path = bam_name(ps, i, True)
    with open(ps.allele_index_log, "w") as f:
        sp.run(
            ["bwa", "index", fasta_path],
            stderr=f,
            check=True,
        )
    # with open(ps.motif_index_log, "w") as f:
    #     sp.run(
    #         ["bowtie-build", "-q", fasta_path, fasta_path],
    #         stderr=f,
    #         check=True,
    # )
    with open(allele_path, "wb") as fo, open(ps.allele_align_log, "w") as fa:
        mem = sp.Popen(
            ["bwa", "mem", "-k", "14", "-a", "-T", "10", fasta_path, ps.alleles_in],
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
    # with open(rss_path, "wb") as fo, open(ps.motif_align_log, "w") as fa:
    #     bowtie = sp.Popen(
    #         ["bowtie", "-a", "-x", fasta_path, "-f", ps.motifs_in, "-S"],
    #         stderr=fa,
    #         stdout=sp.PIPE,
    #     )
    #     sp.run(
    #         ["samtools", "view", "-b"],
    #         stdin=bowtie.stdout,
    #         stdout=fo,
    #         check=True,
    #     )
    #     mem.wait()
    #     if mem.returncode != 0:
    #         exit(1)


def chunky_list(n, xs):
    return [xs[i : i + n] for i in range(0, len(xs), n)]


def merge_bam(bams: list[Path], out):
    sp.run(["samtools", "merge", "-o", "-", *bams], check=True, stdout=out)


def main(smk):
    ps = Paths(
        alleles_in=Path(smk.input["alleles"]),
        # motifs_in=Path(smk.input["motifs"]),
        bam_in=Path(smk.input["bam"]),
        allele_index_log=Path(smk.log["allele_index"]),
        # motif_index_log=Path(smk.log["motif_index"]),
        allele_align_log=Path(smk.log["allele_align"]),
        # motif_align_log=Path(smk.log["motif_align"]),
        fa_out=Path(smk.output["fasta"]),
        sam_out=Path(smk.output["sam"]),
    )
    alleles_merged_out = Path(smk.output["allele_merged"])
    # motifs_merged_out = Path(smk.output["motif_merged"])

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

    def merge_output(path, ns):
        with open(path, "wb") as f:
            if n < BAMCHUNK:
                merge_bam(ns, f)
            else:
                chunked = [(c, NamedTemporaryFile()) for c in chunky_list(BAMCHUNK, ns)]
                for cs, t in chunked:
                    merge_bam(cs, t)
                merge_bam([t.name for _, t in chunked], f)
                for _, t in chunked:
                    t.close()

    merge_output(alleles_merged_out, [bam_name(ps, i, False) for i in range(n)])
    # merge_output(motifs_merged_out, [bam_name(ps, i, True) for i in range(n)])


main(snakemake)  # noqa: F821
