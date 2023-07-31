from pathlib import Path


def main(smk):
    out = Path(smk.output[0])
    out.mkdir(exist_ok=True, parents=True)

    with open(smk.input[0], "r") as fi:
        n = 0
        while header := next(fi, None):
            with open(out / f"{n}.fasta", "w") as fo:
                seq = next(fi)
                fo.write(f">{header.strip()}\n")
                fo.write(seq.strip() + "\n")
            next(fi)
            next(fi)
            n += 1


main(snakemake)
