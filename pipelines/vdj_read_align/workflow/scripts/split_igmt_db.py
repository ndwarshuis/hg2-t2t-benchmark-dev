def fmt_gene(gene, regiontype):
    if (
        gene == "IGHM"
        or gene == "IGHD"
        or gene.startswith("IGHE")
        or gene.startswith("IGHPE")
        or gene.startswith("IGHA")
        or gene.startswith("IGHG")
    ):
        return f"{gene.replace('IGH', 'IGHC')}-{regiontype}"
    elif (
        gene.startswith("TRAC")
        or gene.startswith("TRBC")
        or gene.startswith("TRGC")
        or gene.startswith("TRDC")
    ):
        r = "EX2" if regiontype in ["EX2T", "EX2R"] else regiontype
        return f"{gene}-{r}"
    else:
        return gene


def main(smk):
    locus = snakemake.wildcards["locus"]
    loci = ["TRA", "TRD"] if locus == "TRA" else [locus]
    with open(smk.input[0], "r") as fi, open(smk.output["fa"], "w") as fo, open(
        smk.output["header"], "w"
    ) as ho:
        header = next(fi, None)
        while header:
            hs = header.strip().split("|")
            gene, allele = hs[1].split("*")
            function = hs[3]
            regiontype = hs[4]
            orientation = "-" if hs[14].strip() == "rev-compl" else "+"
            this_locus = gene[:3]
            if this_locus in loci:
                newgene = fmt_gene(gene, regiontype)
                fo.write(f">{newgene}*{allele}\n")
                ho.write("\t".join([newgene, allele, function, orientation]) + "\n")
            while seq := next(fi, None):
                if seq.startswith(">"):
                    break
                if this_locus in loci:
                    fo.write(seq)
            header = seq


main(snakemake)
