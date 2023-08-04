def main(smk):
    locus = smk.wildcards["locus"]
    loci = ["TRA", "TRD"] if locus == "TRA" else [locus]
    with open(smk.input[0], "r") as fi, open(smk.output["fa"], "w") as fo, open(
        smk.output["header"], "w"
    ) as ho:
        header = next(fi, None)
        while header:
            hs = header.strip().split("|")
            gene, allele = hs[1].split("*")
            function = hs[3]
            orientation = "-" if hs[14].strip() == "rev-compl" else "+"
            this_locus = gene[:3]
            if this_locus in loci:
                fo.write(f">{gene}*{allele}\n")
                ho.write("\t".join([gene, allele, function, orientation]) + "\n")
            while seq := next(fi, None):
                if seq.startswith(">"):
                    break
                if this_locus in loci:
                    fo.write(seq)
            header = seq


main(snakemake)
