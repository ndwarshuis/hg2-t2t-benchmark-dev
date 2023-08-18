import Bio.SeqIO.FastaIO as F  # type: ignore


def get_spacer(gene) -> tuple[bool | None, bool | None]:
    major = gene[:3]
    minor = gene[3]
    if major == "IGH":
        if minor == "V":
            return (None, True)
        elif minor == "D":
            return (False, False)
        elif minor == "J":
            return (True, None)
        else:
            assert False, "you dummy"
    elif major == "IGK":
        if minor == "V":
            return (None, False)
        elif minor == "J":
            return (True, None)
        else:
            assert False, "you dummy"
    elif major == "IGL":
        if minor == "V":
            return (None, True)
        elif minor == "J":
            return (False, None)
        else:
            assert False, "you dummy"
    elif major in ["TRA", "TRG"]:
        if minor == "V":
            return (None, True)
        elif minor == "J":
            return (False, None)
        else:
            assert False, "you dummy"
    elif major in ["TRB", "TRD"]:
        if minor == "V":
            return (None, True)
        elif minor == "D":
            return (False, True)
        elif minor == "J":
            return (False, None)
        else:
            assert False, "you dummy"
    else:
        assert False, "you dummy"


def get_rss(orientation, onleft) -> tuple[str, str]:
    if orientation == "+":
        # return ("GGTTTTTGT", "CACTGTG") if onleft else ("CACAGTG", "ACAAAAACC")
        return "CACTGTG" if onleft else "CACAGTG"
    else:
        # return ("ACAAAAACC", "ACAGTGT") if onleft else ("ACACTGT", "GGTTTTTGT")
        return "ACAGTGT" if onleft else "ACACTGT"


def add_seq_rss(seq, gene, orientation):
    s1, s2 = get_spacer(gene)
    if orientation == "-":
        s1, s2 = s2, s1

    def rss(spacer, onleft):
        if spacer is None:
            return ""
        else:
            return get_rss(orientation, onleft)
            # r1, r2 = get_rss(orientation, onleft)
            # return r1 + "N" * (23 if spacer else 12) + r2

    return rss(s1, True) + seq + rss(s2, False)


def main(smk):
    locus = smk.wildcards["locus"]
    loci = ["TRA", "TRD"] if locus == "TRA" else [locus]
    with open(smk.input[0], "r") as fi, open(smk.output["fa"], "w") as fo, open(
        smk.output["header"], "w"
    ) as ho:
        for r in F.FastaIterator(fi):
            hs = r.description.strip().split("|")
            gene, allele = hs[1].split("*")
            regiontype = hs[4]
            function = hs[3]
            orientation = "-" if hs[14].strip() == "rev-compl" else "+"
            this_locus = gene[:3]
            seq = str(r.seq)
            if this_locus in loci:
                ho.write("\t".join([gene, allele, function, orientation]) + "\n")
                fo.write(f">{gene}-noRSS*{allele}\n")
                fo.write(seq + "\n")
                fo.write(f">{gene}-RSS*{allele}\n")
                fo.write(add_seq_rss(seq, gene, orientation) + "\n")


main(snakemake)  # type: ignore
