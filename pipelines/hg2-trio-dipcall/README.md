# hg2 t2t parental alignments

this pipeline contains dipcall alignments for four asm combinations
* hg3 diploid -> hg2 paternal haplotype ref
* hg4 diploid -> hg2 maternal haplotype ref
* hg2 diploid -> hg2 paternal haplotype ref
* hg2 diploid -> hg2 maternal haplotype ref

In all cases, the diploid asms are from Hifiasm (see snakemake source for links)
and the ref haplotypes are from verkko (ditto for links).

The point of the hg3/4 -> hg2 alignments is for polishing the hg2 t2t asm, and
the hg2 -> hg2 alignments are for finding errors in hifiasm (which may
indirectly help with polishing the hg2 t2t asm).

## output

These were uploaded to `s3://giab-data/hg2-t2t/parental_alignments/<run_date>/`.
In each subpath are four directories (corresponding to the combinations listed
above) which each contain dipcall results. Note that in the case of the final
vcf, the final vcf is called `hg2.fixed.vcf.gz` rather than `hg2.dip.vcf.gz`
(for technical reasons explained in the pipeline source). The "fixed" vcf and
pair vcf also have tbi files for indexing and easy import into IGV.

