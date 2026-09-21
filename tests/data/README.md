# Test data

`1000genomes-phase3-subset.vcf` contains two biallelic SNPs and four samples
from the public 1000 Genomes Project Phase 3 chromosome 20 callset (GRCh37).
The source file is hosted by the [AWS Registry of Open Data][registry]:

```text
https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

It was retrieved on 2026-09-21 with bcftools 1.21:

```bash
bcftools view \
  -r 20:10000117-10000598 \
  -s HG00096,HG00097,HG00099,NA12878 \
  -m2 -M2 -v snps \
  -i 'POS==10000117 || POS==10000598' \
  https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  -Ov
```

Unused metadata declarations were removed from the header after subsetting;
the variant and sample fields are otherwise the direct bcftools output. The
fixture deliberately retains phased genotypes and global population-frequency
INFO fields from the public callset.

[registry]: https://registry.opendata.aws/1000-genomes/
