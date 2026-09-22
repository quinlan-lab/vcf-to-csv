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

Unused metadata declarations were removed from the header after subsetting.
The fixture deliberately retains phased genotypes and global
population-frequency INFO fields from the public callset.

The subset was annotated on 2026-09-22 with fastVEP 0.3.0 at commit
`0d31fe0ce0d9be9901944dbd9c3f646363766d75`, using Ensembl's frozen GRCh37.87
chromosome-20 gene models and reference sequence:

```bash
fastvep annotate \
  --input 1000genomes-phase3-subset.unannotated.vcf \
  --output 1000genomes-phase3-subset.vcf \
  --gff3 Homo_sapiens.GRCh37.87.chromosome.20.gff3 \
  --fasta Homo_sapiens.GRCh37.dna.chromosome.20.fa \
  --symbol \
  --canonical \
  --hgvs \
  --no-progress
```

The annotation resources came from the Ensembl GRCh37 release-115 archive:

```text
https://ftp.ensembl.org/pub/grch37/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chromosome.20.gff3.gz
https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz
```

[registry]: https://registry.opendata.aws/1000-genomes/
