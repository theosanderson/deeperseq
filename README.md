# deeperseq

### https://deeperseq.genomium.org


This uses [fastq-dl](https://github.com/rpetit3/fastq-dl), minimap2, and more!

### Linking directly to a dataset

Analyses can be started from query parameters:

- `?acc=SRR10903401` – an SRA/ENA run accession
- `?url=https://example.com/reads.fastq.gz` – a FASTQ file (gzipped or not); add `&url2=...` for the second file of paired-end reads
- `&ref=NC_001802` – the reference: one of the ids in the dropdown (e.g. `ASM985889v3` for SARS-CoV-2), or any other GenBank ID
