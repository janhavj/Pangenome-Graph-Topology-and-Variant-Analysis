# Data Sources

## 1. Synthetic Pangenome Benchmark

The synthetic GFA graph was generated specifically for this
project using `src/01_generate_test_graph.py`.

It contains reference, insertion, deletion and inversion
paths for controlled algorithm validation.

## 2. VG Benchmark Dataset

Source: VG Toolkit official GitHub repository

Repository:
https://github.com/vgteam/vg

Reference FASTA:
https://github.com/vgteam/vg/blob/master/test/small/x.fa

Variant VCF:
https://github.com/vgteam/vg/blob/master/test/small/x.vcf.gz

Dataset characteristics:

- Reference sequence length: 1,001 bp
- Original VCF records: 75
- SNVs: 61
- Insertions: 8
- Deletions: 5
- Multi-nucleotide variants: 1

The dataset is a small software-testing benchmark,
not a population-scale biological pangenome.

## 3. Software Environment

- Python 3.11.0
- VG v1.78.0 (Sharptooth)
- Docker platform: linux/amd64

Exact Docker image:

quay.io/vgteam/vg@sha256:febc2401166d4acd81cbf07d89158feccee8c13de70613f72e641d48ceea4ca8

## 4. File Integrity

SHA-256 checksums are provided in `data_checksums.sha256`.

## 5. Attribution and Redistribution

The upstream benchmark files originate from the VG Toolkit
repository. Dataset-specific redistribution permissions
must be verified before distributing copies.

## 6. Upstream Software Licence

The VG Toolkit is distributed under the MIT Licence.

Copyright (c) 2014 Erik Garrison.

The upstream software licence and associated third-party
licence notices are preserved in VG_UPSTREAM_LICENSE.txt.

The original benchmark FASTA and VCF files and generated
VG/GFA graphs are excluded from this repository.

Dataset-specific licensing has not been independently
confirmed. The software licence does not establish
separate redistribution rights for every benchmark file.
