from pathlib import Path
import pysam

data_dir = Path("data/raw/vg_example")

reference = {}

with pysam.FastxFile(str(data_dir / "x.fa")) as fasta:
    for record in fasta:
        reference[record.name] = record.sequence.upper()

total = 0
matched = 0
mismatches = []

with pysam.VariantFile(str(data_dir / "x.vcf.gz")) as vcf:
    for record in vcf:
        total += 1

        chrom = record.contig
        start = record.pos - 1
        ref = record.ref.upper()

        if chrom not in reference:
            mismatches.append((chrom, record.pos, "Missing contig"))
            continue

        observed = reference[chrom][start:start + len(ref)]

        if observed == ref:
            matched += 1
        else:
            mismatches.append(
                (chrom, record.pos, ref, observed)
            )

print("\n=== REFERENCE VALIDATION ===")
print(f"Total VCF records: {total}")
print(f"Matching REF alleles: {matched}")
print(f"Mismatches: {len(mismatches)}")

if mismatches:
    print("\nFirst mismatches:")
    for item in mismatches[:10]:
        print(item)

if matched == total:
    print("\nPASS: All REF alleles match the FASTA.")
else:
    print("\nFAIL: Reference inconsistencies detected.")
    raise SystemExit(1)
