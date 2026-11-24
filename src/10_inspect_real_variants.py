from pathlib import Path
import pysam

data_dir = Path("data/raw/vg_example")

fasta_file = data_dir / "x.fa"
vcf_file = data_dir / "x.vcf.gz"

print("\n=== REFERENCE GENOME ===")

references = {}

with pysam.FastxFile(str(fasta_file)) as fasta:
    for record in fasta:
        references[record.name] = record.sequence

        print(f"Sequence: {record.name}")
        print(f"Length: {len(record.sequence)} bp")

print("\n=== VCF VARIANTS ===")

variant_counts = {}
total = 0

with pysam.VariantFile(str(vcf_file)) as vcf:
    print("Samples:", list(vcf.header.samples))

    for record in vcf:
        total += 1

        ref = record.ref
        alts = record.alts or ()

        print(
            f"{record.contig}:{record.pos} "
            f"REF={ref} ALT={','.join(alts)}"
        )

        for alt in alts:
            if alt.startswith("<") or "[" in alt or "]" in alt:
                variant_type = "SYMBOLIC_OR_BREAKEND"
            elif len(ref) == len(alt) == 1:
                variant_type = "SNV"
            elif len(ref) < len(alt):
                variant_type = "INSERTION"
            elif len(ref) > len(alt):
                variant_type = "DELETION"
            else:
                variant_type = "COMPLEX"

            variant_counts[variant_type] = (
                variant_counts.get(variant_type, 0) + 1
            )

print("\n=== VARIANT SUMMARY ===")
print(f"Total VCF records: {total}")

for variant_type, count in sorted(variant_counts.items()):
    print(f"{variant_type}: {count}")

print("\nInspection completed.")
