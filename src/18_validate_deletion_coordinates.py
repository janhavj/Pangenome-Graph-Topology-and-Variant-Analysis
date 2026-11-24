from pathlib import Path
import pandas as pd
import pysam

gfa_path = Path("data/processed/vg_example/x.gfa")
skips_path = Path("results/tables/reference_skipping_edges.csv")
output_path = Path("results/tables/deletion_coordinate_validation.csv")

segments = {}
reference_steps = []

with gfa_path.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "P" and fields[1] == "x":
            reference_steps = fields[2].split(",")

assert all(step.endswith("+") for step in reference_steps)

reference_nodes = [step[:-1] for step in reference_steps]

positions = {}
coordinate = 1

for node in reference_nodes:
    positions[node] = coordinate
    coordinate += len(segments[node])

reference_length = coordinate - 1
assert reference_length == 1001, (
    f"Unexpected reference path length: {reference_length}"
)

skips = pd.read_csv(skips_path, dtype={
    "start_node": str,
    "end_node": str,
    "skipped_nodes": str,
})

graph_deletions = []

for _, row in skips.iterrows():
    skipped_nodes = row["skipped_nodes"].split(",")

    start = positions[skipped_nodes[0]]

    deleted_sequence = "".join(
        segments[node] for node in skipped_nodes
    )

    graph_deletions.append({
        "deleted_start_1based": start,
        "deleted_sequence": deleted_sequence,
        "deleted_length": len(deleted_sequence),
    })

vcf_deletions = []

with pysam.VariantFile(
    "data/raw/vg_example/x.vcf.gz"
) as vcf:
    for record in vcf:
        for alt in record.alts or ():
            if len(record.ref) > len(alt):
                # For the simple left-anchored deletions
                # in this example dataset.
                assert record.ref.startswith(alt)
                vcf_deletions.append({
                    "deleted_start_1based": (
                        record.pos + len(alt)
                    ),
                    "deleted_sequence": record.ref[len(alt):],
                    "deleted_length": len(record.ref) - len(alt),
                })

graph_df = pd.DataFrame(graph_deletions)
vcf_df = pd.DataFrame(vcf_deletions)

merged = graph_df.merge(
    vcf_df,
    on=[
        "deleted_start_1based",
        "deleted_sequence",
        "deleted_length"
    ],
    how="outer",
    indicator=True
)

merged.to_csv(output_path, index=False)

print("\n=== DELETION COORDINATE VALIDATION ===")
print(merged.to_string(index=False))

matches = (merged["_merge"] == "both").sum()

print(f"\nExact matches: {matches}")
print(f"Graph candidates: {len(graph_df)}")
print(f"VCF deletions: {len(vcf_df)}")

if (
    matches == len(graph_df)
    and matches == len(vcf_df)
):
    print("\nPASS: All deletion coordinates and sequences match.")
else:
    print("\nFAIL: Some deletion events do not match.")
    raise SystemExit(1)

print("\nSaved:", output_path)
