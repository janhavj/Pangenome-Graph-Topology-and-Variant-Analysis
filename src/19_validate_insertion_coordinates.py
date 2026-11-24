from pathlib import Path
import pandas as pd
import pysam

gfa = Path("data/processed/vg_example/x.gfa")
candidates = pd.read_csv(
    "results/tables/bubble_sequence_validation.csv",
    dtype={"start_node": str, "alternative_node": str, "end_node": str}
)

segments = {}
reference = []

with gfa.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if fields[0] == "S":
            segments[fields[1]] = fields[2]
        elif fields[0] == "P" and fields[1] == "x":
            reference = fields[2].split(",")

assert all(step.endswith("+") for step in reference)
reference_nodes = [step[:-1] for step in reference]

positions = {}
coordinate = 1
for node in reference_nodes:
    positions[node] = coordinate
    coordinate += len(segments[node])

assert coordinate - 1 == 1001

graph_insertions = []
for _, row in candidates.iterrows():
    if row["classification"] != "INSERTION_LIKE":
        continue

    start = row["start_node"]
    inserted = row["alternative_sequence"]

    graph_insertions.append({
        "anchor_position": positions[start] + len(segments[start]) - 1,
        "inserted_sequence": inserted,
        "inserted_length": len(inserted)
    })

vcf_insertions = []
with pysam.VariantFile("data/raw/vg_example/x.vcf.gz") as vcf:
    for record in vcf:
        for alt in record.alts or ():
            if len(alt) > len(record.ref):
                assert alt.startswith(record.ref)
                inserted = alt[len(record.ref):]
                vcf_insertions.append({
                    "anchor_position": record.pos + len(record.ref) - 1,
                    "inserted_sequence": inserted,
                    "inserted_length": len(inserted)
                })

graph_df = pd.DataFrame(graph_insertions)
vcf_df = pd.DataFrame(vcf_insertions)

comparison = graph_df.merge(
    vcf_df,
    on=["anchor_position", "inserted_sequence", "inserted_length"],
    how="outer",
    indicator=True
)

output = Path("results/tables/insertion_coordinate_validation.csv")
comparison.to_csv(output, index=False)

print("\n=== INSERTION COORDINATE VALIDATION ===")
print(comparison.to_string(index=False))
print("\nExact matches:", (comparison["_merge"] == "both").sum())
print("Graph candidates:", len(graph_df))
print("VCF insertions:", len(vcf_df))

if (
    len(graph_df) == len(vcf_df)
    and (comparison["_merge"] == "both").all()
):
    print("\nPASS: All insertion coordinates and sequences match.")
else:
    print("\nREVIEW: Some insertion events do not match.")

print("\nSaved:", output)
