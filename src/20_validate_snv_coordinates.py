from pathlib import Path
import pandas as pd
import pysam

gfa = Path("data/processed/vg_example/x.gfa")

candidates = pd.read_csv(
    "results/tables/bubble_sequence_validation.csv",
    dtype={
        "start_node": str,
        "alternative_node": str,
        "end_node": str
    }
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

positions = {}
coordinate = 1

for step in reference:
    node = step[:-1]
    positions[node] = coordinate
    coordinate += len(segments[node])

assert coordinate - 1 == 1001

graph_snvs = []

for _, row in candidates.iterrows():

    if row["classification"] != "SUBSTITUTION_LIKE":
        continue

    if row["reference_length"] != 1:
        continue

    if row["alternative_length"] != 1:
        continue

    start = row["start_node"]

    graph_snvs.append({
        "position": positions[start] + len(segments[start]),
        "ref": row["reference_sequence"],
        "alt": row["alternative_sequence"]
    })

vcf_snvs = []

with pysam.VariantFile(
    "data/raw/vg_example/x.vcf.gz"
) as vcf:

    for record in vcf:
        for alt in record.alts or ():

            if len(record.ref) == 1 and len(alt) == 1:
                vcf_snvs.append({
                    "position": record.pos,
                    "ref": record.ref,
                    "alt": alt
                })

graph_df = pd.DataFrame(graph_snvs)
vcf_df = pd.DataFrame(vcf_snvs)

comparison = graph_df.merge(
    vcf_df,
    on=["position", "ref", "alt"],
    how="outer",
    indicator=True
).sort_values("position")

output = Path(
    "results/tables/snv_coordinate_validation.csv"
)

comparison.to_csv(output, index=False)

matches = (comparison["_merge"] == "both").sum()

print("\n=== SNV COORDINATE VALIDATION ===")
print("Graph SNV candidates:", len(graph_df))
print("Original VCF SNVs:", len(vcf_df))
print("Exact matches:", matches)

print("\n=== UNMATCHED EVENTS ===")

unmatched = comparison[
    comparison["_merge"] != "both"
]

if unmatched.empty:
    print("None")
else:
    print(unmatched.to_string(index=False))

print("\n=== MATCH SUMMARY ===")
print(comparison["_merge"].value_counts())

print("\nSaved:", output)
