from pathlib import Path
import pandas as pd

gfa_path = Path("data/processed/vg_example/x.gfa")
candidate_path = Path(
    "results/tables/simple_bubble_candidates.csv"
)

segments = {}
reference_steps = []

with gfa_path.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "P" and fields[1] == "x":
            reference_steps = [
                step[:-1]
                for step in fields[2].split(",")
            ]

positions = {
    node: i for i, node in enumerate(reference_steps)
}

df = pd.read_csv(candidate_path, dtype={
    "start_reference_node": str,
    "alternative_node": str,
    "end_reference_node": str
})

results = []

for _, row in df.iterrows():
    start = row["start_reference_node"]
    alt = row["alternative_node"]
    end = row["end_reference_node"]

    start_idx = positions[start]
    end_idx = positions[end]

    # Reference sequence between the two anchors.
    ref_internal_nodes = reference_steps[
        start_idx + 1:end_idx
    ]

    ref_sequence = "".join(
        segments[node] for node in ref_internal_nodes
    )

    alt_sequence = segments[alt]

    length_change = len(alt_sequence) - len(ref_sequence)

    if length_change > 0:
        classification = "INSERTION_LIKE"
    elif length_change < 0:
        classification = "DELETION_LIKE"
    elif alt_sequence != ref_sequence:
        classification = "SUBSTITUTION_LIKE"
    else:
        classification = "IDENTICAL"

    results.append({
        "start_node": start,
        "alternative_node": alt,
        "end_node": end,
        "reference_sequence": ref_sequence,
        "alternative_sequence": alt_sequence,
        "reference_length": len(ref_sequence),
        "alternative_length": len(alt_sequence),
        "length_change": length_change,
        "classification": classification
    })

result_df = pd.DataFrame(results)

output = Path(
    "results/tables/bubble_sequence_validation.csv"
)

result_df.to_csv(output, index=False)

print("\n=== CANDIDATE BRANCH CLASSIFICATION ===")
print(result_df["classification"].value_counts())

print("\n=== ONE-STEP BRANCHES ===")
print(
    result_df.loc[
        result_df["reference_length"] == 0
    ].to_string(index=False)
)

print("\nSaved:", output)
