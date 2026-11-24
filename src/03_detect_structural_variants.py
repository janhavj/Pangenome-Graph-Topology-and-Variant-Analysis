from pathlib import Path
from difflib import SequenceMatcher
import pandas as pd

gfa_file = Path("data/raw/synthetic_pangenome.gfa")
output_file = Path("results/tables/structural_variants.csv")

segments = {}
paths = {}

with gfa_file.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "P":
            paths[fields[1]] = [
                (item[:-1], item[-1])
                for item in fields[2].split(",")
            ]

reference = paths["reference"]
results = []

# Reference coordinates use 0-based half-open intervals.
ref_positions = [0]
for node, orientation in reference:
    ref_positions.append(
        ref_positions[-1] + len(segments[node])
    )

for sample_name, sample_path in paths.items():
    if sample_name == "reference":
        continue

    matcher = SequenceMatcher(
        None,
        reference,
        sample_path,
        autojunk=False
    )

    for operation, i1, i2, j1, j2 in matcher.get_opcodes():
        if operation == "equal":
            continue

        ref_nodes = reference[i1:i2]
        alt_nodes = sample_path[j1:j2]

        ref_length = sum(
            len(segments[node]) for node, _ in ref_nodes
        )
        alt_length = sum(
            len(segments[node]) for node, _ in alt_nodes
        )

        if operation == "insert":
            variant_type = "INS"
        elif operation == "delete":
            variant_type = "DEL"
        elif (
            len(ref_nodes) == len(alt_nodes)
            and len(ref_nodes) > 0
            and alt_nodes == [
                (node, "-" if orientation == "+" else "+")
                for node, orientation in reversed(ref_nodes)
            ]
        ):
            variant_type = "INV"
        else:
            variant_type = "COMPLEX"

        results.append({
            "sample": sample_name,
            "variant_type": variant_type,
            "ref_start": ref_positions[i1],
            "ref_end": ref_positions[i2],
            "ref_length": ref_length,
            "alt_length": alt_length,
            "length_change": alt_length - ref_length,
            "ref_nodes": ",".join(
                n + o for n, o in ref_nodes
            ),
            "alt_nodes": ",".join(
                n + o for n, o in alt_nodes
            )
        })

df = pd.DataFrame(results)
output_file.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_file, index=False)

print("\n=== STRUCTURAL VARIANT CANDIDATES ===")
print(df.to_string(index=False))
print(f"\nTotal candidates: {len(df)}")
print(f"Saved: {output_file}")
