from pathlib import Path
import csv

gfa = Path("data/processed/vg_example/x.gfa")
output = Path("results/tables/reference_skipping_edges.csv")

segments = {}
edges = []
reference = []

with gfa.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "L":
            edges.append((
                fields[1], fields[2],
                fields[3], fields[4]
            ))

        elif fields[0] == "P" and fields[1] == "x":
            reference = fields[2].split(",")

assert all(step.endswith("+") for step in reference), \
    "This script requires a forward-oriented reference path."

reference_nodes = [step[:-1] for step in reference]
positions = {
    node: i for i, node in enumerate(reference_nodes)
}

results = []

for source, source_orient, target, target_orient in edges:
    if source_orient != "+" or target_orient != "+":
        continue

    if source not in positions or target not in positions:
        continue

    start = positions[source]
    end = positions[target]

    if end <= start + 1:
        continue

    skipped = reference_nodes[start + 1:end]

    results.append({
        "start_node": source,
        "end_node": target,
        "skipped_nodes": ",".join(skipped),
        "skipped_reference_bp": sum(
            len(segments[node]) for node in skipped
        ),
        "skipped_reference_steps": len(skipped)
    })

output.parent.mkdir(parents=True, exist_ok=True)

columns = [
    "start_node",
    "end_node",
    "skipped_nodes",
    "skipped_reference_bp",
    "skipped_reference_steps"
]

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(results)

print("\n=== REFERENCE-SKIPPING EDGES ===")
print("Candidate deletion-like edges:", len(results))

for result in results:
    print(result)

print("\nSaved:", output)
