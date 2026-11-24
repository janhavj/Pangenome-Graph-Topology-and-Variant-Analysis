from pathlib import Path
from collections import defaultdict
import csv

gfa_path = Path("data/processed/vg_example/x.gfa")
output_path = Path("results/tables/simple_bubble_candidates.csv")

segments = {}
edges = defaultdict(set)
reference_steps = []

with gfa_path.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "L":
            source, source_orient = fields[1], fields[2]
            target, target_orient = fields[3], fields[4]

            if source_orient == "+" and target_orient == "+":
                edges[source].add(target)

        elif fields[0] == "P" and fields[1] == "x":
            reference_steps = fields[2].split(",")

reference_nodes = [step[:-1] for step in reference_steps]
reference_set = set(reference_nodes)
reference_positions = {
    node: i for i, node in enumerate(reference_nodes)
}

predecessors = defaultdict(set)
for source, targets in edges.items():
    for target in targets:
        predecessors[target].add(source)

candidates = []

for alt_node in segments:
    if alt_node in reference_set:
        continue

    upstream = predecessors[alt_node] & reference_set
    downstream = edges[alt_node] & reference_set

    for start in upstream:
        for end in downstream:
            if reference_positions[start] >= reference_positions[end]:
                continue

            candidates.append({
                "start_reference_node": start,
                "alternative_node": alt_node,
                "end_reference_node": end,
                "alternative_length_bp": len(segments[alt_node]),
                "reference_step_span": (
                    reference_positions[end]
                    - reference_positions[start]
                ),
            })

output_path.parent.mkdir(parents=True, exist_ok=True)

columns = [
    "start_reference_node",
    "alternative_node",
    "end_reference_node",
    "alternative_length_bp",
    "reference_step_span",
]

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(candidates)

print("\n=== SIMPLE BUBBLE CANDIDATES ===")
print("Candidates detected:", len(candidates))

for candidate in candidates[:10]:
    print(candidate)

print("\nSaved:", output_path)
