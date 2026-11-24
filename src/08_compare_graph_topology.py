from pathlib import Path
from itertools import combinations
import pandas as pd

gfa_file = Path("data/raw/synthetic_pangenome.gfa")
output_dir = Path("results/tables")
output_dir.mkdir(parents=True, exist_ok=True)

paths = {}

with gfa_file.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "P":
            paths[fields[1]] = fields[2].split(",")


def reverse_complement_path(path):
    return [
        node[:-1] + ("-" if node[-1] == "+" else "+")
        for node in reversed(path)
    ]


def canonical_adjacency(left, right):
    forward = (left, right)
    reverse = tuple(reverse_complement_path([left, right]))
    return min(forward, reverse)


def get_adjacencies(path):
    return {
        canonical_adjacency(left, right)
        for left, right in zip(path, path[1:])
    }


results = []

for name1, name2 in combinations(paths, 2):
    edges1 = get_adjacencies(paths[name1])
    edges2 = get_adjacencies(paths[name2])

    shared = edges1 & edges2
    union = edges1 | edges2

    similarity = len(shared) / len(union) if union else 1.0

    results.append({
        "path_1": name1,
        "path_2": name2,
        "shared_adjacencies": len(shared),
        "union_adjacencies": len(union),
        "adjacency_jaccard": round(similarity, 4),
        "unique_adjacencies_1": len(edges1 - edges2),
        "unique_adjacencies_2": len(edges2 - edges1),
    })

df = pd.DataFrame(results)

output = output_dir / "graph_topology_comparison.csv"
df.to_csv(output, index=False)

print("\n=== GRAPH ADJACENCY COMPARISON ===")
print(df.to_string(index=False))

print(f"\nSaved: {output}")
