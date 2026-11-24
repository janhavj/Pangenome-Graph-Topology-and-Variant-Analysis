from pathlib import Path
from itertools import combinations
import pandas as pd

gfa_file = Path("data/raw/synthetic_pangenome.gfa")
output_dir = Path("results/tables")
output_dir.mkdir(parents=True, exist_ok=True)

segments = {}
paths = {}

with gfa_file.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "P":
            paths[fields[1]] = fields[2].split(",")

results = []

for name1, name2 in combinations(paths, 2):
    path1 = set(paths[name1])
    path2 = set(paths[name2])

    shared = path1 & path2
    union = path1 | path2

    jaccard = len(shared) / len(union) if union else 1.0

    results.append({
        "path_1": name1,
        "path_2": name2,
        "shared_oriented_segments": len(shared),
        "union_oriented_segments": len(union),
        "jaccard_similarity": round(jaccard, 4),
        "unique_to_path_1": ",".join(sorted(path1 - path2)),
        "unique_to_path_2": ",".join(sorted(path2 - path1)),
    })

df = pd.DataFrame(results)
df.to_csv(
    output_dir / "pairwise_path_comparison.csv",
    index=False
)

print("\n=== PAIRWISE PATH COMPARISON ===")
print(df.to_string(index=False))

common = set.intersection(
    *(set(path) for path in paths.values())
)

print("\n=== CONSERVED ORIENTED SEGMENTS ===")
print(sorted(common))

print("\n=== PATH LENGTHS ===")
for name, path in paths.items():
    length = sum(
        len(segments[item[:-1]])
        for item in path
    )
    print(f"{name}: {length} bp")

print("\nSaved: results/tables/pairwise_path_comparison.csv")
