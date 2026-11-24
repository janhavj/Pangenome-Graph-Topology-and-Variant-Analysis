from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

input_file = Path(
    "results/tables/graph_topology_comparison.csv"
)

output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file)

names = [
    "reference",
    "insertion",
    "deletion",
    "inversion"
]

matrix = pd.DataFrame(
    1.0,
    index=names,
    columns=names
)

for _, row in df.iterrows():
    a = row["path_1"]
    b = row["path_2"]
    similarity = row["adjacency_jaccard"]

    matrix.loc[a, b] = similarity
    matrix.loc[b, a] = similarity

fig, ax = plt.subplots(figsize=(8, 6))

sns.heatmap(
    matrix,
    annot=True,
    fmt=".2f",
    cmap="YlGnBu",
    vmin=0,
    vmax=1,
    square=True,
    linewidths=0.5,
    cbar_kws={
        "label": "Adjacency Jaccard similarity"
    },
    ax=ax
)

ax.set_title(
    "Pangenome Graph Adjacency Similarity",
    fontsize=15,
    fontweight="bold"
)

ax.set_xlabel("Genomic path")
ax.set_ylabel("Genomic path")

plt.tight_layout()

output = output_dir / "graph_topology_heatmap.png"

plt.savefig(
    output,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"Saved: {output}")
print(f"Paths compared: {len(names)}")
print(f"Pairwise comparisons: {len(df)}")
