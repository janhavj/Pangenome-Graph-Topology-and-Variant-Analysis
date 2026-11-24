from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

input_file = Path(
    "results/tables/pairwise_path_comparison.csv"
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
    value = row["jaccard_similarity"]

    matrix.loc[a, b] = value
    matrix.loc[b, a] = value

plt.figure(figsize=(8, 6))

sns.heatmap(
    matrix,
    annot=True,
    fmt=".2f",
    cmap="YlGnBu",
    vmin=0,
    vmax=1,
    square=True,
    linewidths=0.5,
    cbar_kws={"label": "Jaccard similarity"}
)

plt.title(
    "Pangenome Path Similarity",
    fontsize=15,
    fontweight="bold"
)

plt.xlabel("Genomic path")
plt.ylabel("Genomic path")
plt.tight_layout()

output = output_dir / "path_similarity_heatmap.png"
plt.savefig(output, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {output}")
print(f"Paths compared: {len(names)}")
print(f"Pairwise comparisons: {len(df)}")
