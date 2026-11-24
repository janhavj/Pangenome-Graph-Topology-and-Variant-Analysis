from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

table = Path("results/tables/real_graph_degree_distribution.csv")
output = Path("results/figures/real_graph_degree_distribution.png")

output.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(table).sort_values("degree")

fig, ax = plt.subplots(figsize=(8, 5))

ax.bar(
    df["degree"].astype(str),
    df["count"],
    color="#4472C4",
    edgecolor="black",
    linewidth=0.7
)

for i, count in enumerate(df["count"]):
    ax.text(i, count + 2, str(count), ha="center", fontsize=10)

ax.set_xlabel("Segment degree", fontsize=11)
ax.set_ylabel("Number of segments", fontsize=11)
ax.set_title(
    "Degree Distribution of a VCF-derived Variation Graph",
    fontsize=13,
    fontweight="bold"
)

ax.set_ylim(0, df["count"].max() * 1.18)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Figure saved: {output}")
