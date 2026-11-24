from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

gfa_file = Path("data/raw/synthetic_pangenome.gfa")
output_dir = Path("results/figures")
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

colors = {
    "reference": "#64748b",
    "insertion": "#16a34a",
    "deletion": "#ea580c",
    "inversion": "#7c3aed",
}

fig, ax = plt.subplots(figsize=(13, 6))

for row, (name, nodes) in enumerate(paths.items()):
    y = len(paths) - row
    x = 0

    for oriented_node in nodes:
        node = oriented_node[:-1]
        orientation = oriented_node[-1]
        width = len(segments[node])

        if name == "insertion" and node == "E":
            colour = colors["insertion"]
        elif name == "inversion" and orientation == "-":
            colour = colors["inversion"]
        else:
            colour = "#bfdbfe"

        rect = Rectangle(
            (x, y - 0.3),
            width,
            0.6,
            facecolor=colour,
            edgecolor="#334155",
            linewidth=1.2
        )
        ax.add_patch(rect)

        ax.text(
            x + width / 2,
            y,
            f"{node}{orientation}",
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold"
        )

        x += width

    ax.text(
        -2, y, name.capitalize(),
        ha="right", va="center",
        fontsize=11, fontweight="bold"
    )

ax.set_xlim(-15, 65)
ax.set_ylim(0.4, len(paths) + 0.8)
ax.set_yticks([])
ax.set_xlabel("Cumulative path position (bp)")
ax.set_title(
    "Synthetic Pangenome: Structural Variant Paths",
    fontsize=15,
    fontweight="bold"
)

legend = [
    Patch(facecolor="#bfdbfe", label="Forward segment"),
    Patch(facecolor=colors["insertion"], label="Inserted segment"),
    Patch(facecolor=colors["inversion"], label="Reverse-oriented segment"),
]

ax.legend(handles=legend, loc="upper right")
ax.spines[["top", "right", "left"]].set_visible(False)

plt.tight_layout()

output = output_dir / "structural_variant_paths.png"
plt.savefig(output, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {output}")
print(f"Paths visualised: {len(paths)}")
