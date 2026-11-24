from pathlib import Path
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

gfa = Path("data/processed/vg_example/x.gfa")

segments = {}
graph = nx.DiGraph()
reference = []

with gfa.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]
            graph.add_node(fields[1])

        elif fields[0] == "L":
            if fields[2] == "+" and fields[4] == "+":
                graph.add_edge(fields[1], fields[3])

        elif fields[0] == "P" and fields[1] == "x":
            reference = [
                step[:-1]
                for step in fields[2].split(",")
            ]

reference_set = set(reference)

coordinates = {}
position = 1

for node in reference:
    coordinates[node] = (
        position,
        position + len(segments[node]) - 1
    )
    position += len(segments[node])

affected = [
    node for node in reference
    if coordinates[node][0] <= 273
    and coordinates[node][1] >= 272
]

first = reference.index(affected[0])
last = reference.index(affected[-1])

upstream = reference[first - 1]
downstream = reference[last + 1]

paths = list(nx.all_simple_paths(
    graph,
    source=upstream,
    target=downstream,
    cutoff=6
))

local_nodes = {
    node for path in paths for node in path
}

local_graph = graph.subgraph(local_nodes).copy()

assert nx.is_directed_acyclic_graph(local_graph)

for layer, nodes in enumerate(
    nx.topological_generations(local_graph)
):
    for node in nodes:
        local_graph.nodes[node]["layer"] = layer

layout = nx.multipartite_layout(
    local_graph,
    subset_key="layer",
    align="vertical"
)

colours = [
    "#547AA5" if node in reference_set
    else "#319177"
    for node in local_graph.nodes
]

labels = {}

for node in local_graph.nodes:
    sequence = segments[node]

    if len(sequence) > 10:
        sequence = sequence[:7] + "..."

    labels[node] = f"{node}\n{sequence}"

fig, ax = plt.subplots(figsize=(12, 6))

nx.draw_networkx_nodes(
    local_graph,
    layout,
    node_color=colours,
    node_size=2500,
    edgecolors="white",
    linewidths=2,
    ax=ax
)

nx.draw_networkx_edges(
    local_graph,
    layout,
    arrows=True,
    arrowsize=22,
    width=1.8,
    edge_color="#64748B",
    node_size=2500,
    ax=ax
)

nx.draw_networkx_labels(
    local_graph,
    layout,
    labels=labels,
    font_color="white",
    font_size=10,
    font_weight="bold",
    ax=ax
)

ax.legend(
    handles=[
        Patch(color="#547AA5", label="Reference segment"),
        Patch(color="#319177", label="Alternative segment")
    ],
    loc="upper right",
    frameon=False
)

ax.set_title(
    "Local Variation Graph Topology\n"
    "Multi-nucleotide Variant at Positions 272–273",
    fontsize=15,
    fontweight="bold",
    pad=20
)

fig.text(
    0.5, 0.04,
    "Reference allele: TA | Alternative allele: CG",
    ha="center",
    fontsize=11
)

ax.set_axis_off()
plt.tight_layout(rect=[0, 0.08, 1, 1])

output = Path(
    "results/figures/local_variation_graph.png"
)

fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print("\n=== LOCAL VARIATION GRAPH ===")
print("Reference anchors:", upstream, downstream)
print("Local segments:", local_graph.number_of_nodes())
print("Local edges:", local_graph.number_of_edges())
print("Traversable paths:", len(paths))

for path in paths:
    sequence = "".join(
        segments[node] for node in path[1:-1]
    )
    print(" -> ".join(path), "| Allele:", sequence)

print("\nFigure saved:", output)
