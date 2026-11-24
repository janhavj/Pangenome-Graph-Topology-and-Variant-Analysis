from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

gfa_file = Path("data/raw/synthetic_pangenome.gfa")
output_dir = Path("results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

graph = nx.DiGraph()
paths = {}

with gfa_file.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            graph.add_node(fields[1], length=len(fields[2]))

        elif fields[0] == "L":
            graph.add_edge(fields[1], fields[3])

        elif fields[0] == "P":
            paths[fields[1]] = [
                item[:-1]
                for item in fields[2].split(",")
            ]

positions = {
    "A": (0, 0),
    "B": (1, 0),
    "C": (2, 0),
    "D": (3, 0),
    "E": (2.5, 1),
    "F": (4, 0),
}

plt.figure(figsize=(11, 5))

nx.draw_networkx_nodes(
    graph, positions,
    node_color="#dbeafe",
    node_size=1700,
    edgecolors="#1d4ed8"
)

nx.draw_networkx_labels(
    graph, positions,
    font_size=12,
    font_weight="bold"
)

nx.draw_networkx_edges(
    graph, positions,
    edge_color="#94a3b8",
    arrows=True,
    arrowsize=20,
    node_size=1700,
    connectionstyle="arc3,rad=0.08"
)

colors = {
    "reference": "#2563eb",
    "insertion": "#16a34a",
    "deletion": "#ea580c",
}

for name, nodes in paths.items():
    path_edges = list(zip(nodes, nodes[1:]))

    nx.draw_networkx_edges(
        graph, positions,
        edgelist=path_edges,
        edge_color=colors[name],
        width=2.5,
        alpha=0.8,
        arrows=True,
        arrowsize=18,
        node_size=1700,
        connectionstyle="arc3,rad=0.08"
    )

plt.title("Synthetic Pangenome Graph", fontsize=16)
plt.axis("off")
plt.tight_layout()

output = output_dir / "synthetic_pangenome_graph.png"
plt.savefig(output, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {output}")
print(f"Nodes: {graph.number_of_nodes()}")
print(f"Edges: {graph.number_of_edges()}")
print(f"Paths visualised: {len(paths)}")
