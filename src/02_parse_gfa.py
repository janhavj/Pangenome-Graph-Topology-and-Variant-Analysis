from pathlib import Path
import networkx as nx
import pandas as pd

input_file = Path("data/raw/synthetic_pangenome.gfa")

segments = {}
links = []
paths = {}

with input_file.open() as file:
    for line in file:
        fields = line.strip().split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "L":
            links.append({
                "source": fields[1],
                "source_orientation": fields[2],
                "target": fields[3],
                "target_orientation": fields[4],
            })

        elif fields[0] == "P":
            paths[fields[1]] = [
                (node[:-1], node[-1])
                for node in fields[2].split(",")
            ]

graph = nx.MultiDiGraph()

for node, sequence in segments.items():
    graph.add_node(
        node,
        sequence=sequence,
        length=len(sequence)
    )

for link in links:
    graph.add_edge(
        link["source"],
        link["target"],
        source_orientation=link["source_orientation"],
        target_orientation=link["target_orientation"]
    )

print("\n=== GFA GRAPH SUMMARY ===")
print(f"Segments: {graph.number_of_nodes()}")
print(f"Links: {graph.number_of_edges()}")
print(f"Paths: {len(paths)}")

print("\n=== PATHS ===")
for name, nodes in paths.items():
    formatted = " -> ".join(
        f"{node}{orientation}"
        for node, orientation in nodes
    )
    print(f"{name}: {formatted}")

print("\n=== GRAPH STATISTICS ===")
print(f"Total segment length: {sum(len(s) for s in segments.values())} bp")
print(f"Weakly connected components: {nx.number_weakly_connected_components(graph)}")

output_dir = Path("results/tables")
output_dir.mkdir(parents=True, exist_ok=True)

pd.DataFrame([
    {
        "segment": node,
        "sequence": sequence,
        "length": len(sequence),
        "in_degree": graph.in_degree(node),
        "out_degree": graph.out_degree(node)
    }
    for node, sequence in segments.items()
]).to_csv(output_dir / "graph_segments.csv", index=False)

print("\nSaved: results/tables/graph_segments.csv")
