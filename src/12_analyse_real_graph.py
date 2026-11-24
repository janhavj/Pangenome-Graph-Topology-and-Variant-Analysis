from pathlib import Path
from collections import Counter
import pandas as pd
import networkx as nx

gfa_path = Path("data/processed/vg_example/x.gfa")
output_dir = Path("results/tables")
output_dir.mkdir(parents=True, exist_ok=True)

segments = {}
links = []
paths = {}

with gfa_path.open() as handle:
    for line in handle:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments[fields[1]] = fields[2]

        elif fields[0] == "L":
            links.append({
                "from": fields[1],
                "from_orientation": fields[2],
                "to": fields[3],
                "to_orientation": fields[4],
            })

        elif fields[0] == "P":
            paths[fields[1]] = fields[2].split(",")

# Undirected segment-level projection for connectivity analysis.
# This does not preserve GFA endpoint orientation.
G = nx.Graph()
G.add_nodes_from(segments)

for link in links:
    G.add_edge(link["from"], link["to"])

reference_steps = paths.get("x", [])

reference_nodes = {
    step[:-1] for step in reference_steps
}

alternative_nodes = set(segments) - reference_nodes

degrees = dict(G.degree())
degree_counts = Counter(degrees.values())

summary = {
    "total_segments": len(segments),
    "total_links": len(links),
    "reference_path_steps": len(reference_steps),
    "unique_reference_segments": len(reference_nodes),
    "segments_not_on_reference": len(alternative_nodes),
    "connected_components": nx.number_connected_components(G),
    "segment_degree_3_or_more": sum(
        degree >= 3 for degree in degrees.values()
    ),
    "maximum_segment_degree": max(degrees.values(), default=0),
}

print("\n=== REAL VARIATION GRAPH TOPOLOGY ===")

for key, value in summary.items():
    print(f"{key}: {value}")

pd.DataFrame(
    [{"metric": key, "value": value}
     for key, value in summary.items()]
).to_csv(
    output_dir / "real_graph_topology_summary.csv",
    index=False
)

pd.DataFrame(
    [{"degree": degree, "count": count}
     for degree, count in sorted(degree_counts.items())]
).to_csv(
    output_dir / "real_graph_degree_distribution.csv",
    index=False
)

print("\nResults saved to results/tables/")
