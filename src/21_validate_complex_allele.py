from pathlib import Path
import networkx as nx
import pandas as pd
import pysam

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
            reference = fields[2].split(",")

assert all(step.endswith("+") for step in reference)

reference_nodes = [step[:-1] for step in reference]

positions = {}
coordinate = 1

for node in reference_nodes:
    positions[node] = coordinate
    coordinate += len(segments[node])

assert coordinate - 1 == 1001

results = []

with pysam.VariantFile(
    "data/raw/vg_example/x.vcf.gz"
) as vcf:

    for record in vcf:
        for alt in record.alts or ():

            if not (
                len(record.ref) == len(alt)
                and len(alt) > 1
            ):
                continue

            start = record.pos
            end = start + len(record.ref) - 1

            affected = [
                node for node in reference_nodes
                if positions[node] <= end
                and positions[node] + len(segments[node]) - 1 >= start
            ]

            first = reference_nodes.index(affected[0])
            last = reference_nodes.index(affected[-1])

            assert first > 0
            assert last < len(reference_nodes) - 1

            upstream = reference_nodes[first - 1]
            downstream = reference_nodes[last + 1]

            ref_sequence = "".join(
                segments[node] for node in affected
            )

            assert ref_sequence == record.ref

            print("\n=== COMPLEX ALLELE PATH VALIDATION ===")
            print("Position:", start)
            print("Reference allele:", record.ref)
            print("Expected alternative:", alt)
            print("Upstream anchor:", upstream)
            print("Downstream anchor:", downstream)

            matched_paths = []

            for path in nx.all_simple_paths(
                graph,
                source=upstream,
                target=downstream,
                cutoff=len(alt) + 1
            ):
                sequence = "".join(
                    segments[node] for node in path[1:-1]
                )

                print("Path:", " -> ".join(path))
                print("Internal sequence:", sequence)

                if sequence == alt:
                    matched_paths.append(path)

            status = "MATCH" if matched_paths else "NO_MATCH"

            results.append({
                "position": start,
                "reference_allele": record.ref,
                "alternative_allele": alt,
                "matching_paths": len(matched_paths),
                "validated_path": ";".join(
                    "->".join(path) for path in matched_paths
                ),
                "status": status
            })

            print("\nValidation status:", status)

output = Path(
    "results/tables/complex_allele_validation.csv"
)

pd.DataFrame(results).to_csv(output, index=False)

print("\nSaved:", output)
