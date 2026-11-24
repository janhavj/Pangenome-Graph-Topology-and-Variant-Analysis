from pathlib import Path

output = Path("data/raw/synthetic_pangenome.gfa")
output.parent.mkdir(parents=True, exist_ok=True)

segments = {
    "A": "ATGCGTACGA",
    "B": "TTCGATCGAA",
    "C": "GGCTAACGTT",
    "D": "CCGTAAGCTA",
    "E": "AATTCCGGTA",
    "F": "GCTAGCTAGC",
}

paths = {
    "reference": ["A+", "B+", "C+", "D+", "F+"],
    "insertion": ["A+", "B+", "C+", "E+", "D+", "F+"],
    "deletion": ["A+", "B+", "D+", "F+"],
    "inversion": ["A+", "B+", "C-", "D+", "F+"],
}

# GFA links encode oriented segment connections.
links = set()

for nodes in paths.values():
    for left, right in zip(nodes, nodes[1:]):
        links.add((
            left[:-1], left[-1],
            right[:-1], right[-1]
        ))

with output.open("w") as handle:
    handle.write("H\tVN:Z:1.0\n")

    for name, sequence in segments.items():
        handle.write(f"S\t{name}\t{sequence}\n")

    for source, so, target, to in sorted(links):
        handle.write(
            f"L\t{source}\t{so}\t{target}\t{to}\t0M\n"
        )

    for name, nodes in paths.items():
        overlaps = ",".join("*" for _ in range(len(nodes) - 1))
        handle.write(
            f"P\t{name}\t{','.join(nodes)}\t{overlaps}\n"
        )

print(f"Generated: {output}")
print(f"Segments: {len(segments)}")
print(f"Links: {len(links)}")
print(f"Paths: {len(paths)}")
