from pathlib import Path
from collections import Counter

gfa = Path("data/processed/vg_example/x.gfa")

segments = set()
oriented_links = []
reference_steps = []

with gfa.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        if fields[0] == "S":
            segments.add(fields[1])

        elif fields[0] == "L":
            oriented_links.append((
                fields[1], fields[2],
                fields[3], fields[4]
            ))

        elif fields[0] == "P" and fields[1] == "x":
            reference_steps = fields[2].split(",")

reference_nodes = {step[:-1] for step in reference_steps}

nonreference_nodes = segments - reference_nodes

link_classes = Counter()

for source, source_orient, target, target_orient in oriented_links:
    source_ref = source in reference_nodes
    target_ref = target in reference_nodes

    if source_ref and target_ref:
        category = "reference_to_reference"
    elif not source_ref and not target_ref:
        category = "alternative_to_alternative"
    else:
        category = "reference_to_alternative"

    link_classes[category] += 1

print("\n=== ORIENTATION-AWARE LINK INSPECTION ===")
print("Total segments:", len(segments))
print("Total oriented links:", len(oriented_links))
print("Reference segments:", len(reference_nodes))
print("Non-reference segments:", len(nonreference_nodes))

print("\nLink categories:")
for category, count in sorted(link_classes.items()):
    print(f"{category}: {count}")

print("\nOrientation combinations:")
for orientations, count in sorted(
    Counter(
        source_orient + target_orient
        for _, source_orient, _, target_orient in oriented_links
    ).items()
):
    print(f"{orientations}: {count}")
