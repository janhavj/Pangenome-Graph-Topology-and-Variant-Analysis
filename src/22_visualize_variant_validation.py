from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

tables = Path("results/tables")
figures = Path("results/figures")
figures.mkdir(parents=True, exist_ok=True)

results = []

files = {
    "SNVs": "snv_coordinate_validation.csv",
    "Insertions": "insertion_coordinate_validation.csv",
    "Deletions": "deletion_coordinate_validation.csv"
}

for category, filename in files.items():
    df = pd.read_csv(tables / filename)

    expected = df["_merge"].isin(
        ["both", "right_only"]
    ).sum()

    matched = (df["_merge"] == "both").sum()

    results.append({
        "variant_type": category,
        "vcf_records": expected,
        "validated_records": matched
    })

complex_df = pd.read_csv(
    tables / "complex_allele_validation.csv"
)

results.append({
    "variant_type": "MNVs",
    "vcf_records": len(complex_df),
    "validated_records": (
        complex_df["status"] == "MATCH"
    ).sum()
})

summary = pd.DataFrame(results)

summary["concordance_percent"] = (
    100 * summary["validated_records"]
    / summary["vcf_records"]
)

summary.to_csv(
    tables / "final_variant_validation_summary.csv",
    index=False
)

x = np.arange(len(summary))
width = 0.34

fig, ax = plt.subplots(figsize=(10, 6))

bars1 = ax.bar(
    x - width / 2,
    summary["vcf_records"],
    width,
    label="Original VCF records",
    color="#547AA5"
)

bars2 = ax.bar(
    x + width / 2,
    summary["validated_records"],
    width,
    label="Validated graph representation",
    color="#319177"
)

ax.bar_label(bars1, padding=3)
ax.bar_label(bars2, padding=3)

ax.set_xticks(x)
ax.set_xticklabels(summary["variant_type"])
ax.set_ylabel("Number of variant records")
ax.set_xlabel("Variant category")

ax.set_title(
    "Variant Representation Concordance\n"
    "VCF-derived Variation Graph",
    fontweight="bold",
    pad=16
)

ax.set_ylim(
    0, summary["vcf_records"].max() * 1.25
)

ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", alpha=0.2)
ax.set_axisbelow(True)

total = summary["vcf_records"].sum()
validated = summary["validated_records"].sum()

fig.text(
    0.5, 0.01,
    f"Benchmark: {validated}/{total} records validated "
    f"({100 * validated / total:.1f}% concordance)",
    ha="center",
    fontsize=10
)

plt.tight_layout(rect=[0, 0.04, 1, 1])

output = figures / "variant_representation_concordance.png"

fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print("\n=== FINAL VARIANT VALIDATION SUMMARY ===")
print(summary.to_string(index=False))

print("\nTotal VCF records:", total)
print("Validated records:", validated)
print(f"Concordance: {100 * validated / total:.2f}%")

print("\nFigure saved:", output)
