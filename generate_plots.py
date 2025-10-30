#!/usr/bin/env python3
"""
Generate all plots needed for the APSP performance report
Requires: matplotlib, pandas
Install: pip3 install matplotlib pandas numpy --user
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# (Optional) simple, readable defaults
plt.rcParams["figure.figsize"] = (10, 6)
plt.rcParams["font.size"] = 12

# Read the data
df = pd.read_csv("report_data/speedup_data.csv")


def sorted_ps(series):
    return sorted(series.unique())


# ---------- Figure 1: Speedup vs Sequential Blocked ----------
plt.figure(figsize=(10, 6))
p_values = sorted_ps(df["p"])
for graph in ["A", "B", "C", "D"]:
    data = df[df["graph"] == graph].sort_values("p")
    plt.plot(
        data["p"].values.astype(int),
        data["speedup_vs_seq"].astype(float).values,
        marker="o",
        linewidth=2,
        markersize=8,
        label=f"Graph {graph}",
    )

plt.plot(p_values, p_values, "k--", linewidth=1, alpha=0.5, label="Ideal")
plt.xlabel("Number of Processes (P)")
plt.ylabel("Speedup")
plt.title("Speedup Relative to Sequential Blocked Floyd–Warshall")
plt.grid(True, alpha=0.3)
plt.legend()
plt.xticks(p_values)
plt.tight_layout()
plt.savefig("report_data/speedup_vs_sequential.png", dpi=300, bbox_inches="tight")
print("✓ Generated: speedup_vs_sequential.png")

# ---------- Figure 2: Speedup vs P=1 (Parallel Implementation) ----------
plt.figure(figsize=(10, 6))
p_values = sorted_ps(df["p"])
for graph in ["A", "B", "C", "D"]:
    data = df[df["graph"] == graph].sort_values("p")
    plt.plot(
        data["p"].values.astype(int),
        data["speedup_vs_p1"].astype(float).values,
        marker="s",
        linewidth=2,
        markersize=8,
        label=f"Graph {graph}",
    )

plt.plot(p_values, p_values, "k--", linewidth=1, alpha=0.5, label="Ideal")
plt.xlabel("Number of Processes (P)")
plt.ylabel("Speedup")
plt.title("Speedup Relative to Parallel Implementation (P=1)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.xticks(p_values)
plt.tight_layout()
plt.savefig("report_data/speedup_vs_p1.png", dpi=300, bbox_inches="tight")
print("✓ Generated: speedup_vs_p1.png")

# ---------- Figure 3: Strong Scaling (log-log) ----------
plt.figure(figsize=(10, 6))
p_values = sorted_ps(df["p"])
for graph in ["A", "B", "C", "D"]:
    data = df[df["graph"] == graph].sort_values("p")
    # normalized time ∝ 1 / speedup_vs_p1
    relative_times = 1.0 / data["speedup_vs_p1"].astype(float).values
    plt.loglog(
        data["p"].values.astype(int),
        relative_times,
        marker="o",
        linewidth=2,
        markersize=8,
        label=f"Graph {graph}",
    )

p_range = np.array(p_values, dtype=float)
ideal = 1.0 / p_range
plt.loglog(p_range, ideal, "k--", linewidth=1, alpha=0.5, label="Ideal")
plt.xlabel("Number of Processes (P)")
plt.ylabel("Normalized Execution Time")
plt.title("Strong Scaling Analysis (Log-Log Plot)")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.xticks(p_values, p_values)
plt.tight_layout()
plt.savefig("report_data/strong_scaling.png", dpi=300, bbox_inches="tight")
print("✓ Generated: strong_scaling.png")

# ---------- Figure 4: Parallel Efficiency ----------
plt.figure(figsize=(10, 6))
df_eff = df[df["p"] > 1].copy()
df_eff["efficiency"] = df_eff["efficiency"].astype(float)
p_values_eff = sorted_ps(df_eff["p"])
for graph in ["A", "B", "C", "D"]:
    data = df_eff[df_eff["graph"] == graph].sort_values("p")
    plt.plot(
        data["p"].values.astype(int),
        data["efficiency"].values,
        marker="D",
        linewidth=2,
        markersize=8,
        label=f"Graph {graph}",
    )

plt.axhline(
    y=100, color="k", linestyle="--", linewidth=1, alpha=0.5, label="Perfect Efficiency"
)
plt.xlabel("Number of Processes (P)")
plt.ylabel("Parallel Efficiency (%)")
plt.title("Parallel Efficiency Analysis")
plt.grid(True, alpha=0.3)
plt.legend()
plt.xticks(p_values_eff)
plt.ylim(0, 105)  # full range
plt.tight_layout()
plt.savefig("report_data/parallel_efficiency.png", dpi=300, bbox_inches="tight")
print("✓ Generated: parallel_efficiency.png")

# ---------- Figure 5: Speedup Comparison (All Graphs) ----------
plt.figure(figsize=(12, 6))
bars = [4, 9, 16, 25]
x = np.arange(len(bars))
width = 0.2
for i, graph in enumerate(["A", "B", "C", "D"]):
    data = df[(df["graph"] == graph) & (df["p"].isin(bars))].sort_values("p")
    speedups = data["speedup_vs_seq"].astype(float).values
    plt.bar(x + i * width, speedups, width, label=f"Graph {graph}")

plt.xlabel("Number of Processes (P)")
plt.ylabel("Speedup (vs Sequential)")
plt.title("Speedup Comparison Across All Graphs")
plt.xticks(x + width * 1.5, list(map(str, bars)))
plt.legend()
plt.grid(True, axis="y", alpha="0.3")
plt.tight_layout()
plt.savefig("report_data/speedup_comparison.png", dpi=300, bbox_inches="tight")
print("✓ Generated: speedup_comparison.png")

print("\nAll plots generated successfully in report_data/")
print("\nGenerated files:")
print("  1. speedup_vs_sequential.png")
print("  2. speedup_vs_p1.png")
print("  3. strong_scaling.png")
print("  4. parallel_efficiency.png")
print("  5. speedup_comparison.png")
