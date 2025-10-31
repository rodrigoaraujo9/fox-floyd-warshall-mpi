#!/usr/bin/env python3
import csv
import os
from collections import defaultdict
import statistics

results_dir = "results"
report_dir = "report_data"
os.makedirs(report_dir, exist_ok=True)

# Store all results
data = defaultdict(lambda: defaultdict(list))
read_errors = []
ingested = 0

def is_mac_resource_fork(name: str) -> bool:
    return name.startswith("._")

# Read all CSV files
for filename in os.listdir(results_dir):
    if not filename.endswith('.csv') or is_mac_resource_fork(filename):
        continue
    filepath = os.path.join(results_dir, filename)
    try:
        with open(filepath, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                alg = row['algorithm']
                n = int(row['matrix_size'])
                time_sec = float(row['time_sec'])
                p = int(row.get('processes', 1))
                key = f"{alg}_n{n}_p{p}"
                data[key]['times'].append(time_sec * 1000.0)  # ms
                data[key]['n'] = n
                data[key]['p'] = p
                data[key]['alg'] = alg
                ingested += 1
    except Exception as e:
        read_errors.append(f"{filename}: {e}")

# Aggregate
results = {}
for key, values in data.items():
    if values['times']:
        results[key] = {
            'n': values['n'],
            'p': values['p'],
            'alg': values['alg'],
            'median': statistics.median(values['times']),
            'mean': statistics.mean(values['times']),
            'stdev': statistics.stdev(values['times']) if len(values['times']) > 1 else 0.0,
        }

print(f"# Extract Data Summary")
print(f"# Entries ingested: { ingested }")
if read_errors:
    print(f"# Warnings while reading result CSVs ({len(read_errors)}):")
    for w in read_errors[:10]:
        print(f"#  - {w}")
    if len(read_errors) > 10:
        print(f"#  ... {len(read_errors)-10} more")

graph_names = {300: 'A', 600: 'B', 900: 'C', 1200: 'D'}

def get_res(alg, n, p):
    return results.get(f"{alg}_n{n}_p{p}")

# === TABLE 1: Sequential Baseline Performance ===
print("\n=== TABLE 1: Sequential Baseline Performance ===")
print("Algorithm & Graph A (300) & Graph B (600) & Graph C (900) & Graph D (1200) \\\\")
print("\\hline")

algs = [('rs', 'Repeated Squaring'),
        ('fw', 'Standard Floyd-Warshall'),
        ('blocked', 'Blocked Floyd-Warshall')]

for alg_key, alg_name in algs:
    row = [alg_name]
    for n in [300, 600, 900, 1200]:
        r = get_res(alg_key, n, 1)
        row.append(f"{r['median']:.1f}" if r else "—")
    print(" & ".join(row) + " \\\\")

# === TABLE 2: Parallel Performance - Blocking ===
print("\n=== TABLE 2: Parallel Performance - Blocking ===")
print("Graph & P=1 & P=4 & P=9 & P=16 & P=25 \\\\")
print("\\hline")

for n in [300, 600, 900, 1200]:
    row = [f"{graph_names[n]} ({n})"]
    for p in [1, 4, 9, 16, 25]:
        r = get_res("mpi", n, p)
        row.append(f"{r['median']:.1f}" if r else "—")
    print(" & ".join(row) + " \\\\")

# === TABLE 3: Parallel Performance - Non-Blocking ===
print("\n=== TABLE 3: Parallel Performance - Non-Blocking ===")
print("Graph & P=1 & P=4 & P=9 & P=16 & P=25 \\\\")
print("\\hline")

for n in [300, 600, 900, 1200]:
    row = [f"{graph_names[n]} ({n})"]
    for p in [1, 4, 9, 16, 25]:
        r = get_res("mpi-nb", n, p)
        row.append(f"{r['median']:.1f}" if r else "—")
    print(" & ".join(row) + " \\\\")

# === TABLE 4: Speedup Relative to Sequential Blocked ===
print("\n=== TABLE 4: Speedup Relative to Sequential Blocked ===")
print("Graph & P=1 & P=4 & P=9 & P=16 & P=25 \\\\")
print("\\hline")

for n in [300, 600, 900, 1200]:
    seq = get_res("blocked", n, 1)
    if not seq:
        continue
    seq_time = seq['median']
    row = [f"{graph_names[n]} ({n})"]
    for p in [1, 4, 9, 16, 25]:
        r = get_res("mpi", n, p)
        row.append(f"{(seq_time / r['median']):.2f}" if r else "—")
    print(" & ".join(row) + " \\\\")

# === TABLE 5: Parallel Efficiency ===
print("\n=== TABLE 5: Parallel Efficiency ===")
print("Graph & P=4 & P=9 & P=16 & P=25 \\\\")
print("\\hline")

for n in [300, 600, 900, 1200]:
    seq = get_res("blocked", n, 1)
    if not seq:
        continue
    seq_time = seq['median']
    row = [f"{graph_names[n]} ({n})"]
    for p in [4, 9, 16, 25]:
        r = get_res("mpi", n, p)
        if r:
            speedup = seq_time / r['median']
            efficiency = (speedup / p) * 100.0
            row.append(f"{efficiency:.1f}\\%")
        else:
            row.append("—")
    print(" & ".join(row) + " \\\\")

# === TABLE 6: Communication Overhead (% of total time) ===
# Overhead% = 100 - efficiency%
print("\n=== TABLE 6: Communication Overhead (% of total time) ===")
print("Graph & P=4 & P=9 & P=16 & P=25 \\\\")
print("\\hline")

for n in [300, 600, 900, 1200]:
    seq = get_res("blocked", n, 1)
    if not seq:
        continue
    seq_time = seq['median']
    row = [f"{graph_names[n]} ({n})"]
    for p in [4, 9, 16, 25]:
        r = get_res("mpi", n, p)
        if r:
            speedup = seq_time / r['median']
            efficiency = (speedup / p) * 100.0
            overhead = max(0.0, 100.0 - efficiency)
            row.append(f"{overhead:.1f}\\%")
        else:
            row.append("—")
    print(" & ".join(row) + " \\\\")

# === CSV for plotting ===
out_csv = os.path.join(report_dir, 'speedup_data.csv')
with open(out_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['graph', 'n', 'p', 'speedup_vs_seq', 'speedup_vs_p1', 'efficiency'])
    for n in [300, 600, 900, 1200]:
        seq = get_res("blocked", n, 1)
        p1  = get_res("mpi", n, 1)
        if not (seq and p1):
            continue
        seq_time = seq['median']
        p1_time  = p1['median']
        for p in [1, 4, 9, 16, 25]:
            r = get_res("mpi", n, p)
            if not r:
                continue
            par_time = r['median']
            speedup_seq = seq_time / par_time
            speedup_p1  = p1_time  / par_time
            efficiency  = (speedup_seq / p) * 100.0
            writer.writerow([
                graph_names[n], n, p,
                f"{speedup_seq:.3f}",
                f"{speedup_p1:.3f}",
                f"{efficiency:.2f}"
            ])

print(f"\nData saved to {out_csv}")
print("\nAll tables generated successfully!")
