#!/bin/bash
set -euo pipefail

# Comprehensive benchmark script for APSP Performance Evaluation
# Runs all tests needed to populate the tables and figures in the report

ITERATIONS=5
RESULTS_DIR="results"
REPORT_DIR="report_data"

mkdir -p "$RESULTS_DIR" "$REPORT_DIR"

echo "=========================================="
echo "  APSP Comprehensive Performance Testing"
echo "=========================================="
echo ""
echo "This will run all tests needed for:"
echo "  - Table: Sequential Performance (Table 1)"
echo "  - Table: Blocking MPI Performance (Table 2)"
echo "  - Table: Non-blocking MPI Performance (Table 3)"
echo "  - Table: Speedup Analysis (Table 4)"
echo "  - Table: Parallel Efficiency (Table 5)"
echo "  - Table: Communication Overhead (Table 6)"
echo "  - Figure: Speedup vs Sequential"
echo "  - Figure: Speedup vs P=1"
echo "  - Figure: Strong Scaling"
echo ""
echo "Total iterations per configuration: $ITERATIONS"
echo ""

# Helper: integer sqrt + checks (P perfect square, N divisible by sqrt(P))
is_perfect_square_and_divisible() {
  local P="$1" N="$2"
  local Q
  Q=$(printf "%d" "$(echo "scale=12; sqrt($P)" | bc -l)")
  # Must be perfect square
  if [ $((Q*Q)) -ne "$P" ]; then
    return 1
  fi
  # And N must be divisible by Q
  if [ $((N % Q)) -ne 0 ]; then
    return 1
  fi
  echo "$Q"
  return 0
}

# =============================================================================
# SECTION 1: SEQUENTIAL BASELINE PERFORMANCE (Table 1)
# =============================================================================
echo "=========================================="
echo "SECTION 1: Sequential Baseline Tests"
echo "=========================================="

GRAPH_SIZES=(300 600 900 1200)
GRAPH_NAMES=("A" "B" "C" "D")

for idx in "${!GRAPH_SIZES[@]}"; do
  N=${GRAPH_SIZES[$idx]}
  GRAPH=${GRAPH_NAMES[$idx]}

  echo ""
  echo "Testing Graph $GRAPH (N=$N)"
  echo "----------------------------"

  # Repeated Squaring
  echo "  Running Repeated Squaring..."
  for i in $(seq 1 "$ITERATIONS"); do
    ./apsp rs "matrix_examples/input$N" "matrix_examples/output$N" > /dev/null 2>&1
  done

  # Standard Floyd-Warshall
  echo "  Running Standard Floyd-Warshall..."
  for i in $(seq 1 "$ITERATIONS"); do
    ./apsp fw "matrix_examples/input$N" "matrix_examples/output$N" > /dev/null 2>&1
  done

  # Blocked Floyd-Warshall (sequential)
  echo "  Running Blocked Floyd-Warshall..."
  for i in $(seq 1 "$ITERATIONS"); do
    ./apsp blocked "matrix_examples/input$N" "matrix_examples/output$N" > /dev/null 2>&1
  done
done

# =============================================================================
# SECTION 2: PARALLEL PERFORMANCE - BLOCKING (Table 2)
# =============================================================================
echo ""
echo "=========================================="
echo "SECTION 2: Parallel Blocking MPI Tests"
echo "=========================================="

PROCESS_COUNTS=(1 4 9 16 25)

for idx in "${!GRAPH_SIZES[@]}"; do
  N=${GRAPH_SIZES[$idx]}
  GRAPH=${GRAPH_NAMES[$idx]}

  echo ""
  echo "Testing Graph $GRAPH (N=$N)"
  echo "----------------------------"

  for P in "${PROCESS_COUNTS[@]}"; do
    Q=$(is_perfect_square_and_divisible "$P" "$N") || continue
    echo "  Running MPI Blocking with P=$P..."
    for i in $(seq 1 "$ITERATIONS"); do
      mpirun -np "$P" --hostfile hostfile --map-by node \
        ./apsp mpi "matrix_examples/input$N" "matrix_examples/output$N" > /dev/null 2>&1
    done
  done
done

# =============================================================================
# SECTION 3: PARALLEL PERFORMANCE - NON-BLOCKING (Table 3)
# =============================================================================
echo ""
echo "=========================================="
echo "SECTION 3: Parallel Non-Blocking MPI Tests"
echo "=========================================="

for idx in "${!GRAPH_SIZES[@]}"; do
  N=${GRAPH_SIZES[$idx]}
  GRAPH=${GRAPH_NAMES[$idx]}

  echo ""
  echo "Testing Graph $GRAPH (N=$N)"
  echo "----------------------------"

  for P in "${PROCESS_COUNTS[@]}"; do
    Q=$(is_perfect_square_and_divisible "$P" "$N") || continue
    echo "  Running MPI Non-Blocking with P=$P..."
    for i in $(seq 1 "$ITERATIONS"); do
      mpirun -np "$P" --hostfile hostfile --map-by node \
        ./apsp mpi-nb "matrix_examples/input$N" "matrix_examples/output$N" > /dev/null 2>&1
    done
  done
done

# =============================================================================
# SECTION 4: DATA EXTRACTION AND REPORT GENERATION
# =============================================================================
echo ""
echo "=========================================="
echo "SECTION 4: Extracting Results"
echo "=========================================="

cat > extract_data.py << 'PYTHON_SCRIPT'
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
PYTHON_SCRIPT

chmod +x extract_data.py
python3 extract_data.py > "$REPORT_DIR/tables.txt"

echo ""
echo "=========================================="
echo "Testing Complete!"
echo "=========================================="
echo ""
echo "Results saved to:"
echo "  - Raw data: $RESULTS_DIR/"
echo "  - Report tables: $REPORT_DIR/tables.txt"
echo "  - Plot data: $REPORT_DIR/speedup_data.csv"
echo ""
echo "Next steps:"
echo "  1. Review $REPORT_DIR/tables.txt for LaTeX table data"
echo "  2. Use $REPORT_DIR/speedup_data.csv to generate plots (see generate_plots.py)"
echo "  3. Update your report with the actual measurements"
echo ""
