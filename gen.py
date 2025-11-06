import random
import numpy as np
import subprocess
import os


def generate_sparse_graph(n=3000, edge_probability=0.01, max_weight=100, seed=42):
    """
    Generate a sparse directed weighted graph for APSP testing.

    Args:
        n: Matrix dimension
        edge_probability: Probability of edge existence between nodes
        max_weight: Maximum edge weight (weights are 1 to max_weight)
        seed: Random seed for reproducibility
    """
    random.seed(seed)
    np.random.seed(seed)

    # Initialize adjacency matrix with zeros (no edges)
    matrix = np.zeros((n, n), dtype=int)

    # Set diagonal to 0 (distance to self)
    for i in range(n):
        matrix[i][i] = 0

    # Add random edges
    for i in range(n):
        for j in range(n):
            if i != j and random.random() < edge_probability:
                matrix[i][j] = random.randint(1, max_weight)

    # Ensure graph is connected by creating a path through all vertices
    for i in range(n - 1):
        if matrix[i][i + 1] == 0:
            matrix[i][i + 1] = random.randint(1, max_weight)

    # Add some random backward edges for cycles
    for _ in range(n // 10):
        i = random.randint(1, n - 1)
        j = random.randint(0, i - 1)
        matrix[i][j] = random.randint(1, max_weight)

    return matrix


def save_matrix(matrix, filename):
    """Save matrix to file in the required format."""
    n = len(matrix)
    with open(filename, "w") as f:
        f.write(f"{n}\n")
        for i in range(n):
            f.write(" ".join(map(str, matrix[i])) + "\n")


def save_matrix_no_header(matrix, filename):
    """Save matrix without dimension header (for output file)."""
    with open(filename, "w") as f:
        for i in range(len(matrix)):
            f.write(" ".join(map(str, matrix[i])) + "\n")


def run_parallel_fw(input_file, output_file, num_processes=4):
    """Run parallel Floyd-Warshall using MPI to generate output."""
    print(f"\nRunning parallel Floyd-Warshall with {num_processes} processes...")
    print("This will compute the shortest paths using your MPI implementation.")

    # Create a temporary dummy output file for verification
    temp_output = "temp_output_3000.txt"

    # Create dummy output (all zeros) - the algorithm doesn't actually verify during computation
    # We'll use the result it produces as the expected output
    n = 3000
    with open(temp_output, "w") as f:
        for i in range(n):
            f.write(" ".join(["0"] * n) + "\n")

    # Run MPI program to compute the actual result
    cmd = f"mpirun -np {num_processes} ./bin/apsp mpi {input_file} {temp_output}"
    print(f"Executing: {cmd}")

    try:
        result = subprocess.run(
            cmd, shell=True, check=True, capture_output=True, text=True
        )
        print(result.stdout)
        if result.stderr:
            print("Stderr:", result.stderr)
        print("Computation completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running MPI: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    finally:
        if os.path.exists(temp_output):
            os.remove(temp_output)


# Generate input graph
print("Generating 3000x3000 sparse graph...")
input_matrix = generate_sparse_graph(n=3000, edge_probability=0.01, max_weight=100)

print("Saving input matrix...")
save_matrix(input_matrix, "input_3000.txt")

print("\nGraph statistics:")
print(f"  Nodes: 3000")
print(f"  Edges: {np.count_nonzero(input_matrix)}")
print(f"  Density: {np.count_nonzero(input_matrix) / (3000 * 3000) * 100:.2f}%")

print("\n" + "=" * 60)
print("NEXT STEP: Run your parallel Floyd-Warshall to generate output")
print("=" * 60)
print("\nOption 1 - Manual (recommended):")
print("  1. Compile your code: make")
print("  2. Run with MPI: mpirun -np 4 ./bin/apsp mpi input_3000.txt output_3000.txt")
print("     (Use any number of processes that divides sqrt(3000) evenly)")
print("\nOption 2 - Automatic (uncomment code below):")
print("  Uncomment the run_parallel_fw() call in the script")
print("\nNote: Make sure your code writes the result matrix to a file!")
print("=" * 60)

# Uncomment to run automatically:
# run_parallel_fw('input_3000.txt', 'output_3000.txt', num_processes=4)

print("\nInput file created: input_3000.txt")
