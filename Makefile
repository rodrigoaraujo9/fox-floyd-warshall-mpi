# ===== Compiler & Build =====
CC       := mpicc
CFLAGS   := -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
INCLUDE  := -Iinclude
SRC_COMMON := src/io.c src/matrix.c src/algorithms.c src/comm.c src/blocked_fw.c
BIN      := apsp

# ===== MPI Runtime (tuned for the labs) =====
# Choose hostfile at run time: `make HOSTFILE=hostfile_100 test-mpi-1200-25`
HOSTFILE ?= hostfile_100

# Quiet + robust OpenMPI/PMIx defaults
MPI_ENV   := PMIX_MCA_gds=hash OMPI_MCA_pmix_base_async_modex=0
MCA_FLAGS := --mca plm_rsh_num_concurrent 8 --mca plm_rsh_no_tree_spawn 1 --mca btl self,tcp

# mpirun wrapper (unset DISPLAY to silence X auth noise)
MPIRUN := env -u DISPLAY $(MPI_ENV) mpirun $(MCA_FLAGS) --hostfile $(HOSTFILE) --map-by node

# ===== Python / Reporting =====
PYTHON   ?= python3
RESULTS_DIR := results
REPORT_DIR  := report_data

# ===== Default =====
all: $(BIN)

# ===== Build binary =====
$(BIN): src/main.c $(SRC_COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ -lm

# ===== Cleanup =====
clean:
	rm -f $(BIN)

distclean: clean
	rm -rf $(RESULTS_DIR) $(REPORT_DIR)

results-clean:
	rm -rf $(RESULTS_DIR) && mkdir -p $(RESULTS_DIR)

report-clean:
	rm -rf $(REPORT_DIR) && mkdir -p $(REPORT_DIR)

# ===== Small tests (sequential + MPI) =====
test-fw-small: $(BIN)
	./$(BIN) fw matrix_examples/input5 matrix_examples/output5

test-rs-small: $(BIN)
	./$(BIN) rs matrix_examples/input5 matrix_examples/output5

test-blocked-small: $(BIN)
	./$(BIN) blocked matrix_examples/input6 matrix_examples/output6

test-mpi-small: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi matrix_examples/input6 matrix_examples/output6

test-mpi-nb-small: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi-nb matrix_examples/input6 matrix_examples/output6

# ===== Medium (300 / 600) =====
test-fw-300: $(BIN)
	./$(BIN) fw matrix_examples/input300 matrix_examples/output300

test-rs-300: $(BIN)
	./$(BIN) rs matrix_examples/input300 matrix_examples/output300

test-mpi-300-4: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi matrix_examples/input300 matrix_examples/output300

test-mpi-nb-300-4: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi-nb matrix_examples/input300 matrix_examples/output300

test-fw-600: $(BIN)
	./$(BIN) fw matrix_examples/input600 matrix_examples/output600

test-rs-600: $(BIN)
	./$(BIN) rs matrix_examples/input600 matrix_examples/output600

test-mpi-600-4: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi matrix_examples/input600 matrix_examples/output600

test-mpi-600-9: $(BIN)
	$(MPIRUN) -np 9 ./$(BIN) mpi matrix_examples/input600 matrix_examples/output600

test-mpi-nb-600-4: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi-nb matrix_examples/input600 matrix_examples/output600

test-mpi-nb-600-9: $(BIN)
	$(MPIRUN) -np 9 ./$(BIN) mpi-nb matrix_examples/input600 matrix_examples/output600

# ===== Large (900 / 1200) =====
test-fw-900: $(BIN)
	./$(BIN) fw matrix_examples/input900 matrix_examples/output900

test-mpi-900-9: $(BIN)
	$(MPIRUN) -np 9 ./$(BIN) mpi matrix_examples/input900 matrix_examples/output900

test-mpi-900-16: $(BIN)
	$(MPIRUN) -np 16 ./$(BIN) mpi matrix_examples/input900 matrix_examples/output900

test-mpi-nb-900-9: $(BIN)
	$(MPIRUN) -np 9 ./$(BIN) mpi-nb matrix_examples/input900 matrix_examples/output900

test-mpi-nb-900-16: $(BIN)
	$(MPIRUN) -np 16 ./$(BIN) mpi-nb matrix_examples/input900 matrix_examples/output900

test-fw-1200: $(BIN)
	./$(BIN) fw matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-4: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-9: $(BIN)
	$(MPIRUN) -np 9 ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-16: $(BIN)
	$(MPIRUN) -np 16 ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-25: $(BIN)
	$(MPIRUN) -np 25 ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-100: $(BIN)
	$(MPIRUN) -np 100 ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-4: $(BIN)
	$(MPIRUN) -np 4 ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-9: $(BIN)
	$(MPIRUN) -np 9 ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-16: $(BIN)
	$(MPIRUN) -np 16 ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-25: $(BIN)
	$(MPIRUN) -np 25 ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-100: $(BIN)
	$(MPIRUN) -np 100 ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200

# ===== Compare groups =====
compare-300: test-fw-300 test-rs-300 test-mpi-300-4 test-mpi-nb-300-4

compare-600: test-fw-600 test-rs-600 test-mpi-600-4 test-mpi-600-9 test-mpi-nb-600-4 test-mpi-nb-600-9

compare-900: test-fw-900 test-mpi-900-9 test-mpi-900-16 test-mpi-nb-900-9 test-mpi-nb-900-16

compare-1200: test-fw-1200 test-mpi-1200-4 test-mpi-1200-9 test-mpi-1200-16 test-mpi-1200-25 \
              test-mpi-nb-1200-4 test-mpi-nb-1200-9 test-mpi-nb-1200-16

# ===== Scaling suites =====
scale-600: test-mpi-600-4 test-mpi-600-9 test-mpi-nb-600-4 test-mpi-nb-600-9

scale-1200: test-mpi-1200-4 test-mpi-nb-1200-4 test-mpi-1200-9 test-mpi-nb-1200-9 \
            test-mpi-1200-16 test-mpi-nb-1200-16 test-mpi-1200-25 test-mpi-nb-1200-25 \
            test-mpi-1200-100 test-mpi-nb-1200-100

# ===== Run everything =====
test-all-small: test-fw-small test-rs-small test-blocked-small test-mpi-small test-mpi-nb-small
test-all-medium: compare-300 compare-600
test-all-large:  compare-900 compare-1200
test-all: test-all-small test-all-medium test-all-large

# ===== Benchmark / Report pipeline =====
# Run the big benchmark script (honors HOSTFILE)
bench: $(BIN)
	HOSTFILE=$(HOSTFILE) bash ./comprehensive_benchmark.sh

# Extract tables from CSVs
tables: $(REPORT_DIR)/tables.txt
$(REPORT_DIR)/tables.txt: extract_data.py
	mkdir -p $(REPORT_DIR)
	$(PYTHON) extract_data.py > $(REPORT_DIR)/tables.txt

# Generate plots (reads report_data/speedup_data.csv)
plots: generate_plots.py
	mkdir -p $(REPORT_DIR)
	$(PYTHON) generate_plots.py

# Convenience target to prep everything for LaTeX
report: tables plots
	@echo "Tables at:   $(REPORT_DIR)/tables.txt"
	@echo "Figures at:  $(REPORT_DIR)/speedup_vs_sequential.png"
	@echo "             $(REPORT_DIR)/speedup_vs_p1.png"
	@echo "             $(REPORT_DIR)/strong_scaling.png"
	@echo "             $(REPORT_DIR)/parallel_efficiency.png"
	@echo "             $(REPORT_DIR)/speedup_comparison.png"

# ===== Python dependencies quick check =====
deps-check:
	$(PYTHON) - <<'PY'
import sys, importlib
missing=[]
for p in ("pandas","matplotlib","numpy"):
    try: importlib.import_module(p)
    except ImportError: missing.append(p)
if missing:
    sys.exit("Missing: "+", ".join(missing)+"  -> run: pip3 install --user "+ " ".join(missing))
print("Python deps OK")
PY

# ===== Helpers =====
slots:
	@awk '{for(i=1;i<=NF;i++) if($$i ~ /^slots=/){split($$i,a,"="); s+=a[2]}} END{print "TOTAL SLOTS:", s+0}' $(HOSTFILE)

who:
	@$(MPIRUN) -np 8 hostname -s | sort | uniq -c

.PHONY: all clean distclean results-clean report-clean \
        test-fw-small test-rs-small test-blocked-small test-mpi-small test-mpi-nb-small \
        test-fw-300 test-rs-300 test-mpi-300-4 test-mpi-nb-300-4 \
        test-fw-600 test-rs-600 test-mpi-600-4 test-mpi-600-9 test-mpi-nb-600-4 test-mpi-nb-600-9 \
        test-fw-900 test-mpi-900-9 test-mpi-900-16 test-mpi-nb-900-9 test-mpi-nb-900-16 \
        test-fw-1200 test-mpi-1200-4 test-mpi-1200-9 test-mpi-1200-16 test-mpi-1200-25 test-mpi-1200-100 \
        test-mpi-nb-1200-4 test-mpi-nb-1200-9 test-mpi-nb-1200-16 test-mpi-nb-1200-25 test-mpi-nb-1200-100 \
        compare-300 compare-600 compare-900 compare-1200 \
        scale-600 scale-1200 test-all-small test-all-medium test-all-large test-all \
        bench tables plots report deps-check slots who
