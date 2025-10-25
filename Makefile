CC      := mpicc
CFLAGS  := -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
INCLUDE := -Iinclude
SRC_COMMON := src/io.c src/matrix.c src/algorithms.c src/comm.c src/blocked_fw.c
BIN     := apsp

# Default target
all: $(BIN)

# Build binary
$(BIN): src/main.c $(SRC_COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ -lm

clean:
	rm -f $(BIN)

#small
test-fw-small:
	./$(BIN) fw matrix_examples/input5 matrix_examples/output5

test-rs-small:
	./$(BIN) rs matrix_examples/input5 matrix_examples/output5

test-blocked-small:
	./$(BIN) blocked matrix_examples/input6 matrix_examples/output6

test-mpi-small:
	mpirun -np 4 ./$(BIN) mpi matrix_examples/input6 matrix_examples/output6

test-mpi-nb-small:
	mpirun -np 4 ./$(BIN) mpi-nb matrix_examples/input6 matrix_examples/output6


#medium
test-fw-300:
	./$(BIN) fw matrix_examples/input300 matrix_examples/output300

test-rs-300:
	./$(BIN) rs matrix_examples/input300 matrix_examples/output300

test-mpi-300-4:
	mpirun -np 4 ./$(BIN) mpi matrix_examples/input300 matrix_examples/output300

test-mpi-nb-300-4:
	mpirun -np 4 ./$(BIN) mpi-nb matrix_examples/input300 matrix_examples/output300

test-fw-600:
	./$(BIN) fw matrix_examples/input600 matrix_examples/output600

test-rs-600:
	./$(BIN) rs matrix_examples/input600 matrix_examples/output600

test-mpi-600-4:
	mpirun -np 4 ./$(BIN) mpi matrix_examples/input600 matrix_examples/output600

test-mpi-600-9:
	mpirun -np 9 --oversubscribe ./$(BIN) mpi matrix_examples/input600 matrix_examples/output600

test-mpi-nb-600-4:
	mpirun -np 4 ./$(BIN) mpi-nb matrix_examples/input600 matrix_examples/output600

test-mpi-nb-600-9:
	mpirun -np 9 --oversubscribe ./$(BIN) mpi-nb matrix_examples/input600 matrix_examples/output600



#large
test-fw-900:
	./$(BIN) fw matrix_examples/input900 matrix_examples/output900

test-mpi-900-9:
	mpirun -np 9 --oversubscribe ./$(BIN) mpi matrix_examples/input900 matrix_examples/output900

test-mpi-900-16:
	mpirun -np 16 --oversubscribe ./$(BIN) mpi matrix_examples/input900 matrix_examples/output900

test-mpi-nb-900-9:
	mpirun -np 9 --oversubscribe ./$(BIN) mpi-nb matrix_examples/input900 matrix_examples/output900

test-mpi-nb-900-16:
	mpirun -np 16 --oversubscribe ./$(BIN) mpi-nb matrix_examples/input900 matrix_examples/output900

test-fw-1200:
	./$(BIN) fw matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-4:
	mpirun -np 4 ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-9:
	mpirun -np 9 --oversubscribe ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-1200-16:
	mpirun -np 16 --oversubscribe ./$(BIN) mpi matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-4:
	mpirun -np 4 ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-9:
	mpirun -np 9 --oversubscribe ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200

test-mpi-nb-1200-16:
	mpirun -np 16 --oversubscribe ./$(BIN) mpi-nb matrix_examples/input1200 matrix_examples/output1200


#compare
compare-300: test-fw-300 test-rs-300 test-mpi-300-4 test-mpi-nb-300-4

compare-600: test-fw-600 test-rs-600 test-mpi-600-4 test-mpi-600-9 test-mpi-nb-600-4 test-mpi-nb-600-9

compare-900: test-fw-900 test-mpi-900-9 test-mpi-900-16 test-mpi-nb-900-9 test-mpi-nb-900-16

compare-1200: test-fw-1200 test-mpi-1200-4 test-mpi-1200-9 test-mpi-1200-16 test-mpi-nb-1200-4 test-mpi-nb-1200-9 test-mpi-nb-1200-16

#scale (processes)
scale-600: test-mpi-600-4 test-mpi-600-9 test-mpi-nb-600-4 test-mpi-nb-600-9

scale-1200: test-mpi-1200-4 test-mpi-1200-9 test-mpi-1200-16 test-mpi-nb-1200-4 test-mpi-nb-1200-9 test-mpi-nb-1200-16

#test all
test-all-small: test-fw-small test-rs-small test-blocked-small test-mpi-small test-mpi-nb-small

test-all-medium: compare-300 compare-600

test-all-large: compare-900 compare-1200

test-all: test-all-small test-all-medium test-all-large

.PHONY: all clean help \
        test-fw-small test-rs-small test-blocked-small test-mpi-small \
        test-fw-300 test-rs-300 test-mpi-300-4 \
        test-fw-600 test-rs-600 test-mpi-600-4 test-mpi-600-9 \
        test-fw-900 test-mpi-900-9 test-mpi-900-16 \
        test-fw-1200 test-mpi-1200-4 test-mpi-1200-9 test-mpi-1200-16 \
        compare-300 compare-600 compare-900 compare-1200 \
        scale-600 scale-1200 \
        test-all-small test-all-medium test-all-large test-all
