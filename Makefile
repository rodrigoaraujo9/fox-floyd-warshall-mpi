CC      := mpicc
CFLAGS  := -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
INCLUDE := -Iinclude
SRC_COMMON := src/io.c src/matrix.c src/algorithms.c src/blocked_fw.c src/blocked_rs.c
BIN     := apsp
BIN_RS  := apsp_rs

all: $(BIN)

$(BIN): src/main.c $(SRC_COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ -lm

$(BIN_RS): src/main_rs.c $(SRC_COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ -lm

clean:
	rm -f $(BIN) $(BIN_RS)

run600-9: $(BIN)
	mpirun -np 9 ./$(BIN) matrix_examples/input600 matrix_examples/output600

run600-9-rs: $(BIN_RS)
	mpirun -np 4 ./$(BIN_RS) matrix_examples/input600 matrix_examples/output600
