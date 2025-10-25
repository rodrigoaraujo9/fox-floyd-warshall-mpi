CC      := mpicc
CFLAGS  := -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
INCLUDE := -Iinclude
SRC_COMMON := src/io.c src/matrix.c src/algorithms.c src/comm.c src/blocked_fw.c
BIN     := apsp

all: $(BIN)

$(BIN): src/main.c $(SRC_COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ -lm

clean:
	rm -f $(BIN) $(BIN_RS)

run600-4:
	mpirun -np 4 ./apsp matrix_examples/input600 matrix_examples/output600

run1200-9:
	mpirun -np 9 --oversubscribe ./apsp matrix_examples/input1200 matrix_examples/output1200
