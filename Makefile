CC      := mpicc
CFLAGS  := -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
INCLUDE := -Iinclude
SRC     := src/main.c src/io.c src/matrix.c src/algorithms.c src/blocked_fw.c
BIN     := apsp

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) -lm

clean:
	rm -f $(BIN)

run600-4:
	mpirun -np 4 ./apsp matrix_examples/input600 matrix_examples/output600
