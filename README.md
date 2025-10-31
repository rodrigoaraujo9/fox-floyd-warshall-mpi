    	mpicc main.c && mpirun -np 8 --oversubscribe a.out

     changed comm blocks to be contiguous and have less comm calls cause they are expensive timewise.
