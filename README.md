    	mpicc main.c && mpirun -np 8 --oversubscribe a.out

    	https://moorejs.github.io/APSP-in-parallel/

     changed comm blocks to be contiguous and have less comm calls cause they are expensive timewise.
