# define the C/C++ compiler to use,default here is clang
CC = gcc-5

all:
	$(CC) -o sequential vptree_sequential.c main.c -lm;
	./sequential
	$(CC) -o pthreads vptree_pthreads.c main.c -lm -pthread;
	./pthreads
	$(CC) -o cilk vptree_cilk.c main.c -lm -fcilkplus;
	./cilk
	$(CC) -o openmp vptree_openmp.c main.c -lm -fopenmp;
	./openmp
