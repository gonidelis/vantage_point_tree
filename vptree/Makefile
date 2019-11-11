# define the C/C++ compiler to use,default here is clang
CC = gcc-5

.PHONY: lib

lib:
	cd src; $(CC) -c vptree_sequential.c; cd ..
	cd src; $(CC) -c vptree_pthreads.c; cd ..
	cd src; $(CC) -c -fcilkplus vptree_cilk.c; cd ..
	cd src; $(CC) -c vptree_openmp.c; cd ..
	cd src; ar rcs ../lib/vptree_sequential.a vptree_sequential.o; cd ..
	cd src; ar rcs ../lib/vptree_pthreads.a vptree_pthreads.o; cd ..
	cd src; ar rcs ../lib/vptree_cilk.a vptree_cilk.o; cd ..
	cd src; ar rcs ../lib/vptree_openmp.a vptree_openmp.o; cd ..
