CFLAGS=-std=c++17
OPT=-O2


all:
	g++ -o multigrid -DMULTIGRID test_multigrid.cpp $(CFLAGS) $(OPT)
	./multigrid > twolvl.txt
	g++ -o multigrid test_multigrid.cpp $(CFLAGS) $(OPT)
	./multigrid > gseidel.txt

coarse:
	g++ -o coarsening test_coarsening.cpp $(CFLAGS)


convergence:
	g++ -o convergence test_convergence.cpp $(CFLAGS)
	./convergence
