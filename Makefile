CFLAGS=-std=c++17


all:
	g++ -o main main.cpp $(CFLAGS)
	./main

coarse:
	g++ -o coarsening test_coarsening.cpp $(CFLAGS)


multigrid:
	g++ -o multigrid test_multigrid.cpp $(CFLAGS)


convergence:
	g++ -o convergence test_convergence.cpp $(CFLAGS)
	./convergence
