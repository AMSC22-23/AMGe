CFLAGS=-std=c++17
OPT=-O2


all:
	g++ -o multigrid -DMULTIGRID test_multigrid.cpp $(CFLAGS) $(OPT)
	./multigrid > twolvl.txt
	g++ -o multigrid test_multigrid.cpp $(CFLAGS) $(OPT)
	./multigrid > gseidel.txt
	octave --persist report.m


multigrid_test:
	g++ -o mg_constructor test_multigrid_class.cpp $(CFLAGS)
	./mg_constructor


coarse:
	g++ -o coarsening test_coarsening.cpp $(CFLAGS)


convergence:
	g++ -o convergence test_convergence.cpp $(CFLAGS)
	./convergence


parallelism:
	g++ -fopenmp -o parallel test_parallelism.cpp -O2 $(CFLAGS)
	./parallel

benchmark:
	g++ -o bench test_jacobi_parallelism.cpp $(CFLAGS) -lbenchmark -lpthread -fopenmp
	./bench
