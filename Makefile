CFLAGS=-std=c++17
OPT=-O2


all:
	g++ -o mg_constructor test_multigrid_class.cpp $(CFLAGS)
	./mg_constructor


test:
	g++ -o multigrid -DMULTIGRID test_multigrid.cpp $(CFLAGS) $(OPT)
	./multigrid > twolvl.txt
	g++ -o multigrid test_multigrid.cpp $(CFLAGS) $(OPT)
	./multigrid > gseidel.txt
	octave --persist report.m

coarse:
	g++ -o coarsening test_coarsening.cpp $(CFLAGS)


convergence:
	g++ -o convergence test_convergence.cpp $(CFLAGS)
	./convergence

parallelism:
    g++ -fopenmp -o parallel test_parallelism.cpp -O2
	./parallel

