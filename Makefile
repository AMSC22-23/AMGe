CFLAGS=-std=c++17


all:
	g++ -o test main.cpp $(CFLAGS)


coarse:
	g++ -o coarsening test_coarsening.cpp $(CFLAGS)


multigrid:
	g++ -o multigrid test_multigrid.cpp $(CFLAGS)
