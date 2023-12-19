INCLUDE_DIR = ../include
SRC_DIR     = ../src
TEST_DIR    = ../test/


# compiler used
CXX = g++


CXXINC = -I $(INCLUDE_DIR)


# standard used
CXXFLAGS += -std=c++17


# warnings
# CXXWARNS = -Wall


# optimization flags
OPT_LVL = -O2


TARGET = main
TARGET_TEST = convergence_levels


all:
	g++ -o main ../src/main.cpp $(CXXFLAGS) -I../include


#
#all:
#	g++ -o multigrid test_multigrid.cpp $(CFLAGS) $(OPT)
#	./multigrid > gs.txt
#	g++ -o multigrid -DMULTIGRID2 test_multigrid.cpp $(CFLAGS) $(OPT)
#	./multigrid > due.txt
#	g++ -o multigrid -DMULTIGRID3 test_multigrid.cpp $(CFLAGS) $(OPT)
#	./multigrid > tre.txt
#	g++ -o multigrid -DMULTIGRID4 test_multigrid.cpp $(CFLAGS) $(OPT)
#	./multigrid > quattro.txt
#	octave --persist report.m
#
#
#multigrid_test:
#	g++ -o mg_constructor test_multigrid_class.cpp $(CFLAGS)
#	./mg_constructor
#
#
#coarse:
#	g++ -o coarsening test_coarsening.cpp $(CFLAGS)
#
#
#convergence:
#	g++ -o convergence test_convergence.cpp $(CFLAGS)
#	./convergence
#
#
#parallelism:
#	g++ -fopenmp -o parallel test_parallelism.cpp -O2 $(CFLAGS)
#	./parallel
#
#
#benchmark:
#	g++ -o bench test_jacobi_parallelism.cpp $(CFLAGS) -lbenchmark -lpthread -fopenmp
#	./bench
#
#
#benchmark_levels:
#	g++ -o bench_levels test_convergence_time_levels.cpp $(CFLAGS) $(OPT) -lbenchmark -lpthread -fopenmp
#	./bench
