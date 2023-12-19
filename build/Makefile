BUILD_DIR = .
INC_DIR = ../include
SRC_DIR = ../src
TST_DIR = ../test


CXX = g++
CXXFLAGS += -std=c++17
INCFLAGS = -I $(INC_DIR)
OPTFLAGS = -O2


SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))


all: main


main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -o $@ $^


$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -c -o $@ $<


clean:
	rm -f main
	rm -f *.o


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
