BUILD_DIR = .
INC_DIR = ../include
SRC_DIR = ../src
TEST_DIR = ../test


CXX = g++
CXXFLAGS += -std=c++17
INCFLAGS = -I $(INC_DIR)
OPTFLAGS = -O2
PARALLELFLAGS = -fopenmp
BENCHFLAGS = -lbenchmark -lpthread


# hack
CXXFLAGS += $(PARALLELFLAGS)


SRCS      = $(wildcard $(SRC_DIR)/*.cpp)
SRCS_TEST = $(wildcard $(TEST_DIR)/*.cpp)
OBJS      = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))
OBJS_TEST = $(patsubst $(TEST_DIR)/%.cpp,$(BUILD_DIR),$(SRCS_TEST))


all: main


main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -o $@ ../main.cpp $^


$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -c -o $@ $<


test: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -o test_coarsening $(TEST_DIR)/test_coarsening.cpp $^
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -o test_convergence $(TEST_DIR)/test_convergence.cpp $^
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -o test_convergence_levels $(TEST_DIR)/test_convergence_levels.cpp $^
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -o bench_multilevel $(TEST_DIR)/benchmark_multilevel.cpp $^ $(BENCHFLAGS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(OPTFLAGS) -o bench_jacobi $(TEST_DIR)/benchmark_jacobi.cpp $^ $(BENCHFLAGS)


benchmark:
	./bench_jacobi --benchmark_format=csv > jacobi.csv
	python plot.py -f jacobi.csv


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
