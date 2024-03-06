CXX        = g++
CXXFLAGS   = -std=c++17
INCLUDE    = -I ./include
OPTFLAGS   = -O2
OMPFLAGS   = -fopenmp
BENCHFLAGS = -lbenchmark -lpthread


sources = $(wildcard src/*.cpp)
objects = $(patsubst src/%.cpp,build/%.o,$(sources))


test         = $(wildcard test/*.cpp)
test_objects = $(patsubst test/%.cpp,build/test/%.o,$(test))
test_targets = $(patsubst build/test/%.o,build/test/%,$(test_objects))


bench         = $(wildcard benchmark/*.cpp)
bench_targets = $(patsubst benchmark/%.cpp,build/benchmark/%,$(bench))


all: $(test_targets) $(bench_targets) controllo


controllo: build/test/test_multigrid2
	./$^ | tee build/temp_multigrid2.txt
	diff -s build/temp_multigrid2.txt build/reference_multigrid2.txt


build/test/%: build/test/%.o $(objects)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) $(OMPFLAGS) -o $@ $^


build/test/%.o: test/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) $(OMPFLAGS) -o $@ $^


build/bench/%: build/bench/%.o $(objects)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) $(OMPFLAGS) -o $@ $^ $(BENCHFLAGS)


build/bench/%.o: bench/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) $(OMPFLAGS) $(BENCHFLAGS) -o $@ $^


build/%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) -o $@ $^


.PHONY: folder
folder:
	mkdir -p build
	mkdir -p build/test
	mkdir -p build/bench


.PHONY: clean
clean:
	rm -f $(objects) $(test_objects) $(test_targets) $(bench_targets)
