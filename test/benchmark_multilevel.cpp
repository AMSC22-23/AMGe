#include <benchmark/benchmark.h>
#include <iostream>
#include <vector>
#include "Lattice.hpp"
#include "Multigrid.hpp"


#define PRE_SMOOTHING_STEPS 10
#define POST_SMOOTHING_STEPS 10


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


static void BM_TwoLevel(benchmark::State& state) {
	int n = state.range(0);
	Lattice fine(0.0, 0.0, 1.0, 1.0, n+1);
	Multigrid solver(fine, PRE_SMOOTHING_STEPS, POST_SMOOTHING_STEPS, 2);


	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());


	fine.evaluate_boundary_conditions(u, g);
	fine.evaluate_forcing_term(b, f);


	for (auto _ : state) {
		solver.step(u, b);
	}
}


static void BM_ThreeLevel(benchmark::State& state) {
	int n = state.range(0);
	Lattice fine(0.0, 0.0, 1.0, 1.0, n+1);
	Multigrid solver(fine, PRE_SMOOTHING_STEPS, POST_SMOOTHING_STEPS, 3);


	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());


	fine.evaluate_boundary_conditions(u, g);
	fine.evaluate_forcing_term(b, f);


	for (auto _ : state) {
		solver.step(u, b);
	}
}


static void BM_FourLevel(benchmark::State& state) {
	int n = state.range(0);
	Lattice fine(0.0, 0.0, 1.0, 1.0, n+1);
	Multigrid solver(fine, PRE_SMOOTHING_STEPS, POST_SMOOTHING_STEPS, 4);


	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());


	fine.evaluate_boundary_conditions(u, g);
	fine.evaluate_forcing_term(b, f);


	for (auto _ : state) {
		solver.step(u, b);
	}
}


BENCHMARK(BM_TwoLevel)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK(BM_ThreeLevel)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK(BM_FourLevel)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK_MAIN();
