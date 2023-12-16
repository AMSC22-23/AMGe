#include <benchmark/benchmark.h>
#include <iostream>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"
#include "Smoothers.hpp"


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


static void BM_JacobiSerial(benchmark::State& state) {
	int n = state.range(0);
	Lattice fine(0.0, 0.0, 1.0, 1.0, n);


	std::vector<double> old(fine.numel());
	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());


	fine.evaluate_boundary_conditions(u, g);
	fine.evaluate_boundary_conditions(old, g);
	fine.evaluate_forcing_term(b, f);


	for (auto _ : state) {
		jacobi(fine, u, old, b);
	}
}


static void BM_JacobiParallelNaive(benchmark::State& state) {
	int n = state.range(0);
	Lattice fine(0.0, 0.0, 1.0, 1.0, n);


	std::vector<double> old(fine.numel());
	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());


	fine.evaluate_boundary_conditions(u, g);
	fine.evaluate_boundary_conditions(old, g);
	fine.evaluate_forcing_term(b, f);


	for (auto _ : state) {
		jacobi_parallel_naive(fine, u, old, b);
	}
}


static void BM_JacobiParallel(benchmark::State& state) {
	int n = state.range(0);
	Lattice fine(0.0, 0.0, 1.0, 1.0, n);


	std::vector<double> old(fine.numel());
	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());


	fine.evaluate_boundary_conditions(u, g);
	fine.evaluate_boundary_conditions(old, g);
	fine.evaluate_forcing_term(b, f);


	for (auto _ : state) {
		jacobi_parallel(fine, u, old, b);
	}
}


BENCHMARK(BM_JacobiSerial)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK(BM_JacobiParallelNaive)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK(BM_JacobiParallel)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK_MAIN();
