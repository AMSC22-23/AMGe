#include <benchmark/benchmark.h>
#include <iostream>
#include <vector>
#include "Multigrid.hpp"


#define PRE_STEP 10
#define POST_STEP 10


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


static void BM_TwoLevel(benchmark::State& state) {
    int n = state.range(0);
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, n+1);


    Multigrid solver(mesh, PRE_STEP, POST_STEP, 2);


    std::vector<double> u(mesh.numel());
    std::vector<double> b(mesh.numel());
    std::vector<double> r(mesh.numel());


	for (auto _ : state) {
		int it = 0;
		mesh.evaluate_boundary_conditions(u, g);
		mesh.evaluate_forcing_term(b, f);


		do {
			solver.step(u, b);
			residual(mesh, u, b, r);
			++it;
		} while(norm(r) > 1.e-9 and it < 10000);
	}
}


static void BM_FiveLevel(benchmark::State& state) {
    int n = state.range(0);
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, n+1);


    Multigrid solver(mesh, PRE_STEP, POST_STEP, 5);


    std::vector<double> u(mesh.numel());
    std::vector<double> b(mesh.numel());
    std::vector<double> r(mesh.numel());


	for (auto _ : state) {
		int it = 0;
		mesh.evaluate_boundary_conditions(u, g);
		mesh.evaluate_forcing_term(b, f);


		do {
			solver.step(u, b);
			residual(mesh, u, b, r);
			++it;
		} while(norm(r) > 1.e-9 and it < 10000);
	}
}


static void BM_SevenLevel(benchmark::State& state) {
    int n = state.range(0);
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, n+1);


    Multigrid solver(mesh, PRE_STEP, POST_STEP, 7);


    std::vector<double> u(mesh.numel());
    std::vector<double> b(mesh.numel());
    std::vector<double> r(mesh.numel());


	for (auto _ : state) {
		int it = 0;
		mesh.evaluate_boundary_conditions(u, g);
		mesh.evaluate_forcing_term(b, f);


		do {
			solver.step(u, b);
			residual(mesh, u, b, r);
			++it;
		} while(norm(r) > 1.e-9 and it < 10000);
	}
}


BENCHMARK(BM_TwoLevel)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK(BM_FiveLevel)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK(BM_SevenLevel)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024);
BENCHMARK_MAIN();
