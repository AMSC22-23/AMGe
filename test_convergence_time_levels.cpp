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


static void DUE(benchmark::State& state) {
    int n = state.range(0);
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, n);


    Multigrid due(mesh, PRE_STEP, POST_STEP, 2);


    std::vector<double> b(mesh.numel());
    std::vector<double> u(mesh.numel());
    std::vector<double> r(mesh.numel());


    mesh.evaluate_forcing_term(b, f);
    mesh.evaluate_boundary_conditions(u, g);


    do{
        due.step(u, b);
        residual(mesh, u, b, r);
    }while(norm(r) > 1.e-16);
}


static void TRE(benchmark::State& state) {
    int n = state.range(0);
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, n);


    Multigrid tre(mesh, PRE_STEP, POST_STEP, 2);


    std::vector<double> b(mesh.numel());
    std::vector<double> u(mesh.numel());
    std::vector<double> r(mesh.numel());


    mesh.evaluate_forcing_term(b, f);
    mesh.evaluate_boundary_conditions(u, g);


    do{
        tre.step(u, b);
        residual(mesh, u, b, r);
    }while(norm(r) > 1.e-16);
}


static void QUATTRO(benchmark::State& state) {
    int n = state.range(0);
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, n);


    Multigrid quattro(mesh, PRE_STEP, POST_STEP, 2);


    std::vector<double> b(mesh.numel());
    std::vector<double> u(mesh.numel());
    std::vector<double> r(mesh.numel());


    mesh.evaluate_forcing_term(b, f);
    mesh.evaluate_boundary_conditions(u, g);


    do{
        quattro.step(u, b);
        residual(mesh, u, b, r);
    }while(norm(r) > 1.e-16);
}


BENCHMARK(DUE)->Arg(33)->Arg(65)->Arg(129)->Arg(257)->Arg(513)->Arg(1025);
BENCHMARK(TRE)->Arg(33)->Arg(65)->Arg(129)->Arg(257)->Arg(513)->Arg(1025);
BENCHMARK(QUATTRO)->Arg(33)->Arg(65)->Arg(129)->Arg(257)->Arg(513)->Arg(1025);
BENCHMARK_MAIN();