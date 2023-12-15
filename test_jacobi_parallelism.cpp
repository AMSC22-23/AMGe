#include <iostream>
#include <cmath>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"
#include "Smoothers.hpp"
#include <chrono>


using namespace std;


void set_initial_guess(Lattice &mesh, std::vector<double> &u, double (*g)(double x, double y)) {
	mesh.evaluate_function(u, g);

	for (Index i : mesh.get_inner_nodes()) {
		u[i] = 0.0;
	}
}


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


#define ITERATIONS_PER_MEASURE 10


int main (int argc, char *argv[]) {
	/*
		* @TODO(big): measure the speedup of openmp parallelization
		*	. use std::chrono::high_resolution_clock as we did in amsc laboratories
		*	. measure let's say 100 iterations of the method and average the time (both serial and parallel versions)
		*	. do a scalability test: for meshes with increasingly number of points (from 3 to 1000 for example) I want to
		*				know the execution times of one iteration (of course averaged across 100 for example) of serial and parallel jacobi
		*
		*	The initial values of u are not very important, since an iteration is computationally the same for every data
		*	At the end print to std::cout for every experiment
		*		number of points in the mesh, jacobi serial time, jacobi parallel time
		*
		* 
	*/
	

	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;


	std::vector<int> N = {32, 64, 128, 256, 512, 1024, 2048};


	for (int n : N) {
		Lattice fine(0.0, 0.0, 1.0, 1.0, n+1);


		std::vector<double> old(fine.numel());
		std::vector<double> u(fine.numel());
		std::vector<double> b(fine.numel());


		fine.evaluate_boundary_conditions(u, g);
		fine.evaluate_forcing_term(b, f);


		start = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < ITERATIONS_PER_MEASURE; ++i) {
				jacobi_parallel(fine, u, old, b);
			}
		end = std::chrono::high_resolution_clock::now();
		auto parallel_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / ITERATIONS_PER_MEASURE;


		start = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < ITERATIONS_PER_MEASURE; ++i) {
				jacobi(fine, u, old, b);
			}
		end = std::chrono::high_resolution_clock::now();
		auto serial_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / ITERATIONS_PER_MEASURE;


		std::cout << n+1 << "\t" << serial_time << "\t" << parallel_time << std::endl;
	}


	return 0;
}
