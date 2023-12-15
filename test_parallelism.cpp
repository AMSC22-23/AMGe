#include <iostream>
#include <cmath>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"
#include "Smoothers.hpp"


using namespace std;


#ifndef PRE_SMOOTHING_STEPS
#define PRE_SMOOTHING_STEPS 10
#endif

#ifndef POST_SMOOTHING_STEPS
#define POST_SMOOTHING_STEPS 10
#endif


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
	const int N = 257;


	Lattice fine(0.0, 0.0, 1.0, 1.0, N);


	std::vector<double> u_old(fine.numel());
	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());
	std::vector<double> r(fine.numel());


	set_initial_guess(fine, u_old, g);
	set_initial_guess(fine, u, g);
	fine.evaluate_forcing_term(b, f);


	time_t start = time(nullptr);

	for (int i = 0; i < 10000; ++i) {
		for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
			jacobi(fine, u, u_old, b);
			std::swap(u, u_old);
		}
		residual(fine, u, b, r);
	}


	time_t end = time(nullptr);
	double time_elapsed = difftime(end, start);
	std::cout << "Il programma ha impiegato " << time_elapsed << " secondi.\n";





	set_initial_guess(fine, u_old, g);
	set_initial_guess(fine, u, g);
	fine.evaluate_forcing_term(b, f);


	time_t start2 = time(nullptr);

	for (int i = 0; i < 10000; ++i) {
		for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
			jacobi_parallel(fine, u, u_old, b);
			std::swap(u, u_old);
		}
		residual(fine, u, b, r);
	}


	time_t end2 = time(nullptr);
	time_elapsed = difftime(end2, start2);
	std::cout << "Il programma parallelo ha impiegato " << time_elapsed << " secondi.\n";



	return 0;
}
