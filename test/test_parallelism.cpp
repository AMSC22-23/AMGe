#include <iostream>
#include <cmath>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"
#include "Smoothers.hpp"
#include <chrono>


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
	

	std::chrono::high_resolution_clock::time_point t0;
	std::chrono::high_resolution_clock::time_point t1;
	std::int64_t medium_time_elapsed;
	std::int64_t time_elapsed;



	std::vector<int> mesh_size;

	mesh_size.push_back(33);
	mesh_size.push_back(65);
	mesh_size.push_back(129);
	mesh_size.push_back(257);




	for(int d : mesh_size){
			Lattice fine(0.0, 0.0, 1.0, 1.0, d);


			std::vector<double> u_old(fine.numel());
			std::vector<double> u(fine.numel());
			std::vector<double> b(fine.numel());
			std::vector<double> r(fine.numel());


			for(int n_iterations = 100 ; n_iterations < 10001; n_iterations *= 10){
					
					set_initial_guess(fine, u_old, g);
					set_initial_guess(fine, u, g);
					fine.evaluate_forcing_term(b, f);


					medium_time_elapsed = 0;


					for (int i = 0; i < n_iterations; ++i) {
						t0 = std::chrono::high_resolution_clock::now();
					
						for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
								jacobi(fine, u, u_old, b);
								std::swap(u, u_old);
						}
						
						t1 = std::chrono::high_resolution_clock::now();
						time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
						medium_time_elapsed += time_elapsed;
					}


					medium_time_elapsed = medium_time_elapsed / n_iterations;
					std::cout << "Jacobi seriale con mesh di size: " << d << " in " << n_iterations << " iterazioni termina in " << medium_time_elapsed << " microsecondi.\n";




					set_initial_guess(fine, u_old, g);
					set_initial_guess(fine, u, g);
					fine.evaluate_forcing_term(b, f);

					medium_time_elapsed = 0;


					for (int i = 0; i < n_iterations; ++i) {
						t0 = std::chrono::high_resolution_clock::now();
						
						for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
								jacobi_parallel(fine, u, u_old, b);
								std::swap(u, u_old);
						}

						t1 = std::chrono::high_resolution_clock::now();
						time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
						medium_time_elapsed += time_elapsed;
					}


					medium_time_elapsed = medium_time_elapsed / n_iterations;
					std::cout << "Jacobi parallelo con mesh di size: " << d << " in " << n_iterations << " iterazioni termina in " << medium_time_elapsed << " microsecondi.\n";
					std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
			}


			std::cout << " " << std::endl;
			std::cout << " " << std::endl;
	}


	std::cout << " " << std::endl;
	std::cout << " " << std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "--------------------------------------------------------REVERSE-------------------------------------------------" << std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << std::endl;


	for(int d : mesh_size){
			Lattice fine(0.0, 0.0, 1.0, 1.0, d);


			std::vector<double> u_old(fine.numel());
			std::vector<double> u(fine.numel());
			std::vector<double> b(fine.numel());
			std::vector<double> r(fine.numel());


			for(int n_iterations = 100 ; n_iterations < 10001; n_iterations *= 10){
					
					set_initial_guess(fine, u_old, g);
					set_initial_guess(fine, u, g);
					fine.evaluate_forcing_term(b, f);


					medium_time_elapsed = 0;


					for (int i = 0; i < n_iterations; ++i) {
						t0 = std::chrono::high_resolution_clock::now();
					
						for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
								jacobi_parallel(fine, u, u_old, b);
								std::swap(u, u_old);
						}
						
						t1 = std::chrono::high_resolution_clock::now();
						time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
						medium_time_elapsed += time_elapsed;
					}


					medium_time_elapsed = medium_time_elapsed / n_iterations;
					std::cout << "Jacobi parallelo con mesh di size: " << d << " in " << n_iterations << " iterazioni termina in " << medium_time_elapsed << " microsecondi.\n";




					set_initial_guess(fine, u_old, g);
					set_initial_guess(fine, u, g);
					fine.evaluate_forcing_term(b, f);

					medium_time_elapsed = 0;


					for (int i = 0; i < n_iterations; ++i) {
						t0 = std::chrono::high_resolution_clock::now();
						
						for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
								jacobi(fine, u, u_old, b);
								std::swap(u, u_old);
						}

						t1 = std::chrono::high_resolution_clock::now();
						time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
						medium_time_elapsed += time_elapsed;
					}


					medium_time_elapsed = medium_time_elapsed / n_iterations;
					std::cout << "Jacobi seriale con mesh di size: " << d << " in " << n_iterations << " iterazioni termina in " << medium_time_elapsed << " microsecondi.\n";
					std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
			}


			std::cout << " " << std::endl;
			std::cout << " " << std::endl;
	}


	return 0;
}
