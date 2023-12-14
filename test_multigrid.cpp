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

void two_level(Lattice &fine, Lattice &coarse, std::vector<double> &u, std::vector<double> &b){
	std::vector<double> r_fine(fine.numel());
	std::vector<double> e_fine(fine.numel());
	std::vector<double> r_coarse(coarse.numel());
	std::vector<double> e_coarse(coarse.numel());

	for (int pre = 0; pre < PRE_SMOOTHING_STEPS; ++pre) {
			gseidel(fine, u, b);
			//jacobi(fine, u, b);
		}
	
		residual(fine, u, b, r_fine);

		fine.project_on_coarse(coarse, r_fine, r_coarse);


		for (auto &x : e_coarse) {
			x = 0.0;
		}

		for (int coarse_it = 0; coarse_it < 300; ++coarse_it) {
			gseidel(coarse, e_coarse, r_coarse);
			//jacobi(coarse, e_coarse, r_coarse);
		}

		fine.interpolate_on_fine(coarse, e_fine, e_coarse);

		for (Index i : fine.get_inner_nodes()) {
			u[i] -= e_fine[i];
		}

		for (int post = 0; post < POST_SMOOTHING_STEPS; ++post) {
			gseidel(fine, u, b);
			//jacobi(fine, u, b);
		}
}


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


int main (int argc, char *argv[]) {
	const int N = 129;
	
	
	Lattice fine(0.0, 0.0, 1.0, 1.0, N);
	Lattice coarse = fine.build_coarse();


	std::vector<double> u_fine(fine.numel());
	std::vector<double> b(fine.numel());
	std::vector<double> r_fine(fine.numel());
	std::vector<double> e_fine(fine.numel());


	std::vector<double> r_coarse(coarse.numel());
	std::vector<double> e_coarse(coarse.numel());


	set_initial_guess(fine, u_fine, g);
	fine.evaluate_forcing_term(b, f);
	



	for (int i = 0; i < 1000; ++i) {
		#ifdef MULTIGRID
			two_level(fine, coarse, u_fine, b);
		#else
			for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
				gseidel(fine, u_fine, b);
				//jacobi(fine,u_fine,b);
			}
		#endif

		residual(fine, u_fine, b, r_fine);
		std::cout << norm(r_fine) << std::endl;
	}


	return 0;
}
