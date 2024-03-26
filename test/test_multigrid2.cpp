#include <iostream>
#include <cmath>
#include <vector>
#include "../include/Multigrid2.hpp"


#ifndef PRE_SMOOTHING_STEPS
#define PRE_SMOOTHING_STEPS 10
#endif

#ifndef POST_SMOOTHING_STEPS
#define POST_SMOOTHING_STEPS 10
#endif


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


int main (int argc, char *argv[]) {
	const int N = 129;
	
	
	Graph fine = Graph();
	Multigrid2 two(fine, PRE_SMOOTHING_STEPS, POST_SMOOTHING_STEPS, 2);
	Multigrid2 three(fine, PRE_SMOOTHING_STEPS, POST_SMOOTHING_STEPS, 3);
	Multigrid2 four(fine, PRE_SMOOTHING_STEPS, POST_SMOOTHING_STEPS, 4);


	std::vector<double> u_fine(fine.numel());
	std::vector<double> b(fine.numel());
	std::vector<double> r_fine(fine.numel());


	fine.evaluate_boundary_conditions(u_fine, g);
	fine.evaluate_forcing_term(b, f);
	

	for (int i = 0; i < 300; ++i) {
		#ifdef MULTIGRID2
			two.step(u_fine, b);
		#elif MULTIGRID3
			three.step(u_fine, b);
		#elif MULTIGRID4
			four.step(u_fine, b);
		#else
			for (int steps = 0; steps < (PRE_SMOOTHING_STEPS + POST_SMOOTHING_STEPS); ++steps) {
				gseidel(fine, u_fine, b);
			}
		#endif

		residual(fine, u_fine, b, r_fine);
		std::cout << i<< " " << norm(r_fine) << std::endl;
	}

	return 0;
}
