#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include <array>


template <unsigned int pre_smoothing_steps, unsigned int post_smoothing_steps, unsigned int levels>
class Multigrid {
public:
	Multigrid() {
		// qui ci vanno un po' di allocazioni
	}

	~Multigrid() {

	}

	void step();
	void solve(double tol = 1e-9, int max_iter = 10000);
	std::vector<double> get_solution();

private:
	std::array<std::vector<double>, levels> u;
	std::array<std::vector<double>, levels> err;
	std::array<std::vector<double>, levels-1> res;
};

#endif
