#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include <vector>


class Multigrid {
public:
	Multigrid(
		  unsigned int pre_smoothing_steps_
		, unsigned int post_smoothing_steps_
		, unsigned int levels_
	) :
		  pre_smoothing_steps(pre_smoothing_steps_)
		, post_smoothing_steps(post_smoothing_steps_)
		, levels(levels_)
	{
		// qui ci vanno un po' di allocazioni
	}

	~Multigrid() {

	}

	void step();
	void solve(double tol = 1e-9, int max_iter = 10000);
	std::vector<double> get_solution();

private:
	const unsigned int pre_smoothing_steps;
	const unsigned int post_smoothing_steps;
	const unsigned int levels;

	std::vector<std::vector<double>> u;
	std::vector<std::vector<double>> err;
	std::vector<std::vector<double>> res;
};

#endif
