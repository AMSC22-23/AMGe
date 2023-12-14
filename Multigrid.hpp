#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include <iostream>
#include <vector>
#include "Lattice.hpp"


class Multigrid {
public:
	Multigrid(
		  Lattice mesh
		, unsigned int pre_smoothing_steps_
		, unsigned int post_smoothing_steps_
		, unsigned int levels_
	) :
		  pre_smoothing_steps(pre_smoothing_steps_)
		, post_smoothing_steps(post_smoothing_steps_)
		, levels(levels_)
	{
		// qui ci vanno un po' di allocazioni, inizio dalle mesh
		for (int lvl = 0; lvl < levels; ++lvl) {
			if (lvl == 0) {
				meshes.push_back(mesh);
			}
			else {
				meshes.push_back(meshes[lvl-1].build_coarse());
			}
		}
	}

	~Multigrid() {

	}


	void step();
	void solve(double tol = 1e-9, int max_iter = 10000);
	std::vector<double> get_solution();


	void test_allocation() {
		std::cout << "Testing the mesh allocation in multigrid class" << std::endl;
		std::cout << "-----------------------------------------------------" << std::endl;
		for (int lvl = 0; lvl < levels; ++lvl) {
			std::cout << "level " << lvl << std::endl;
			meshes[lvl].test_constructor();
			std::cout << "========================" << std::endl;
		}
	}

private:
	const unsigned int pre_smoothing_steps;
	const unsigned int post_smoothing_steps;
	const unsigned int levels;

	std::vector<std::vector<double>> u;
	std::vector<std::vector<double>> err;
	std::vector<std::vector<double>> res;
	std::vector<Lattice> meshes;
};

#endif
