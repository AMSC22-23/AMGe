#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include <iostream>
#include <vector>
#include "Lattice.hpp"
#include "Smoothers.hpp"


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
		if (levels < 1) {
			std::cerr << "Multigrid method has to work with at least one level" << std::endl;
			exit(EXIT_FAILURE);
		}

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


	void step(std::vector<double> &u, const std::vector<double> &b, int lvl = 0) {
		Lattice &fine = meshes[lvl];

		if (lvl == levels) {
			// solve well the problem at the coarsest grid, for now we do a fixed amount of iterations
			for (int it = 0; it < 300; ++it) {
				gseidel(fine, u, b);
			}
		}
		else {
			for (int pre = 0; pre < pre_smoothing_steps; ++pre) {
				gseidel(fine, u, b);
			}


			residual(fine, u, b, res[lvl]);
			Lattice &coarse = meshes[lvl+1];
			// projection
			// recursive call
			// prolongation
			// correction



			for (int post = 0; post < post_smoothing_steps; ++post) {
				gseidel(meshes[lvl], u, b);
			}
		}
	}


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
