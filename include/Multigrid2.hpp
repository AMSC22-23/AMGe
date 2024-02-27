#ifndef __MULTIGRID2_HPP__
#define __MULTIGRID2_HPP__


#include <iostream>
#include <vector>
#include "Graph.hpp"
#include "Smoothers.hpp"


class Multigrid2 {
public:
	Multigrid2(
		//@note: pass by const reference to avoid copy
		  Graph mesh
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

		// i vettori di errore, residuo
		for (int lvl = 0; lvl < levels; ++lvl) {
			err.push_back(std::vector<double>(meshes[lvl].numel()));
			res.push_back(std::vector<double>(meshes[lvl].numel()));

			if (lvl == 0) {
				u_internal.push_back(std::vector<double>());
				b_internal.push_back(std::vector<double>());
			}
			else {
				u_internal.push_back(std::vector<double>(meshes[lvl].numel()));
				b_internal.push_back(std::vector<double>(meshes[lvl].numel()));
			}
		}
	}

	//@note: nothing wrong with this, but is useless
	~Multigrid2() {

	}


	void step(std::vector<double> &u, const std::vector<double> &b, int lvl = 0) {
		Graph &fine = meshes[lvl];

		if (lvl == (levels-1)) {
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
			Graph &coarse = meshes[lvl+1];

			// projection
			fine.project_on_coarse(coarse, res[lvl], b_internal[lvl+1]);

			// recursive call
			coarse.evaluate_zero(u_internal[lvl+1]);
			step(u_internal[lvl+1], b_internal[lvl+1], lvl+1);

			// prolongation
			fine.interpolate_on_fine(coarse, err[lvl], u_internal[lvl+1]);

			// correction
			for (Index i : fine.get_inner_nodes()) {
				u[i] -= err[lvl][i];
			}

			for (int post = 0; post < post_smoothing_steps; ++post) {
				gseidel(fine, u, b);
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

	std::vector<std::vector<double>> u_internal;
	std::vector<std::vector<double>> b_internal;
	std::vector<std::vector<double>> err;
	std::vector<std::vector<double>> res;
	std::vector<Graph> meshes;
};

#endif
