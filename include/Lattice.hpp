#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__


#include <iostream>
#include <vector>
#include <tuple>
#include "Utils.hpp"


using Index = int;


class Lattice {
public:
	Lattice(double x_corner_, double y_corner_, double width_, double height_, unsigned int N_);

	/* 
	 *	number of total nodes (boundary and internal of the mesh) used for allocating the solution vector like:
	 *		std::vector<double> u(mesh.numel());
	*/
	int numel();

	
	const std::vector<Index>& get_inner_nodes();
	std::tuple<int, int, int, int> get_cardinal_neighbours(Index i);
	std::tuple<int, int, int, int> get_diagonal_neighbours(Index i);
	Index index(int i, int j);
	std::pair<int, int> inverse_index(Index i);
	Lattice build_coarse();
	void project_on_coarse(Lattice &coarse, const std::vector<double> &u, std::vector<double> &v);
	void interpolate_on_fine(Lattice &coarse, std::vector<double> &u_fine, const std::vector<double> &u_coarse);
	void evaluate_zero(std::vector<double> &u);
	void evaluate_forcing_term(std::vector<double> &b, double (*f)(double x, double y));
	void evaluate_boundary_conditions(std::vector<double> &u, double (*f)(double x, double y));
	void evaluate_function(std::vector<double> &u, double (*f)(double x, double y));
	void print_vector(const std::vector<double> &u, const char *name);


	// test functions
	void test_constructor();
	void test_inner_nodes();
	void test_index();
	void test_inverse_index();


private:
	double  x_corner,
			y_corner,
			width,
			height;

	double h;

	int N;
	std::vector<Index> inner_nodes;
	std::vector<Index> boundary_nodes;

	bool minimal;
};


#endif
