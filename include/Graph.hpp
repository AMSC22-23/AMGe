#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__


#include <iostream>
#include <vector>
#include <tuple>
#include "Utils.hpp"


using Index = int;


class Graph {
public:
	Graph(double x_corner_, double y_corner_, double width_, double height_, unsigned int N_);

	/* 
	 *	number of total nodes (boundary and internal of the mesh) used for allocating the solution vector like:
	 *		std::vector<double> u(mesh.numel());
	*/
	int numel();

	
	// There was the intention to return an iterator to the inner nodes of the mesh, but was inconvenient for the parallelization
	const std::vector<Index>& get_inner_nodes();


	// The solution vector is a 1D vector, this map and its inverse helps converting 2D indeces (natural for a 2D mesh) into 1D equivalents
	Index index(int i, int j);
	std::pair<int, int> inverse_index(Index i);


	// This functions are safe when called on indeces in inner_nodes, the special cases won't be implemented
	std::tuple<int, int, int, int> get_cardinal_neighbours(Index i);
	std::tuple<int, int, int, int> get_diagonal_neighbours(Index i);


	Graph build_coarse();

    int num_neighbours(Index i){};

	// no check whatsoever of mesh size, assume that coarse == this.build_coarse()
	void project_on_coarse(Graph &coarse, const std::vector<double> &u, std::vector<double> &v);
	void interpolate_on_fine(Graph &coarse, std::vector<double> &u_fine, const std::vector<double> &u_coarse);


	void evaluate_zero               (std::vector<double> &u);
	void evaluate_forcing_term       (std::vector<double> &b, double (*f)(double x, double y));
	void evaluate_boundary_conditions(std::vector<double> &u, double (*g)(double x, double y));
	void evaluate_function           (std::vector<double> &u, double (*f)(double x, double y));


	void print_vector(const std::vector<double> &u, const char *name);


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
    std::vector<Index> nodes;
    std::vector<Index> neighbours;
    std::vector<double> weights;
    std::vector<bool> boundary_bool;

	bool minimal;
};


#endif
