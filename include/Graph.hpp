#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__


#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include "Utils.hpp"
#include "Triplet.hpp"


using Index = int;


class Graph {
public:
	Graph(double x_corner_, double y_corner_, double width_, double height_, unsigned int N_);

	/* 
	 *	number of total nodes (boundary and internal of the mesh) used for allocating the solution vector like:
	 *		std::vector<double> u(mesh.numel());
	*/
	int numel();

	


	// return a vector where each element rapresent a node and points to its first neighbour in vector neighbours
	const std::vector<Triplet>& get_nodes();


	// return a vector of nodes identified by their index
	const std::vector<std::pair<Index, double>>& get_neighbours();


	//return index of nodes in the boundary
	const std::vector<Index>& get_boundary();


	// return a vector of the same size of nodes, it specifies if the element is in the boundary or not
	const std::vector<bool>& get_bool_boundary();


	// return a vector of the same size of neighbours, it specifies if the neighbour is cardinal or  not
	const std::vector<bool>& get_bool_cardinal_neighbour();


	// return a vector correspondence between an index of the node and its coordinates in the system (coordinates are the key)
	const std::map<std::pair<int,int>, Index>& get_position_coordinates();


	// The solution vector is a 1D vector, this map and its inverse helps converting 2D indeces (natural for a 2D mesh) into 1D equivalents
	Index index(int i, int j);
	std::pair<int, int> inverse_index(Index i);


	// This functions are safe when called on indeces in inner_nodes, the special cases won't be implemented
	std::vector<int> get_cardinal_neighbours(Index i);
	std::vector<int> get_diagonal_neighbours(Index i);
	std::pair<Index, int> get_node_neighbours(Index i);


	Graph build_coarse();


	//return number of neighbours
    int num_neighbours(Index i);

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
    std::vector<Triplet> nodes;
	std::vector<std::pair<Index, double>> neighbours;
    std::vector<bool> boundary_bool;
	std::vector<bool> cardinal_neighbour_bool;
	std::map<std::pair<int, int>, Index> position_coordinates;
	std::vector<Index> boundary;


	bool minimal;
};


#endif