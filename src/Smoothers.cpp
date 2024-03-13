#include "Smoothers.hpp"


void residual(Lattice &mesh, const std::vector<double> &u, const std::vector<double> &b, std::vector<double> &r) {
	// @TODO: add also parallelization here
	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		r[i] = 4.0 * u[i] - u[nord] - u[sud] - u[ovest] - u[est] - b[i];
	}
}


void residual(Graph &graph, const std::vector<double> &u, const std::vector<double> &b, std::vector<double> &r) {  
	// @TODO: add also parallelization here
	std::pair<Index, double> neighbour;
	for (int i = 0; i < graph.get_nodes().size(); i++) {

		if(graph.get_bool_boundary().at(i) == false){

			r[i] = 4.0 * u[i] - b[i];
			
			for(int j = 0; j < graph.num_neighbours(i); j++){

					neighbour = graph.get_neighbours().at(graph.get_nodes().at(i).index_node + j);

					if(neighbour.second == 0.125){
						
						r[i] = r[i] - u[neighbour.first];
					
					}
			}
		}
	}  

}


void gseidel(Lattice &mesh, std::vector<double> &u, const std::vector<double> &b) {
	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + u[nord] + u[sud] + u[ovest] + u[est]);
	}
}


void gseidel(Graph &graph, std::vector<double> &u, const std::vector<double> &b) {
	for (int i = 0; i < graph.get_nodes().size(); i++) {

		if(graph.get_bool_boundary().at(i) == false){

			u[i] = 0.25 * b[i];

			for(int j : graph.get_cardinal_neighbours(i)){
				
				u[i] = u[i] + 0.25 * u[j];
			
			}
		}
	}
}


void jacobi(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}
}


void jacobi_parallel_naive(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	// implementazione naive, forse openmp non riesce a parallelizzare i for di c++11
	#pragma omp parallel for

	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}

	#pragma omp barrier
}


void jacobi_parallel(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	// adattato il for loop alla sua forma tradizionale
	const auto inner_nodes = mesh.get_inner_nodes();

	#pragma omp parallel for

	for (int j = 0; j < inner_nodes.size() ; j ++) {
		auto i = inner_nodes[j];

		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}

	#pragma omp barrier
}
