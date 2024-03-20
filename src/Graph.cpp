#include "Graph.hpp"
#include "Utils.hpp"


Graph::Graph(double x_corner_, double y_corner_, double width_, double height_, unsigned int N_) {
	x_corner = x_corner_;
	y_corner = y_corner_;
	width    = width_;
	height   = height_;


	if (width != height) {
		//@note: I know that `std::cout` is clunky but is more aligned with C++
		//       Do not worry, C++23 has std::format!
		fprintf(stderr, "Mesh has to be square\n");
		exit(EXIT_FAILURE);
	}


	N = N_;
	h = width / (N-1);

    int node = 0;
    int neighbour = 0;
	index_bool=0;
	int node_neighbours;


	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == 0 || i == N-1 || j == 0 || j == N-1) {
                boundary_bool.push_back(true);
				boundary.push_back(node);
				nodes.push_back(Node(neighbour, i, j));
				position_coordinates[std::make_pair(i,j)] = node;
			} else {
                boundary_bool.push_back(false);
				//boundary_bool.insert(boundary_bool.begin(), false);
				//nodes.insert(nodes.begin(), Triplet(neighbour, i, j));
				nodes.push_back(Node(neighbour, i, j));
				position_coordinates[std::make_pair(i,j)] = node;
				index_bool++;
			}
            
			node_neighbours = neighbour;


            if(i < N - 1){
				neighbours.push_back(std::make_pair(node + N, 0.125));
				cardinal_neighbour_bool.push_back(true);
                neighbour++;
            }


            if(i > 0){
				neighbours.push_back(std::make_pair(node - N, 0.125));
				cardinal_neighbour_bool.push_back(true);
                neighbour++;
            }


            if(j > 0){
				neighbours.push_back(std::make_pair(node - 1, 0.125));
				cardinal_neighbour_bool.push_back(true);
                neighbour++;
            }

            if(i < N - 1){
				neighbours.push_back(std::make_pair(node + 1, 0.125));
				cardinal_neighbour_bool.push_back(true);
                neighbour++;
            }

            if(i < N - 1 && j < N - 1){
				neighbours.push_back(std::make_pair(node + N + 1, 0.0625));
				cardinal_neighbour_bool.push_back(false);
                neighbour++;
            }

            if(i < N - 1 && j > 0){
				neighbours.push_back(std::make_pair(node + N - 1, 0.0625));
				cardinal_neighbour_bool.push_back(false);
                neighbour++;
            }

            if(i > 0 && j > 0){
				neighbours.push_back(std::make_pair(node - N - 1, 0.0625));
				cardinal_neighbour_bool.push_back(false);
                neighbour++;
            }

            if(i > 0 && j < N - 1){
				neighbours.push_back(std::make_pair(node - N + 1, 0.0625));
				cardinal_neighbour_bool.push_back(false);
                neighbour++;
            }

			nodes.at(node).node_neighbours = neighbour - node_neighbours;

            node++;

		}
	}


	minimal = !is_2nplusone(N);
}


int Graph::numel() {
	return N*N;
}


const std::vector<Node>& Graph::get_nodes(){
	return nodes;
}


const std::vector<std::pair<Index, double>>& Graph::get_neighbours(){
	return neighbours;
}

const std::vector<Index>& Graph::get_boundary(){
	return boundary;
}


const std::vector<bool>& Graph::get_bool_boundary(){
	return boundary_bool;
}


const std::vector<bool>& Graph::get_bool_cardinal_neighbour(){
	return cardinal_neighbour_bool;
}




const std::map<std::pair<int, int>, Index>& Graph::get_position_coordinates(){
	return position_coordinates;
}


std::vector<int> Graph::get_cardinal_neighbours(Index i) {
	std::vector<int> cardinal_neighbours;
	
	for(int j = 0; j < num_neighbours(i); j++){
		if(get_bool_cardinal_neighbour().at(nodes.at(i).index_node + j) == true){
			
			cardinal_neighbours.push_back(neighbours.at(nodes.at(i).index_node + j).first);
		}
	}
	
	return cardinal_neighbours;

}


std::pair<Index, int> Graph::get_node_neighbours(Index i){
	return std::make_pair(nodes.at(i).index_node, num_neighbours(i));
}


std::vector<int> Graph::get_diagonal_neighbours(Index i) {
	std::vector<int> cardinal_neighbours;
	
	for(int j = 0; j < num_neighbours(i); j++){
		if(get_bool_cardinal_neighbour().at(nodes.at(i).index_node + j) == false){
			
			cardinal_neighbours.push_back(neighbours.at(nodes.at(i).index_node + j).first);
		}
	}
	
	return cardinal_neighbours;
}


Index Graph::index(int i, int j) {
	return position_coordinates.at(std::make_pair(i,j));
}


std::pair<int, int> Graph::inverse_index(Index i) {
	return  std::make_pair(nodes.at(i).x, nodes.at(i).y);
}



Graph Graph::build_coarse() {
	if (minimal) {
		fprintf(stderr, "Mesh can't create coarse\n");
		exit(EXIT_FAILURE);
	}


	return Graph(
		  x_corner
		, y_corner
		, width
		, height
		, 1 + ((N - 1) / 2)
	);
}

int Graph::num_neighbours(Index i){
	return nodes.at(i).node_neighbours;
}


void Graph::project_on_coarse(Graph &coarse, const std::vector<double> &u, std::vector<double> &v) {
	std::pair<Index, double> neighbour;

	for (auto &x : v) {
		x = 0.0;
	}

	for (int i = 2; i < (N-1); i = i + 2) {
		for (int j = 2; j < (N-1); j = j + 2) {
			int ii = i / 2;
			int jj = j / 2;


			v[coarse.index(ii, jj)] = 0.25   * (u[index(i,j)]);

            for(int z = 0; z < num_neighbours(index(i,j)); z++){
				
				neighbour = neighbours.at(nodes.at(index(i,j)+z).index_node);
                v[coarse.index(ii, jj)] = u[neighbour.first] * neighbour.second;
            }
		}
	}
}


void Graph::interpolate_on_fine(Graph &coarse, std::vector<double> &u_fine, const std::vector<double> &u_coarse) {
	// devo iterare solo sui nodi interni
	for (int i = 1; i < N-1; ++i) {
		for (int j = 1; j < N-1; ++j) {
			if (i % 2 == 0 and j % 2 == 0) {
				u_fine[index(i,j)] = u_coarse[coarse.index(i/2, j/2)];
			}
			else if (i % 2 == 1 and j % 2 == 0) {
				u_fine[index(i,j)] =
					0.5 * (
							u_coarse[coarse.index(    i/2, j/2)]
							+ u_coarse[coarse.index(1 + i/2, j/2)]
						  );
			}
			else if(i % 2 == 0 and j % 2 == 1) {
				u_fine[index(i,j)] =
					0.5 * (
							u_coarse[coarse.index(i/2,     j/2)]
							+ u_coarse[coarse.index(i/2, 1 + j/2)]
						  );
			}
			else {
				u_fine[index(i,j)] =
					0.25 * (
							u_coarse[coarse.index(    i/2,     j/2)]
							+ u_coarse[coarse.index(    i/2, 1 + j/2)]
							+ u_coarse[coarse.index(1 + i/2,     j/2)]
							+ u_coarse[coarse.index(1 + i/2, 1 + j/2)]
						   );
			}
		}
	}
}


void Graph::evaluate_zero(std::vector<double> &u) {
	for (int i = 0; i < nodes.size(); ++i) {
		u[i] = 0.0;
	}
}


void Graph::evaluate_forcing_term(std::vector<double> &b, double (*f)(double x, double y)) {
	std::pair<int,int> value;

	for(Index z = 0; z < nodes.size(); z++){
		value = inverse_index(z);

		double x = x_corner + value.first * h * width;
		double y = y_corner + value.second * h * height;


		b[index(value.first,value.second)] = h * h * f(x,y);
	}
}


void Graph::evaluate_boundary_conditions(std::vector<double> &u, double (*f)(double x, double y)) {
	for (Index i : boundary){
			const auto [j, k] = inverse_index(i);
			const double x = x_corner + j * h * width;
			const double y = y_corner + k * h * height;
			u[i] = f(x,y);
	}
}


void Graph::evaluate_function(std::vector<double> &u, double (*f)(double x, double y)) {
	std::pair<int,int> value;

	for(Index z = 0; z < nodes.size(); z++){
			value = inverse_index(z);
			
			double x = x_corner + value.first * h * width;
			double y = y_corner + value.second * h * height;


			u[index(value.first,value.second)] = f(x,y);

	}
}


void Graph::print_vector(const std::vector<double> &u, const char *name) {
	std::cout << name << std::endl;

	for (int j = 0; j < N; ++j) {
		for (int i = 0; i < N; ++i) {
			printf("%e\t", u[index(i, j)]);
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}


void Graph::test_constructor() {
	std::cout << "Lattice mesh" << std::endl;

	std::cout
		<< " bottom left corner position "
		<< "("
		<< x_corner
		<< ", "
		<< y_corner
		<< ")"
		<< std::endl;

	std::cout
		<< " number of subdivisions "
		<< N
		<< std::endl;
}


void Graph::test_inner_nodes() {
	for (int i = 0; i < nodes.size(); i++) {
		if(boundary_bool[i] == false){
			std::cout << i << std::endl;
		}
	}
}


void Graph::test_index() {
	for (int j = 0; j < N; ++j) {
		for (int i = 0; i < N; ++i) {
			std::cout << index(i,j) << std::endl;
		}
	}
}


void Graph::test_inverse_index() {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			const auto [k, l] = inverse_index(index(i,j));

			if (k != i || l != j) {
				std::cout << "errore" << std::endl;
			}
		}
	}
}