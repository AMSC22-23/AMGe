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


	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == 0 || i == N-1 || j == 0 || j == N-1) {
				boundary_nodes.push_back(index(i, j));
                boundary_bool.push_back(true);
			} else {
				inner_nodes.push_back(index(i, j));
                boundary_bool.push_back(false);
			}
            nodes.push_back(neighbour);


            if(i < N - 1){
                neighbours.push_back(node + N);
                weights.push_back(0.125);
                neighbour++;
            }


            if(i > 0){
                neighbours.push_back(node - N);
                weights.push_back(0.125);
                neighbour++;
            }


            if(j > 0){
                neighbours.push_back(node - 1);
                weights.push_back(0.125);
                neighbour++;
            }

            if(i < N - 1){
                neighbours.push_back(node + 1);
                weights.push_back(0.125);
                neighbour++;
            }

            if(i < N - 1 && j < N - 1){
                neighbours.push_back(node + N + 1);
                weights.push_back(0.0625);
                neighbour++;
            }

            if(i < N - 1 && j > 0){
                neighbours.push_back(node + N - 1);
                weights.push_back(0.0625);
                neighbour++;
            }

            if(i > 0 && j > 0){
                neighbours.push_back(node - N - 1);
                weights.push_back(0.0625);
                neighbour++;
            }

            if(i > 0 && j < N - 1){
                neighbours.push_back(node - N + 1);
                weights.push_back(0.0625);
                neighbour++;
            }

            node++;

		}
	}


	minimal = !is_2nplusone(N);
}


int Graph::numel() {
	return N*N;
}


const std::vector<Index>& Graph::get_inner_nodes() {
	return inner_nodes;
}


std::tuple<int, int, int, int> Graph::get_cardinal_neighbours(Index i) {
	return {
		/* nord  */   neighbours.at(nodes.at(i))
		/* sud   */ , neighbours.at(nodes.at(i) + 1)
		/* ovest */ , neighbours.at(nodes.at(i) + 2)
		/* est   */ , neighbours.at(nodes.at(i) + 3)
	};
}


std::tuple<int, int, int, int> Graph::get_diagonal_neighbours(Index i) {
	return {
		/* nord est   */   neighbours.at(nodes.at(i) + 4)
		/* nord ovest */ , neighbours.at(nodes.at(i) + 5)
		/* sud  ovest */ , neighbours.at(nodes.at(i) + 6)
		/* sud  est   */ , neighbours.at(nodes.at(i) + 7)
	};
}


Index Graph::index(int i, int j) {
	return i + j * N;
}


std::pair<int, int> Graph::inverse_index(Index i) {
	return {i % N, i / N};
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
    return nodes.at(i + 1) - nodes.at(i);
}


void Graph::project_on_coarse(Graph &coarse, const std::vector<double> &u, std::vector<double> &v) {
	for (auto &x : v) {
		x = 0.0;
	}

	for (int i = 2; i < (N-1); i = i + 2) {
		for (int j = 2; j < (N-1); j = j + 2) {
			int ii = i / 2;
			int jj = j / 2;


			v[coarse.index(ii, jj)] = 0.25   * (u[index(i,j)]);

            for(int z = 0; z < num_neighbours(index(i,j)); z++){

                v[coarse.index(ii, jj)] = u[neighbours.at(nodes.at(index(i,j)+z))] * weights.at(nodes.at(index(i,j)+z));
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
	for (int i = 0; i < numel(); ++i) {
		u[i] = 0.0;
	}
}


void Graph::evaluate_forcing_term(std::vector<double> &b, double (*f)(double x, double y)) {
	for (int j = 0; j < N; ++j) {
		for (int i = 0; i < N; ++i) {
			double x = x_corner + i * h * width;
			double y = y_corner + j * h * height;


			b[index(i,j)] = h * h * f(x,y);
		}
	}
}


void Graph::evaluate_boundary_conditions(std::vector<double> &u, double (*f)(double x, double y)) {
	for (Index i : boundary_nodes) {
		const auto [j, k] = inverse_index(i);
		const double x = x_corner + j * h * width;
		const double y = y_corner + k * h * height;

		u[i] = f(x,y);
	}


	for (Index i : inner_nodes) {
		u[i] = 0.0;
	}
}


void Graph::evaluate_function(std::vector<double> &u, double (*f)(double x, double y)) {
	for (int j = 0; j < N; ++j) {
		for (int i = 0; i < N; ++i) {
			double x = x_corner + i * h * width;
			double y = y_corner + j * h * height;


			u[index(i,j)] = f(x,y);
		}
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
	for (Index i : inner_nodes) {
		std::cout << i << std::endl;
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
