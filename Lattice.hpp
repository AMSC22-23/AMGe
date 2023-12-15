#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__


#include <iostream>
#include <vector>
#include <tuple>
#include "Utils.hpp"


using Index = int;


class Lattice {
public:
	Lattice(double x_corner_, double y_corner_, double width_, double height_, unsigned int N_) {
		x_corner = x_corner_;
		y_corner = y_corner_;
		width    = width_;
		height   = height_;


		if (width != height) {
			fprintf(stderr, "Mesh has to be square\n");
			exit(EXIT_FAILURE);
		}


		N        = N_;
		h        = width / (N-1);


		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (i == 0 || i == N-1 || j == 0 || j == N-1) {
					boundary_nodes.push_back(index(i, j));
				} else {
					inner_nodes.push_back(index(i, j));
				}
			}
		}


		minimal = !is_2nplusone(N) || width != height;
	}


	int numel() {
		return N*N;
	}


	const std::vector<Index>& get_inner_nodes() {
		return inner_nodes;
	}


	std::tuple<int, int, int, int> get_cardinal_neighbours(Index i) {
		return {
			/* nord  */   i + N
			/* sud   */ , i - N
			/* ovest */ , i - 1
			/* est   */ , i + 1
		};
	}


	std::tuple<int, int, int, int> get_diagonal_neighbours(Index i) {
		return {
			/* nord est   */   i + N + 1
			/* nord ovest */ , i + N - 1
			/* sud  ovest */ , i - N - 1
			/* sud  est   */ , i - N + 1
		};
	}


	Index index(int i, int j) {
		return i + j * N;
	}


	std::pair<int, int> inverse_index(Index i) {
		return {i % N, i / N};
	}


	Lattice build_coarse() {
		return Lattice(
			  x_corner
			, y_corner
			, width
			, height
			, 1 + ((N - 1) / 2)
		);
	}


	void project_on_coarse(Lattice &coarse, const std::vector<double> &u, std::vector<double> &v) {
		for (auto &x : v) {
			x = 0.0;
		}

		for (int i = 2; i < (N-1); i = i + 2) {
			for (int j = 2; j < (N-1); j = j + 2) {
				int ii = i / 2;
				int jj = j / 2;

				const auto [nord, sud, ovest, est] = get_cardinal_neighbours(index(i,j));
				const auto [nordest, nordovest, sudovest, sudest] = get_diagonal_neighbours(index(i,j));

				v[coarse.index(ii, jj)] =
					  0.25   * (u[index(i,j)])
					+ 0.125  * (u[nord] + u[sud] + u[ovest] + u[est])
					+ 0.0625 * (u[nordest] + u[nordovest] + u[sudovest] + u[sudest]);
			}
		}
	}


	void interpolate_on_fine(Lattice &coarse, std::vector<double> &u_fine, const std::vector<double> &u_coarse) {
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


	void evaluate_zero(std::vector<double> &u) {
		for (int i = 0; i < numel(); ++i) {
			u[i] = 0.0;
		}
	}


	void evaluate_forcing_term(std::vector<double> &b, double (*f)(double x, double y)) {
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < N; ++i) {
				double x = x_corner + i * h * width;
				double y = y_corner + j * h * height;


				b[index(i,j)] = h * h * f(x,y);
			}
		}
	}


	void evaluate_boundary_conditions(std::vector<double> &u, double (*f)(double x, double y)) {
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


	void evaluate_function(std::vector<double> &u, double (*f)(double x, double y)) {
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < N; ++i) {
				double x = x_corner + i * h * width;
				double y = y_corner + j * h * height;


				u[index(i,j)] = f(x,y);
			}
		}
	}


	void print_vector(const std::vector<double> &u, const char *name) {
		std::cout << name << std::endl;

		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < N; ++i) {
				printf("%e\t", u[index(i, j)]);
			}

			std::cout << std::endl;
		}

		std::cout << std::endl;
	}


	void test_constructor() {
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

	
	void test_inner_nodes() {
		for (Index i : inner_nodes) {
			std::cout << i << std::endl;
		}
	}


	void test_index() {
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < N; ++i) {
				std::cout << index(i,j) << std::endl;
			}
		}
	}


	void test_inverse_index() {
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				const auto [k, l] = inverse_index(index(i,j));

				if (k != i || l != j) {
					std::cout << "errore" << std::endl;
				}
			}
		}
	}


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
