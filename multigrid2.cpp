
//@note: is this file even used? if not just remove it, it stays in git history
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <vector>
#include <tuple>
#include <stdio.h>


using Index = int;


class Lattice {
public:
	Lattice(double x_corner_, double y_corner_, double width_, double height_, unsigned int N_) {
		x_corner = x_corner_;
		y_corner = y_corner_;
		width = width_;
		height = height_;
		N = N_;
		h = width / (N-1);


		for (int i = 1; i < N-1; ++i) {
			for (int j = 1; j < N-1; ++j) {
				inner_nodes.push_back(index(i, j));
			}
		}
	}


	double get_step() {
		return h;
	}


	int numel() {
		return N*N;
	}


	const std::vector<Index>& get_inner_nodes() {
		return inner_nodes;
	}


	std::tuple<int, int, int, int> get_cardinal_neighbours(Index i) {
		return {
			  i + N
			, i - N
			, i - 1
			, i + 1
		};
	}

	std::tuple<int, int, int, int> get_diagonal_neighbours(Index i) {
		return {
			  i + N + 1
			, i + N - 1
			, i - N - 1
			, i - N + 1
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

				v[coarse.index(ii, jj)] =
					  0.25 * u[index(i,j)]
					+ 0.125 * (u[index(i+1,j)] + u[index(i,j+1)] + u[index(i-1,j)] + u[index(i,j-1)])
					+ 0.0625 * (u[index(i+1,j+1)] + u[index(i-1,j+1)] + u[index(i-1,j-1)] + u[index(i+1,j-1)]);
			}
		}
	}


	void interpolate_on_fine(Lattice &coarse, std::vector<double> &u, const std::vector<double> &v) {
		for (auto &x : u) {
			x = 0.0;
		}


		for (int i = 1; i < (N-1); ++i) {
			for (int j = 1; j < (N-1); ++j) {
				int ii = i / 2;
				int jj = j / 2;

				if ((i % 2 == 0) and (j % 2 == 0)) {
					u[index(i,j)] = v[coarse.index(ii, jj)];
				}
				else {
					u[index(i,j)] =
						0.25 * (v[coarse.index(ii+1,jj)] + v[coarse.index(ii,jj+1)] + v[coarse.index(ii-1,jj)] + v[coarse.index(ii,jj-1)])
						+ 0.125 * (v[coarse.index(ii+1,jj+1)] + v[coarse.index(ii-1,jj+1)] + v[coarse.index(ii-1,jj-1)] + v[coarse.index(ii+1,jj-1)]);
				}
			}
		}
	}
/*
def prolong(x):
    n = x.shape[0]
    N = (n - 1) * 2 + 1
    y = np.zeros((N, N))
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            if (i % 2 == 0) and (j % 2 == 0):
                y[i, j] = x[i // 2, j // 2]
            else:
                ii, jj = i // 2, j // 2
                y[i, j] = +0.25 * (
                    x[ii + 1, jj] + x[ii, jj + 1] + x[ii - 1, jj] + x[ii, jj - 1]
                ) + 0.125 * (
                    x[ii + 1, jj + 1]
                    + x[ii - 1, jj + 1]
                    + x[ii - 1, jj - 1]
                    + x[ii + 1, jj - 1]
                )
    return y
    */


	void evaluate_zero(std::vector<double> &u) {
		for (int i = 0; i < numel(); ++i) {
			u[i] = 0.0;
		}
	}


	void evaluate_function(std::vector<double> &u, double (*f)(double x, double y)) {
		double h = get_step();

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
};


double norm(const std::vector<double> &u) {
	double sum = 0.0;

	for (auto x : u) {
		sum += x * x;
	}

	return std::sqrt(sum);
}


void residual(Lattice &mesh, const std::vector<double> &u, const std::vector<double> &b, std::vector<double> &r) {
	const double h = mesh.get_step();

	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		// implementazione originale
		//r[i] = h*h*b[i] + u[nord] + u[sud] + u[ovest] + u[est] - 4.0 * u[i];
		//

		// implementazione di caldana
		r[i] = 4.0 * u[i] - u[nord] - u[sud] - u[ovest] - u[est] - b[i];
	}
}


double error(const std::vector<double> &u, const std::vector<double> &v) {
	std::vector<double> err(u.size());

	for (int i = 0; i < u.size(); ++i) {
		err[i] = u[i] - v[i];
	}

	return norm(err);

}


void jacobi(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	const double h = mesh.get_step();


	// occhio al bordo
	for (int i = 0; i < mesh.numel(); ++i) {
		u[i] = old[i];
	}


	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (h*h*b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}
}


void gseidel(Lattice &mesh, std::vector<double> &u, std::vector<double> &b) {
	const double h = mesh.get_step();

	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		//u[i] = 0.25 * (h*h*b[i] + u[nord] + u[sud] + u[ovest] + u[est]);
		u[i] = 0.25 * (b[i] + u[nord] + u[sud] + u[ovest] + u[est]);
	}
}


/*
double g(double x,double y){
	return std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
}


double f(double x,double y){
	return 8.0 * M_PI * M_PI * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
}
*/


double g(double x,double y){
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y){
	return -(x * x * x + y * y * y);
}

/*double g(double x,double y){
	//return std::exp(x) * std::exp(-2.0 * y);
	//return 1;
}*/


/*double f(double x,double y){
	//return -5.0 * std::exp(x) * std::exp(-2.0 * y);
	//return 0;
}*/


void set_initial_guess(Lattice &mesh, std::vector<double> &u, double (*f)(double x, double y)) {
	mesh.evaluate_function(u, f);

	for (Index i : mesh.get_inner_nodes()) {
		u[i] = 0.0;
	}
}


int main() {
	Lattice fine(0.0, 0.0, 1.0, 1.0, 33);
	Lattice coarse = fine.build_coarse();


	std::vector<double> u(fine.numel());
	std::vector<double> exact(fine.numel());
	std::vector<double> u_old(fine.numel());

	std::vector<double> b(fine.numel());
	std::vector<double> b2h(coarse.numel());

	std::vector<double> r(fine.numel());
	std::vector<double> r2h(coarse.numel());

	std::vector<double> e(fine.numel());
	std::vector<double> e2h(coarse.numel());
	std::vector<double> e2h_old(coarse.numel());


	set_initial_guess(fine, u, g);

	fine.evaluate_function(exact, g);

	fine.evaluate_function(b, f);
	double h = fine.get_step();

	for (auto &x : b) {
		x = h * h * x;
	}


	for (int i = 0; i < 1000; ++i) {
		for (int pre = 0; pre < 10; ++pre) {
			gseidel(fine, u, b);
		}

		// multigrid
		residual(fine, u, b, r);
		fine.project_on_coarse(coarse, r, b2h);

		coarse.evaluate_zero(e2h);
		coarse.evaluate_zero(e2h_old);

		for (int coarse_it = 0; coarse_it < 300; ++coarse_it) {
			gseidel(coarse, e2h, b2h);
		}

		fine.interpolate_on_fine(coarse, e, e2h);

		for (Index i : fine.get_inner_nodes()) {
			u[i] -= e[i];
		}

		for (int post = 0; post < 10; ++post) {
			gseidel(fine, u, b);
		}

		residual(fine, u, b, r);
		std::cout <<i<<" "<< norm(r) << std::endl;
	}


	return 0;
}