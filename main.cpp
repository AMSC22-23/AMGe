#include <iostream>

#define N 10
#define N_INTERNAL (N-2)

void print_matrix(double A[N_INTERNAL*N_INTERNAL][N_INTERNAL*N_INTERNAL]){
	for(int i = 0; i < N_INTERNAL*N_INTERNAL; ++i){
		for(int j = 0; j < N_INTERNAL*N_INTERNAL; ++j){
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


void print_vector(double b[N_INTERNAL*N_INTERNAL]) {
		for(int i = 0; i < N_INTERNAL*N_INTERNAL; ++i){
			std::cout << b[i] << std::endl;
		}
}


int index(int i, int j) {
	return i * (N_INTERNAL) + j;
}


int local_index(int i, int j) {
	return (i-1) * N_INTERNAL + j - 1;
}


double bordo(int i, int j) {
	return 1.0;
}


double f(int i, int j) {
	return 0.0;
}


int main (int argc, char *argv[]) {
	double A[(N-2)*(N-2)][(N-2)*(N-2)] = {0.0};
	double b[(N-2)*(N-2)]      = {0.0};

	double h = 1.0 / static_cast<double>(N-1);



	int riga = 0;
	for (int i = 1; i < N-1; ++i) {
		for (int j = 1; j < N-1; ++j) {
			// sono un nodo generico
			A[riga][local_index(i,  j)] = -4.0;
			b[riga] = f(i,j);


			if (i != 1) {
				A[riga][local_index(i-1,j)] =  1.0;
			}
			else {
				b[riga] -= bordo(i-1,j);
			}

			if (i != N-2) {
				A[riga][local_index(i+1,j)] =  1.0;
			}
			else {
				b[riga] -= bordo(i+1,j);
			}

			if (j != 1) {
				A[riga][local_index(i,j-1)] =  1.0;
			}
			else {
				b[riga] -= bordo(i,j-1);
			}

			if (j != N-2) {
				A[riga][local_index(i,j+1)] =  1.0;
			}
			else {
				b[riga] -= bordo(i,j+1);
			}

			++riga;
		}
	}


	//print_matrix(A);
	print_vector(b);




	return 0;
}
