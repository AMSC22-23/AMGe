#include <iostream>
#include <cmath>

#define N 10
#define N_INTERNAL (N-2)
#define L 10
#define H 10

template<class T>
class FunctionNode {		//abstract class for mapping
public:
    FunctionNode() = default;
    virtual T getValue(int i,int j)  = 0;
    virtual ~FunctionNode() = default;
};


class  Mesh: public FunctionNode<std::pair<double,double>> {   //map from mesh to omega
public:
    Mesh() : FunctionNode() {
        hx=static_cast<double>(L)/static_cast<double>(N);
        hy=static_cast<double>(H)/static_cast<double>(N);
    };
    std::pair<double, double> getValue(int i,int j) override{
		std::pair<double,double> a;
		a.first=hx*static_cast<double>(i);
		a.second=hy*static_cast<double>(j);
		return a;
	}
	double getX(){
		return hx;
	}
    ~Mesh() override = default;
    private:
        double hx;
        double hy;
};




class  ComputeFunctionNode : public FunctionNode<double> {  //compute mapping of a function from mesh index to R
public:
    ComputeFunctionNode(Mesh a,double (*func)(double, double)) : FunctionNode() {
		mesh=a;
		fa=func;
	};
	double getValue(int i,int j) override{
		std::pair<double,double> index;
		index=mesh.getValue(i,j);
		return fa(index.first,index.second);
	}
    virtual ~ComputeFunctionNode() override = default;
	private: 
	Mesh mesh;
	double (*fa)(double, double);

};





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


/*int index(int i, int j) {
	return i * (N_INTERNAL) + j;
}*/


int local_index(int i, int j) {			//index matrix in a vector
	return (i-1) * N_INTERNAL + j - 1;
}


/*double bordo(int i, int j) {
	return 1.0;
}*/

double g(double x,double y){  //function buondary conditions
	return 1.0;
}

double f(double x,double y){  //main function
	return 0.0;
}

/*double f(int i, int j) {
	return 0.0;
}*/


int main (int argc, char *argv[]) {
	double A[(N-2)*(N-2)][(N-2)*(N-2)] = {0.0};
	double b[(N-2)*(N-2)]      = {0.0};

	double h = 1.0 / static_cast<double>(N-1);
	Mesh a;
	ComputeFunctionNode bordo(a,g);
	ComputeFunctionNode funzione(a,f);

	int riga = 0;

	for (int i = 1; i < N-1; ++i) {
		for (int j = 1; j < N-1; ++j) {
			// sono un nodo generico
			A[riga][local_index(i,  j)] = -4.0;
			//b[riga] = f(i,j);
			b[riga]=funzione.getValue(i,j);


			if (i != 1) {
				A[riga][local_index(i-1,j)] =  1.0;
			}
			else {
				//b[riga] -= bordo(i-1,j);
				b[riga] -= bordo.getValue(i-1,j);
			}

			if (i != N-2) {
				A[riga][local_index(i+1,j)] =  1.0;
			}
			else {
				b[riga] -= bordo.getValue(i+1,j);
			}

			if (j != 1) {
				A[riga][local_index(i,j-1)] =  1.0;
			}
			else {
				b[riga] -= bordo.getValue(i,j-1);
			}

			if (j != N-2) {
				A[riga][local_index(i,j+1)] =  1.0;
			}
			else {
				b[riga] -= bordo.getValue(i,j+1);
			}

			++riga;
		}
	}


	//print_matrix(A);     //print
	//print_vector(b);




	return 0;
}
