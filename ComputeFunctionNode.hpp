
#include "Mesh.hpp"

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
