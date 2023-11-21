#include "Mesh.hpp"
#include "Point.hpp"


class ComputeFunctionNode : public FunctionNode<double> {  //compute mapping of a function from mesh index to R
public:
    ComputeFunctionNode(Mesh *_M, double (*_f)(double, double)) : FunctionNode() {
	M = _M;
	f = _f;
    };


    double getValue(int i, int j) override {
	Point p = M->getValue(i, j);
	return f(p.x, p.y);
    }


    Mesh* getMesh(){
        return M;
    }


    virtual ~ComputeFunctionNode() override = default;


private:
    Mesh *M;
    double (*f)(double, double);
};
