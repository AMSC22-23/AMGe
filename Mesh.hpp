#include "Point.hpp"
#include "FunctionNode.hpp"


class Mesh : public FunctionNode<Point> {   //map from mesh to omega
public:
    Mesh(double _x, double _y, double _width, double _height, unsigned int _N) : FunctionNode() {
	x      = _x;
	y      = _y;
	width  = _width;
	height = _height;
	N      = _N;
    };


    Point getValue(int i, int j) override {
	return Point(
	     x + i * width  / N
	    ,y + j * height / N
	);
    }


    unsigned int getDimension () {
	return N;
    }


    double getDiscretizationStepX() {
	return width / static_cast<double>(N);
    }

    double getDiscretizationStepY() {
	return height / static_cast<double>(N);
    }

    double getWidth(){
        return width;
    }
    
    double getHeight(){
        return height;
    }

    double getXCenter(){
        return x;
    }

    double getYCenter(){
        return y;
    }


    ~Mesh() override = default;


private:
    unsigned int N;
    double x,y,width,height;
};
