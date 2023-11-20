#include "Point.hpp"
#include "FunctionNode.hpp"
#include <string>


class Mesh : public FunctionNode<Point> {   //map from mesh to omega
	public:
		Mesh(double _x, double _y, double _width, double _height, unsigned int _Nx, unsigned int _Ny) : FunctionNode() {
			x      = _x;
			y      = _y;
			width  = _width;
			height = _height;
			Nx     = _Nx;
			Ny     = _Ny;
		};


		Point getValue(int i, int j) override {
			return Point(
					 x + i * getDiscretizationStepX()
					,y + j * getDiscretizationStepY()
				    );
		}


		unsigned int getDimensionX() {
			return Nx;
		}

		unsigned int getDimensionY() {
			return Ny;
		}


		double getDiscretizationStepX() {
			return width / static_cast<double>(Nx);
		}

		double getDiscretizationStepY() {
			return height / static_cast<double>(Ny);
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


		void save(std::ostream &);


		~Mesh() override = default;


	private:
		unsigned int Nx, Ny;
		double x,y,width,height;
};


void Mesh::save(std::ostream &out){
	{
		out
			<< "# name: xx"            << std::endl
			<< "# type: matrix"        << std::endl
			<< "# rows: 1"             << std::endl
			<< "# columns: "    << Nx  << std::endl;


		for (int i = 0; i < Nx; ++i) {
			out << getValue(i, 0).x << " ";
		}


		out << std::endl;
	}


	out << std::endl;


	{
		out
			<< "# name: yy"            << std::endl
			<< "# type: matrix"        << std::endl
			<< "# rows: 1"             << std::endl
			<< "# columns: "    << Ny  << std::endl;


		for (int j = 0; j < Ny; ++j) {
			out << getValue(0, j).y << " ";
		}


		out << std::endl;
	}
}
