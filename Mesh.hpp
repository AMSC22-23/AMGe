#include "Point.hpp"
#include "FunctionNode.hpp"
#include <string>


class Mesh : public FunctionNode<Point> {   //map from mesh to omega
	public:
		Mesh(
			  double _x
			, double _y
			, double _width
			, double _height
			, unsigned int _Nx
			, unsigned int _Ny
		) :
			  x(_x)
			, y(_y)
			, width(_width)
			, height(_height)
			, Nx(_Nx)
			, Ny(_Ny)
			, hx(width / static_cast<double>(Nx))
			, hy(width / static_cast<double>(Ny))
		{}


		Point getValue(int i, int j) override {
			return Point(
					 x + i * hx
					,y + j * hy
				    );
		}


		int index(int i, int j) {
			return i * Ny + j;
		}


		void save(std::ostream &);
		~Mesh() override = default;


		const unsigned int Nx, Ny;
		const double x,y,width,height;
		const double hx,hy;
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
