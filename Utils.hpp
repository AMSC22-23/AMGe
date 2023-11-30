#ifndef __UTILS_HPP__
#define __UTILS_HPP__


double norm(std::vector<double> &v){
	double result = 0.0;

	for (auto x : v) {
		result += x * x;
	}

	return std::sqrt(result);
}


void export_to_matlab(const char *name, std::vector<double> &v, std::ostream &out = std::cout){
	out
		<< "# name: "      << name     << std::endl
		<< "# type: matrix"            << std::endl
		<< "# rows: 1"                 << std::endl
		<< "# columns: "   << v.size() << std::endl;


	for (auto x : v) {
		out << x << ' ';
	}


	out << '\n' << '\n';
}

bool is_2nplusone(int n){
	int index=0;

	while(pow(2,index)+1<=n){
		if(pow(2,index)+1==n){
			return true;
		}
		index++;
	}
	return false;
}



#endif
