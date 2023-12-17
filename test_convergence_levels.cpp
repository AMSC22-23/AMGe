#include <iostream>
#include "Multigrid.hpp"

#define N 65
#define PRE_STEP 10
#define POST_STEP 10



double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}



int main(){ 
    std::vector<Multigrid> tests;   
    
    
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, N);
    
    
    tests.push_back(Multigrid(mesh, PRE_STEP, POST_STEP, 2));
    tests.push_back(Multigrid(mesh, PRE_STEP, POST_STEP, 3));
    tests.push_back(Multigrid(mesh, PRE_STEP, POST_STEP, 4));


    std::vector<double> b(mesh.numel());
    std::vector<double> u(mesh.numel());
    std::vector<double> r(mesh.numel());


    mesh.evaluate_forcing_term(b, f);

    
    int lvl = 2;
    int iterations;


    for(Multigrid m : tests){
        mesh.evaluate_boundary_conditions(u, g);
        iterations = 0;

        do{
            iterations ++;
            m.step(u, b);
            residual(mesh, u, b, r);
       }while(norm(r) > 1.e-14);


        std::cout<< "Multigird a " << lvl << " livelli converge in " << iterations << " iterazioni" << std::endl;
        lvl ++;
    }


    return 0;
}