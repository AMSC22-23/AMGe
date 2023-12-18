#include <iostream>
#include<chrono>
#include "Multigrid.hpp"

#define N 257
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


    std::chrono::high_resolution_clock::time_point t0;
	std::chrono::high_resolution_clock::time_point t1;
	std::int64_t time_elapsed; 
    
    
    Lattice mesh = Lattice(0.0, 0.0, 1.0, 1.0, N);
    
    
    tests.push_back(Multigrid(mesh, PRE_STEP, POST_STEP, 2));
    tests.push_back(Multigrid(mesh, PRE_STEP, POST_STEP, 3));
    tests.push_back(Multigrid(mesh, PRE_STEP, POST_STEP, 4));


    std::vector<double> b(mesh.numel());
    std::vector<double> u(mesh.numel());
    std::vector<double> r(mesh.numel());


    mesh.evaluate_forcing_term(b, f);
    mesh.evaluate_boundary_conditions(u, g);

    
    int lvl = 2;
    t0 = std::chrono::high_resolution_clock::now();
    

    do{
        gseidel(mesh, u, b);
        residual(mesh, u, b, r);
    }while(norm(r) > 1.e-14); 


    t1 = std::chrono::high_resolution_clock::now();
    time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout<< "Gauss-Seidel converge in " << time_elapsed << " millisecondi" << std::endl;



    for(Multigrid m : tests){
        mesh.evaluate_boundary_conditions(u, g);
        t0 = std::chrono::high_resolution_clock::now();

        do{
            m.step(u, b);
            residual(mesh, u, b, r);
       }while(norm(r) > 1.e-14);

        t1 = std::chrono::high_resolution_clock::now();
        time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        std::cout<< "Multigrid a " << lvl << " livelli converge in " << time_elapsed << " millisencondi" << std::endl;
        lvl++;
    }


    return 0;
}