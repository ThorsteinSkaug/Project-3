#include "PenningTrap.cpp"
#include "Particle.cpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>




int main(){
    int t = 100;
    double k_e =  1.38935333*pow(10,5);
    double T = 9.64852558 * 10;
    double V = 9.64852558 * pow(10,7);
    double B_0 = 1. * T;
    double V_0 = 10 * V;
    double d = pow(10,4);
    double V_d_ratio = 9.65;

    int q = 1;
    double m = 4.077;
    arma::vec r = {1., 1., 1.};
    arma::vec v = {2., 1., 3.};
    vector<Particle> pl;
    Particle singly_charged_Calcium = Particle(q, m, r, v);
    PenningTrap trap = PenningTrap(B_0, V_0, d, pl);

    std::cout << singly_charged_Calcium.charge << singly_charged_Calcium.mass << singly_charged_Calcium.position;

    trap.evolve_forward_Euler(0.1);
    trap.evolve_forward_Euler(0.1);
    trap.evolve_forward_Euler(0.1);
    trap.evolve_forward_Euler(0.1);
    std::cout << trap.particle_l[0].mass;
    return 0;
}
