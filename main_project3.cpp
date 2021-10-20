#include "PenningTrap.cpp"
#include "Particle.cpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void analytical_solution(arma::vec& x, arma::vec& y, arma::vec& z, double dt, double w_z, double w_p, double w_m, double A_p, double A_m){
  for(int i=1; i<x.size(); i++){
    x[i] = A_p*cos(w_m*i*dt) + A_p*cos(w_m*i*dt);
    y[i] = A_p*sin(w_m*i*dt) + A_p*sin(w_m*i*dt);
    z[i] = 1./2 * (z[0]*exp(w_z*(i*dt)) + z[0]*exp(-w_z*(i*dt)));
  }
}



int main(){
    int t = 100;
    double k_e =  1.38935333*pow(10,5);
    double T = 9.64852558 * 10;
    double V = 9.64852558 * pow(10,7);
    double B_0 = 1. * T;
    double V_0 = 10 * V;
    double d = pow(10,4);
    double V_d_ratio = 9.65;

    //One particle
    int q = 1.;
    double m = 40.078;
    arma::vec r = {1., 0., 1.};
    arma::vec v = {0., 1., 0.};
    vector<Particle> pl;
    Particle singly_charged_Calcium = Particle(q, m, r, v);
    PenningTrap trap = PenningTrap(B_0, V_0, d, pl);
    trap.add_particle(singly_charged_Calcium);
    double dt = 0.001;


    //Start analytical solution part
    double w_0 = singly_charged_Calcium.charge*B_0/singly_charged_Calcium.mass;
    double w_z = sqrt(2*singly_charged_Calcium.charge*V_0/(singly_charged_Calcium.mass*pow(d,2)));

    double w_p = (w_0 + sqrt(pow(w_0,2)-2*pow(w_z,2)))/2;
    double w_m = (w_0 - sqrt(pow(w_0,2)-2*pow(w_z,2)))/2;

    arma::vec x(int(100/dt)+1);
    x(0) = trap.particle_l[0].position[0];
    arma::vec y(int(100/dt)+1);
    y(0) = trap.particle_l[0].position[1];
    arma::vec z(int(100/dt)+1);
    z(0) = trap.particle_l[0].position[2];

    double A_p = (1+w_m*1.)/(w_m - w_p);
    double A_m = (1+w_p*1.)/(w_m - w_p);

    analytical_solution(x, y, z, dt, w_z, w_p, w_m, A_p, A_m);

    std::ofstream myfile_an;
    myfile_an.open("coordinates_an.txt");
    for(int i=0; i<x.size(); i++){
      myfile_an << std::scientific << dt*(i) << " " << std::scientific << x[i] << " " << std::scientific << y[i] <<" " << std::scientific << z[i] << "\n";
    }
    //End analytical solution part


    //Start RK4 solution part
    std::ofstream myfile;
    myfile.open("coordinates_rk.txt");
    myfile << std::scientific << 0 << " " << std::scientific << trap.particle_l[0].position[0] << " " << std::scientific << trap.particle_l[0].position[1] <<" " << std::scientific << trap.particle_l[0].position[2] << "\n";
    for(int i=0; i<int(100/dt); i++){
      trap.evolve_RK4(dt);
      myfile << std::scientific << dt*(i+1) << " " << std::scientific << trap.particle_l[0].position[0] << " " << std::scientific << trap.particle_l[0].position[1] <<" " << std::scientific << trap.particle_l[0].position[2] << "\n"; //Write the number of iteration needed for convergence to file
    }
    //End RK4 solution part

    //Start euler solution part
    q = 1.;
    m = 40.078;
    r = {1., 0., 1.};
    v = {0., 1., 0.};
    vector<Particle> pl2;
    Particle singly_charged_Calcium2 = Particle(q, m, r, v);
    PenningTrap trap2 = PenningTrap(B_0, V_0, d, pl2);
    trap2.add_particle(singly_charged_Calcium2);

    std::ofstream myfile2;
    myfile2.open("coordinates_euler.txt");
    myfile2 << std::scientific << 0 << " " << std::scientific << trap2.particle_l[0].position[0] << " " << std::scientific << trap2.particle_l[0].position[1] <<" " << std::scientific << trap2.particle_l[0].position[2] << "\n";
    dt = 0.001;
    for(int i=0; i<int(100/dt); i++){
      trap2.evolve_forward_Euler(dt);
      myfile2 << std::scientific << dt*(i+1) << " " << std::scientific << trap2.particle_l[0].position[0] << " " << std::scientific << trap2.particle_l[0].position[1] <<" " << std::scientific << trap2.particle_l[0].position[2] << "\n"; //Write the number of iteration needed for convergence to file
    }
    //End euler solution part




    //Two particles:
    arma::vec r2 = {2., 0.3, 1.2};
    arma::vec v2 = {-0.3, 1., 0.2};
    std::vector<Particle> pl2particles;
    Particle singly_charged_Calcium_nr2 = Particle(q, m, r2, v2);
    PenningTrap trap2particles = PenningTrap(B_0, V_0, d, pl2particles);
    trap2particles.add_particle(singly_charged_Calcium);
    trap2particles.add_particle(singly_charged_Calcium_nr2);

    std::ofstream myfile2particles;
    myfile2particles.open("coordinates_rk2particles.txt");
    myfile2particles << std::scientific << 0 << " " << std::scientific << trap2particles.particle_l[0].position[0] << " " << std::scientific << trap2particles.particle_l[0].position[1] <<" " << std::scientific << trap2particles.particle_l[0].position[2] << " " << std::scientific << trap2particles.particle_l[1].position[0] << " " << std::scientific << trap2particles.particle_l[1].position[1] <<" " << std::scientific << trap2particles.particle_l[1].position[2] << "\n";
    for(int i=0; i<int(100/dt); i++){
      trap2particles.evolve_RK4(dt);
      myfile2particles << std::scientific << dt*(i+1) << " " << std::scientific << trap2particles.particle_l[0].position[0] << " " << std::scientific << trap2particles.particle_l[0].position[1] <<" " << std::scientific << trap2particles.particle_l[0].position[2] << " " << std::scientific << trap2particles.particle_l[1].position[0] << " " << std::scientific << trap2particles.particle_l[1].position[1] <<" " << std::scientific << trap2particles.particle_l[1].position[2] << "\n"; //Write the number of iteration needed for convergence to file
    }




    return 0;
}
