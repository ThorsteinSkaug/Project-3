#include "PenningTrap.cpp"
#include "Particle.cpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int count_part_inside(std::vector<Particle> particle_l, double d){
  int count = 0;
  for(int i = 0; i<particle_l.size(); i++){
    if(sqrt(pow(particle_l[i].position[0], 2) + pow(particle_l[i].position[1], 2) + pow(particle_l[i].position[2], 2)) < d){
      count ++;
    }
  }
  return count;
}


void fill_random_particles(std::vector<Particle>& particle_l, double d){
  for(int i=0; i<particle_l.size(); i++){
    particle_l[i].position = arma::vec(3).randn() * 0.1 * d;
    particle_l[i].velocity = arma::vec(3).randn() * 0.1 * d;
  }
}


int main(){
  arma::arma_rng::set_seed(4);

  double k_e =  1.38935333*pow(10,5);
  double T = 9.64852558 * 10;
  double V = 9.64852558 * pow(10,7);
  double B_0 = 1. * T;
  double V_0 = 0.0025 * V;
  double d = 0.05*pow(10,4);
  int q = 1;
  double m = 40.078;
  arma::vec r = {1., 1., 1.};
  arma::vec v = {1., 1., 1.};

  int n_particles = 100;
  std::vector<Particle> pl;
  for(int part=0; part<n_particles;part++){
    pl.push_back(Particle(q, m, r, v));
  }
  fill_random_particles(pl, d);
  std::vector<double> f = {0.1, 0.4, 0.7};
  std::vector<double> w_V;

  for(int part=0; part<n_particles;part++){
    cout << pl[part].position << "\n";
  }

  bool particle_interaction = false;
  double dt = 0.1;
  for(double w=0.22; w<2.5; w += 0.02){
    PenningTrap trap = PenningTrap(B_0, V_0, d, pl, f[2], w);
    std::vector<Particle> original_pl = pl;
    for(int j=0; j<int(500/dt); j++){
      trap.evolve_RK4(dt, particle_interaction);
    }
    int count_part = count_part_inside(trap.particle_l, d);

    cout << count_part << "\n";
    pl = original_pl;
  }
}
