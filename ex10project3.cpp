#include "PenningTrap.cpp"
#include "Particle.cpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>



int count_part_outside(std::vector<Particle> particle_l, double d){
  int count = 0;
  for(int i = 0; i<particle_l.size(); i++){
    if(sqrt(pow(particle_l[i].position[0], 2) + pow(particle_l[i].position[1], 2) + pow(particle_l[i].position[2], 2)) < d){
      count ++;
    }
  }
  return count;
}


void fill_random_particles(std::vector<Paricles>& particle_l){
  for(int i=0; i<particle_l.size(); i++){
    particle_l(i).position = vec(3).randn() * 0.1 * trap.d_;
    particle_l(i).velocity = vec(3).randn() * 0.1 * trap.d_;
  }
}


int main(){
  arma_rng::set_seed(4)

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
    pl(part) = Particle(q, m, r, v)
  }
  PenningTrap trap = PenningTrap(B_0, V_0, d, pl);

}
