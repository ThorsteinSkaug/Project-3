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
  arma::arma_rng::set_seed(74);

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


/*
  std::ofstream myfile;
  myfile.open("resonance.txt");
  bool particle_interaction = false;
  double dt = 0.1;
  double w_V_step = 0.02;
  for(double w=0.20+w_V_step; w<2.5; w += w_V_step){
    cout << w << "\n";
    PenningTrap trap1 = PenningTrap(B_0, V_0, d, pl, f[0], w);
    PenningTrap trap2 = PenningTrap(B_0, V_0, d, pl, f[1], w);
    PenningTrap trap3 = PenningTrap(B_0, V_0, d, pl, f[2], w);
    std::vector<Particle> original_pl = pl;
    for(int j=0; j<int(500/dt); j++){
      trap1.evolve_RK4(dt, particle_interaction);
      trap2.evolve_RK4(dt, particle_interaction);
      trap3.evolve_RK4(dt, particle_interaction);
    }
    int count_part1 = count_part_inside(trap1.particle_l, d);
    int count_part2 = count_part_inside(trap2.particle_l, d);
    int count_part3 = count_part_inside(trap3.particle_l, d);

    myfile << w << " " << count_part1 << " " << count_part2 << " " << count_part3 << "\n";
    pl = original_pl;
  }
*/

  std::ofstream myfile2;
  myfile2.open("resonance_zoom.txt");
  bool particle_interaction = false;
  double dt = 0.1;
  for(double w=0.37; w<0.52; w += 0.005){
    cout << w << "\n";
    PenningTrap trap_without = PenningTrap(B_0, V_0, d, pl, f[1], w);
    PenningTrap trap_with = PenningTrap(B_0, V_0, d, pl, f[1], w);
    std::vector<Particle> original_pl = pl;
    for(int j=0; j<int(500/dt); j++){
      trap_without.evolve_RK4(dt, false);
      trap_with.evolve_RK4(dt, true);
    }
    int count_part1 = count_part_inside(trap_without.particle_l, d);
    int count_part2 = count_part_inside(trap_with.particle_l, d);

    myfile2 << w << " " << count_part1 << " " << count_part2 << "\n";
    pl = original_pl;
  }
}
