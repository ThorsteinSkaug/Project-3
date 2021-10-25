#include "PenningTrap.cpp"
#include "Particle.cpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;


//A function for counting the number of particles still inside the trap
//Take as its input variables a vector containing Particle objects and a double d which defines the size of the trap
int count_part_inside(std::vector<Particle> particle_l, double d){
  int count = 0; //Variable for counting
  for(int i = 0; i<particle_l.size(); i++){
    if(sqrt(pow(particle_l[i].position[0], 2) + pow(particle_l[i].position[1], 2) + pow(particle_l[i].position[2], 2)) < d){ //The condition for still being inside the trap
      count ++; //Add one to the count if condition is met
    }
  }
  return count;
}

//A function for randomizing the position and velocity of the particles in the trap
//Takes as its input a vector containing particle objects and the size of the trap d
//After running it has updated the position and velocity for every particle in the particle_l vector
void fill_random_particles(std::vector<Particle>& particle_l, double d){
  for(int i=0; i<particle_l.size(); i++){
    particle_l[i].position = arma::vec(3).randn() * 0.1 * d; //Randomize the position
    particle_l[i].velocity = arma::vec(3).randn() * 0.1 * d; //Randomize the velocity
  }
}


int main(){
  arma::arma_rng::set_seed(74); //Set a random seed

  // Defining constants and initial vectors
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

  int n_particles = 100; //Number of particles in the trap
  std::vector<Particle> pl;
  for(int part=0; part<n_particles;part++){
    pl.push_back(Particle(q, m, r, v)); //Add a new Particle to the trap
  }
  fill_random_particles(pl, d); //Fille the Particle vector pl with random positions and velocity for each particle
  std::vector<double> f = {0.1, 0.4, 0.7}; //The different amplitudes


  std::ofstream myfile;
  myfile.open("resonance.txt");
  bool particle_interaction = false; //Set particle interactions to false
  double dt = 0.1; //Set a step size
  double w_V_step = 0.02; //Define how fast we go through the different values for w_V
  for(double w=0.20+w_V_step; w<2.5; w += w_V_step){
    //Define some Penning traps
    PenningTrap trap1 = PenningTrap(B_0, V_0, d, pl, f[0], w);
    PenningTrap trap2 = PenningTrap(B_0, V_0, d, pl, f[1], w);
    PenningTrap trap3 = PenningTrap(B_0, V_0, d, pl, f[2], w);
    std::vector<Particle> original_pl = pl; //store the original particles in a seperate vector
    for(int j=0; j<int(500/dt); j++){
      //Evelove the three traps using runge-kutta
      trap1.evolve_RK4(dt, particle_interaction);
      trap2.evolve_RK4(dt, particle_interaction);
      trap3.evolve_RK4(dt, particle_interaction);
    }
    //Count the number of particles still inside the trap at the end
    int count_part1 = count_part_inside(trap1.particle_l, d);
    int count_part2 = count_part_inside(trap2.particle_l, d);
    int count_part3 = count_part_inside(trap3.particle_l, d);

    myfile << w << " " << count_part1 << " " << count_part2 << " " << count_part3 << "\n";
    pl = original_pl; //Reset the pl vector to the initial values
  }

  std::ofstream myfile2;
  myfile2.open("resonance_zoom.txt");
  dt = 0.1; //Step size
  for(double w=0.37; w<0.52; w += 0.005){
    //define traps
    PenningTrap trap_without = PenningTrap(B_0, V_0, d, pl, f[1], w);
    PenningTrap trap_with = PenningTrap(B_0, V_0, d, pl, f[1], w);
    std::vector<Particle> original_pl = pl;
    for(int j=0; j<int(500/dt); j++){
      //evelove the two traps using runge-kutta
      trap_without.evolve_RK4(dt, false);
      trap_with.evolve_RK4(dt, true);
    }
    //Count the number of particles still inside the trap at the end
    int count_part1 = count_part_inside(trap_without.particle_l, d);
    int count_part2 = count_part_inside(trap_with.particle_l, d);

    myfile2 << w << " " << count_part1 << " " << count_part2 << "\n";
    pl = original_pl;
  }
}
