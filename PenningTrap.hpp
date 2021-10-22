#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__
#include "Particle.hpp"
#include <armadillo>

using namespace std;


class PenningTrap
{

public:
  double B0_, V0_, d_, f_, w_V_; //Storing the magnetic field strength. applied potential and dimension of the penning trap
  std::vector<Particle> particle_l; //Storing a vector with the variables from the particle class.

  // Constructor
  PenningTrap(double B0_in, double V0_in, double d_in, vector<Particle> pl, double f, double w_V);

  // Add a particle to the trap
  void add_particle(Particle p_in);

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r, double t);

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r);

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j);

  // The total force on particle_i from the external fields
  arma::vec total_force_external(int i);

  // The total force on particle_i from the other particles
  arma::vec total_force_particles(int i);

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force(int i, bool particle_interaction);

  // Evolve the system one time step (dt) using Runge-Kutta 4th order
  void evolve_RK4(double dt, bool particle_interaction, bool time_dependence);

  // Evolve the system one time step (dt) using Forward Euler
  void evolve_forward_Euler(double dt, bool particle_interaction);
};


#endif
