#ifndef __Particle_hpp__
#define __Particle_hpp__
#include <armadillo>

using namespace std;

class Particle
{

public:
  double charge, mass; //Storing the charge and mass of the particle
  arma::vec position, velocity; //Storing the position and velocity as vector values

  // Constructor
  Particle(int q, double m, arma::vec r, arma::vec v);

};

#endif
