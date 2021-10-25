#include "Particle.hpp"
#include <armadillo>
using namespace std;

// Definition of the constructor
Particle::Particle(int q, double m, arma::vec r, arma::vec v)
{
  // Use the input variables (q, m, r, v) to assign values to the class memeber variables
  charge = q;
  mass = m;
  position = r;
  velocity = v;
}
