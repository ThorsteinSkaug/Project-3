#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <armadillo>

using namespace std;


// Definition of the constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, vector<Particle> pl_in)
{
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
  particle_l = pl_in;
}


void PenningTrap::add_particle(Particle p_in)
{
  particle_l.push_back(p_in);
}


arma::vec PenningTrap::external_E_field(arma::vec r){
  arma::vec E_field = {(V0_/pow(d_, 2))*r(0), (V0_/pow(d_, 2))*r(1), (V0_/pow(d_, 2))*2*r(2)};
  return E_field;
}


arma::vec PenningTrap::external_B_field(arma::vec r){
    arma::vec B_field = {0,0,B0_};
    return B_field;
}


// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
  arma::vec cross_prod = {particle_l[i].charge*B0_*particle_l[i].velocity[1], particle_l[i].charge*B0_*particle_l[i].velocity[0], 0};  // MOST LIKELY ERROR HERE

  arma::vec F_ex = particle_l[i].charge * external_E_field(particle_l[i].position) + cross_prod;  // MOST LIKELY ERROR HERE

  return F_ex;
}


// The total force on particle_i from the other particles
// MOST LIKELY COMPLETELY WRONG
arma::vec PenningTrap::total_force_particles(int i){
  double k_e =  1.38935333*pow(10,5);
  arma::vec F_part = {k_e, k_e, k_e};
  for(Particle p : particle_l){
    F_part += p.charge * ((particle_l[i].position - p.position) / pow(abs(particle_l[i].position - p.position), 3));
  }
  return F_part;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i){
  return total_force_external(i) + total_force_particles(i);
}



void PenningTrap::evolve_RK4(double dt, vector<Particle>& updated_particle_l){
  arma::mat a;
  int n = particle_l.size();
  for(int i = 0; i < n; i++){
    a[i,0] = total_force(i)[0] / particle_l[i].mass;
    a[i,1] = total_force(i)[1] / particle_l[i].mass;
    a[i,2] = total_force(i)[2] / particle_l[i].mass;
  }


  arma::vec k1x, k2x, k3x, k4x, k1v, k2v, k3v, k4v;

  for(int j = 0; j < n; j++){
    for(int i=0; i<3; i++){
      k1v[i] = a[j,i] * dt;
      k1x[i] = particle_l[j].velocity[i] * dt;
    }

    Particle original_particle = particle_l[j];

    particle_l[j].position += k1x/2;
    particle_l[j].velocity += k1v/2;
    k2v = total_force(j)/particle_l[j].mass * dt;
    k2x = (original_particle.velocity + k1v/2) * dt;

    particle_l[j] = original_particle;
    particle_l[j].position += k2x/2;
    particle_l[j].velocity += k2v/2;
    k3v = total_force(j)/particle_l[j].mass * dt;
    k3x = (original_particle.velocity + k2v/2) * dt;

    particle_l[j] = original_particle;
    particle_l[j].position += k3x;
    particle_l[j].velocity += k3v;
    k4v = total_force(j)/particle_l[j].mass * dt;
    k4x = (original_particle.velocity + k3v) * dt;


    Particle new_particle(original_particle.charge, original_particle.mass, original_particle.position + 1/6 * (k1x+2*k2x+2*k3x+k4x), original_particle.velocity + 1/6 * (k1v+2*k2v+2*k3v+k4v));
    updated_particle_l.push_back(new_particle);
  }

}


//void PenningTrap::evolve_RK4(double dt){
//  arma::mat a;
//  int n = p_.size();
//  for(int i = 0; i < n; i++){
//    a(i) = total_force(i) / p_(i).m;
//  }
//
//  double k1, k2, k3, k4;
//  for(int j = 0; j < n; j++){
//    for(int k=0;k<3;k++){
//      k1 = dt*a(j,k);
//      k2 = dt/2*a(j,k)
//    }
//  }
//}


void PenningTrap::evolve_forward_Euler(double dt){
  arma::mat a;
  int n = particle_l.size();
  for(int i = 0; i < n; i++){
    a[i,0] = total_force(i)[0] / particle_l[i].mass;
    a[i,1] = total_force(i)[1] / particle_l[i].mass;
    a[i,2] = total_force(i)[2] / particle_l[i].mass;
  }

  for(int j = 0; j < n; j++){
    for(int i=0; i<3; i++){
      particle_l[j].position[i] += particle_l[j].velocity[i] *dt;
      particle_l[j].velocity[i] += a[j, i] * dt;
    }
  }
  cout << "hei" <<particle_l[0].position << "\n";

}
