#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <armadillo>

using namespace std;


// Definition of the constructor without time_dependence
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, vector<Particle> pl_in)
{
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
  particle_l = pl_in;
  t_ = 0;

  E_field = &PenningTrap::external_E_field; //The pointer
}

// Definition of the constructor with time_dependence
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, vector<Particle> pl_in, double f, double w_V)
{
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
  particle_l = pl_in;
  f_ = f;
  w_V_ = w_V;
  t_ = 0;
  E_field = &PenningTrap::external_E_field_time_dependent; //The pointer
}


void PenningTrap::add_particle(Particle p_in)
{
  particle_l.push_back(p_in); //Add the new particle to the end of the vector storing the particles
}


arma::vec PenningTrap::external_E_field_time_dependent(arma::vec r, double t){
  double new_V0 = V0_*(1+f_*cos(w_V_*t)); //The time dependent V0
  arma::vec E_field = {(new_V0/pow(d_, 2))*r(0), (new_V0/pow(d_, 2))*r(1), -2*(new_V0/pow(d_, 2))*r(2)}; //Calculate the external E field

  return E_field;
}

arma::vec PenningTrap::external_E_field(arma::vec r, double t){
  arma::vec E_field = {(V0_/pow(d_, 2))*r(0), (V0_/pow(d_, 2))*r(1), -2*(V0_/pow(d_, 2))*r(2)}; //Calculate the external E field

  return E_field;
}

arma::vec PenningTrap::external_B_field(arma::vec r){
    arma::vec B_field = {0,0,B0_}; //The external B field
    return B_field;
}


// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
  arma::vec cross_prod = {particle_l[i].charge*B0_*particle_l[i].velocity[1], -particle_l[i].charge*B0_*particle_l[i].velocity[0], 0}; //Calculate the cross product q*vec(r) x vec(B)

  arma::vec F_ex = particle_l[i].charge * (this->*E_field)(particle_l[i].position, t_)  + cross_prod; //Calculate the force from the external using the correct E_field function.
  return F_ex;
}


// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){
  double k_e =  1.38935333*pow(10,5); //Defining the constant
  arma::vec F_part = {0, 0, 0}; //Defining the end vector
  for(int j=0; j<particle_l.size(); j++){ //Loop through all particles in the system
    if(j != i){ //Make sure we not include particle at index i
      double abs_vec = 0; //Storing the length of the (position vector i - position vector j)
      for(int l=0; l<3; l++){
        abs_vec += pow(particle_l[i].position[l] - particle_l[j].position[l], 2); //Calculate the length of the position vector
      }
      abs_vec = sqrt(abs_vec);
      for(int k=0; k<3; k++){
        F_part[k] += k_e*(particle_l[j].charge * ((particle_l[i].position[k] - particle_l[j].position[k]) / pow(abs_vec, 3))); //Calculate the total force on particle_i from the other particles
      }
    }
  }
  //cout << F_part[0] << '\n';
  return F_part;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i, bool particle_interaction){
  if(sqrt(pow(particle_l[i].position[0], 2) + pow(particle_l[i].position[1], 2) + pow(particle_l[i].position[2], 2)) > d_){ //Check if the particle is outside of the trap
    if(particle_interaction){ //Check if we have particle interactions or not
      return total_force_particles(i); //If we have particle interactions we need to add the force from other particles
    }else{
      arma::vec null_vec(3, arma::fill::zeros); //If not we add 0 force to the particle
      return null_vec;
    }
  }
  else if(particle_interaction){
    return total_force_external(i) + total_force_particles(i); //If particle interaction and inside trap, add both force from trap and from other particles
  }else{
    return total_force_external(i); //Else add only force from trap
  }
}



void PenningTrap::evolve_RK4(double dt, bool particle_interaction, bool time_dependence = false){

  int n = particle_l.size(); //Storing the number of particles
  arma::mat a(n, 3); //Storing the acceleration to particle_i as the i'th row

  arma::mat k1x(n, 3), k2x(n, 3), k3x(n, 3), k4x(n, 3), k1v(n, 3), k2v(n, 3), k3v(n, 3), k4v(n, 3); //Storing the kx and kv values

  std::vector<Particle> original_particle_l = particle_l; //Store the original particles before we take one step

  for(int i = 0; i < n; i++){ //Find the acceleration at the first time step
    a(i,0) = total_force(i, particle_interaction)[0] / particle_l[i].mass;
    a(i,1) = total_force(i, particle_interaction)[1] / particle_l[i].mass;
    a(i,2) = total_force(i, particle_interaction)[2] / particle_l[i].mass;
  }

  for(int j = 0; j < n; j++){ //Calculate the k1x and k1v values to all particles
    for(int i=0; i<3; i++){
      k1v(j,i) = a(j,i) * dt;
      //cout << k1v[j,i];
      k1x(j,i) = particle_l[j].velocity[i] * dt;

      //Update all particles such that we are in time step t+dt/2
      particle_l[j].velocity[i] += k1v(j,i)/2;
      particle_l[j].position[i] += k1x(j,i)/2;
    }
  }

  t_ += dt/2; //Update the global time

  for(int i = 0; i < n; i++){ //Find the acceleration at the second time step
    a(i,0) = total_force(i, particle_interaction)[0] / particle_l[i].mass;
    a(i,1) = total_force(i, particle_interaction)[1] / particle_l[i].mass;
    a(i,2) = total_force(i, particle_interaction)[2] / particle_l[i].mass;
  }

  for(int j = 0; j < n; j++){ //Calculate the k2x and k2v values to all particles
      for(int i=0; i<3; i++){
        k2v(j,i) = a(j,i) * dt;
        k2x(j,i) = particle_l[j].velocity[i] * dt;

        particle_l[j].velocity[i] += k2v(j,i)/2;
        particle_l[j].position[i] += k2x(j,i)/2;
      }
    }

    for(int i = 0; i < n; i++){ //Find the acceleration at the third time step
      a(i,0) = total_force(i, particle_interaction)[0] / particle_l[i].mass;
      a(i,1) = total_force(i, particle_interaction)[1] / particle_l[i].mass;
      a(i,2) = total_force(i, particle_interaction)[2] / particle_l[i].mass;
    }

    for(int j = 0; j < n; j++){ //Calculate the k3x and k3v values to all particles
        for(int i=0; i<3; i++){
          k3v(j,i) = a(j,i) * dt;
          k3x(j,i) = particle_l[j].velocity[i] * dt;

          particle_l[j].velocity[i] += k3v(j,i);
          particle_l[j].position[i] += k3x(j,i);
        }
      }

    t_ += dt/2; //Update the global t value

    for(int i = 0; i < n; i++){ //Find the acceleration at the last time step
      a(i,0) = total_force(i, particle_interaction)[0] / particle_l[i].mass;
      a(i,1) = total_force(i, particle_interaction)[1] / particle_l[i].mass;
      a(i,2) = total_force(i, particle_interaction)[2] / particle_l[i].mass;
    }

    for(int j = 0; j < n; j++){ //Calculate the k4x and k4v values to all particles
        for(int i=0; i<3; i++){
          k4v(j,i) = a(j,i) * dt;
          k4x(j,i) = particle_l[j].velocity[i] * dt;
        }
      }

      for(int j = 0; j < n; j++){ //Finally update all the original positions and velocities with the calculated kx and kv values
          for(int i=0; i < 3; i++){
            particle_l[j].velocity[i] = original_particle_l[j].velocity[i] + 1./6*(k1v(j,i)+2*k2v(j,i)+2*k3v(j,i)+k4v(j,i));
            particle_l[j].position[i] = original_particle_l[j].position[i] + 1./6*(k1x(j,i)+2*k2x(j,i)+2*k3x(j,i)+k4x(j,i));
          }
        }
}



void PenningTrap::evolve_forward_Euler(double dt, bool particle_interaction){
  int n = particle_l.size(); //Storing the number of particles
  arma::mat a(n, 3); //Storing the acceleration to particle_i as the i'th row
  for(int i = 0; i < n; i++){
    a[i,0] = total_force(i, particle_interaction)[0] / particle_l[i].mass;
    a[i,1] = total_force(i, particle_interaction)[1] / particle_l[i].mass;
    a[i,2] = total_force(i, particle_interaction)[2] / particle_l[i].mass;
  }

  for(int j = 0; j < n; j++){
    for(int i=0; i<3; i++){
      particle_l[j].position[i] += particle_l[j].velocity[i] *dt; //Update the position
      particle_l[j].velocity[i] += a[j, i] * dt; //Update the velocity
    }
  }
}
