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
  arma::vec cross_prod = {particle_l[i].charge*B0_*particle_l[i].velocity[1], -particle_l[i].charge*B0_*particle_l[i].velocity[0], 0};
  //cout<< cross_prod;
  arma::vec ma_kalles_noe = {particle_l[i].position[0], particle_l[i].position[1], -2*particle_l[i].position[2]};

  arma::vec F_ex = particle_l[i].charge *(V0_ / pow(d_,2))*ma_kalles_noe  + cross_prod;
  //cout << F_ex << "\n";
  return F_ex;
}


// The total force on particle_i from the other particles
// MOST LIKELY COMPLETELY WRONG
arma::vec PenningTrap::total_force_particles(int i){
  double k_e =  1.38935333*pow(10,5);
  arma::vec F_part = {k_e, k_e, k_e};
  if(particle_l.size()==1){
    F_part = {0,0,0};

  }else{
    for(int j=0; j<particle_l.size(); j++){
      F_part += particle_l[j].charge * ((particle_l[i].position - particle_l[j].position) / pow(abs(particle_l[i].position - particle_l[j].position), 3));
    }
  }
  return F_part;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i){
  return total_force_external(i) + total_force_particles(i);
}



void PenningTrap::evolve_RK4(double dt){

  int n = particle_l.size();
  arma::mat a(n, 3);


  arma::mat k1x(n, 3), k2x(n, 3), k3x(n, 3), k4x(n, 3), k1v(n, 3), k2v(n, 3), k3v(n, 3), k4v(n, 3);

  std::vector<Particle> original_particle_l = particle_l;

  for(int i = 0; i < n; i++){
    a[i,0] = total_force(i)[0] / particle_l[i].mass;
    a[i,1] = total_force(i)[1] / particle_l[i].mass;
    a[i,2] = total_force(i)[2] / particle_l[i].mass;
  }

  //k1 update
  for(int j = 0; j < n; j++){
    for(int i=0; i<3; i++){
      k1v[j,i] = a[j,i] * dt;
      //cout << k1v[j,i];
      k1x[j,i] = particle_l[j].velocity[i] * dt;

      particle_l[j].velocity[i] += k1v[j,i]/2;
      particle_l[j].position[i] += k1x[j,i]/2;
    }
  }

  //k2 update
  for(int i = 0; i < n; i++){
    a[i,0] = total_force(i)[0] / particle_l[i].mass;
    a[i,1] = total_force(i)[1] / particle_l[i].mass;
    a[i,2] = total_force(i)[2] / particle_l[i].mass;
  }

  for(int j = 0; j < n; j++){
      for(int i=0; i<3; i++){
        k2v[j,i] = a[j,i] * dt;
        k2x[j,i] = particle_l[j].velocity[i] * dt;

        particle_l[j].velocity[i] += k2v[j,i]/2;
        particle_l[j].position[i] += k2x[j,i]/2;
      }
    }

    //k3 update
    for(int i = 0; i < n; i++){
      a[i,0] = total_force(i)[0] / particle_l[i].mass;
      a[i,1] = total_force(i)[1] / particle_l[i].mass;
      a[i,2] = total_force(i)[2] / particle_l[i].mass;
    }

    for(int j = 0; j < n; j++){
        for(int i=0; i<3; i++){
          k3v[j,i] = a[j,i] * dt;
          k3x[j,i] = particle_l[j].velocity[i] * dt;

          particle_l[j].velocity[i] += k3v[j,i];
          particle_l[j].position[i] += k3x[j,i];
        }
      }

      //k4 update
    for(int i = 0; i < n; i++){
      a[i,0] = total_force(i)[0] / particle_l[i].mass;
      a[i,1] = total_force(i)[1] / particle_l[i].mass;
      a[i,2] = total_force(i)[2] / particle_l[i].mass;
    }

    for(int j = 0; j < n; j++){
        for(int i=0; i<3; i++){
          k4v[j,i] = a[j,i] * dt;
          k4x[j,i] = particle_l[j].velocity[i] * dt;
        }
      }

      for(int j = 0; j < n; j++){
          for(int i=0; i < 3; i++){
            particle_l[j].velocity[i] = original_particle_l[j].velocity[i] + 1./6*(k1v[j,i]+2*k2v[j,i]+2*k3v[j,i]+k4v[j,i]);
            //cout << 1/6*(k1v[j,i]+2*k2v[j,i]+2*k3v[j,i]+k4v[j,i]);
            particle_l[j].position[i] = original_particle_l[j].position[i] + 1./6*(k1x[j,i]+2*k2x[j,i]+2*k3x[j,i]+k4x[j,i]);
          }
        }
}





void PenningTrap::evolve_forward_Euler(double dt){
  int n = particle_l.size();
  arma::mat a(n, 3);
  for(int i = 0; i < n; i++){
    //arma::vec tf = total_force(i);
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
}
