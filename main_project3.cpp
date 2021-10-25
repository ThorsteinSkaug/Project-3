#include "PenningTrap.cpp"
#include "Particle.cpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//Function for calculating the analytical solution for one particle
void analytical_solution(arma::vec& x, arma::vec& y, arma::vec& z, double dt, double w_z, double w_p, double w_m, double A_p, double A_m){
  for(int i=1; i<x.size(); i++){
    x(i) = A_p*cos(w_p*i*dt) + A_m*cos(w_m*i*dt); //Calculate the x-coordinates at time step (dt*i)
    y(i) = -A_p*sin(w_p*i*dt) - A_m*sin(w_m*i*dt); //Calculate the y-coordinates at time step (dt*i)
    z(i) = z(0) * cos(w_z*(i*dt)); //Calculate the z-coordinates at time step (dt*i)
  }
}


int main(){
    int t = 100; //The end time

    //Defining some constants
    double T = 9.64852558 * 10;
    double V = 9.64852558 * pow(10,7);
    double B_0 = 1. * T;
    double V_0 = 10 * V;
    double d = pow(10,4);

    //One particle
    int q = 1.; //The charge of the particle
    double m = 40.078; //The mass of the particle
    arma::vec r = {1., 0., 1.}; //The position vector
    arma::vec v = {0., 1., 0.}; //The velocity vector
    vector<Particle> pl; //Vector for storing all the particles in the trap
    Particle singly_charged_Calcium = Particle(q, m, r, v); //Make a particle
    PenningTrap trap = PenningTrap(B_0, V_0, d, pl); //Make a trap
    trap.add_particle(singly_charged_Calcium); //Add the particle to the trap
    double dt = 0.001; //The length of each time step


    //Start analytical solution part
    double w_0 = singly_charged_Calcium.charge*B_0/singly_charged_Calcium.mass; //Constant needed for calculating the analytical solution
    double w_z = sqrt(2*singly_charged_Calcium.charge*V_0/(singly_charged_Calcium.mass*pow(d,2))); //Constant needed for calculating the analytical solution

    double w_p = (w_0 + sqrt(pow(w_0,2)-2*pow(w_z,2)))/2; //Constant needed for calculating the analytical solution
    double w_m = (w_0 - sqrt(pow(w_0,2)-2*pow(w_z,2)))/2; //Constant needed for calculating the analytical solution

    bool particle_interaction1particle = false; //Set the partical interactions of when we have only one particle


    double A_p = (singly_charged_Calcium.velocity[1]+w_m*singly_charged_Calcium.position[0])/(w_m - w_p); //Constant needed for calculating the analytical solution
    double A_m = -(singly_charged_Calcium.velocity[1]+w_p*singly_charged_Calcium.position[0])/(w_m - w_p); //Constant needed for calculating the analytical solution

    arma::vec x(int(t/dt)+1); //A vector for storing the analytical solution
    x(0) = trap.particle_l[0].position[0]; //Add the starting point
    arma::vec y(int(t/dt)+1); //A vector for storing the analytical solution
    y(0) = trap.particle_l[0].position[1]; //Add the starting point
    arma::vec z(int(t/dt)+1); //A vector for storing the analytical solution
    z(0) = trap.particle_l[0].position[2]; //Add the starting point
    analytical_solution(x, y, z, dt, w_z, w_p, w_m, A_p, A_m); //Find the analytical solution

    std::ofstream myfile_an;
    myfile_an.open("coordinates_an.txt"); //Write analytical solution to file
    for(int i=0; i<x.size(); i++){
      myfile_an << std::scientific << dt*(i) << " " << std::scientific << x[i] << " " << std::scientific << y[i] <<" " << std::scientific << z[i] << "\n";
    }
    //End analytical solution part


    //Start RK4 solution part single particle
    std::ofstream myfile;
    myfile.open("coordinates_rk.txt"); //Write runge-kutta solution to file
    myfile << std::scientific << 0 << " " << std::scientific << trap.particle_l[0].position[0] << " " << std::scientific << trap.particle_l[0].position[1] <<" " << std::scientific << trap.particle_l[0].position[2] << "\n";
    for(int i=0; i<int(t/dt); i++){
      trap.evolve_RK4(dt, particle_interaction1particle); //Evolve the particle using runge-kutta
      myfile << std::scientific << dt*(i+1) << " " << std::scientific << trap.particle_l[0].position[0] << " " << std::scientific << trap.particle_l[0].position[1] <<" " << std::scientific << trap.particle_l[0].position[2] << "\n"; //Write the number of iteration needed for convergence to file
    }

    std::vector<double> h = {0.1, 0.05, 0.01, 0.005, 0.001}; //Different step size values to test
    for(int i=0; i<5; i++){
      singly_charged_Calcium = Particle(q, m, r, v); //Make a particle
      pl.clear(); //Clear the particle vector such that we can start a new time step at time 0
      trap = PenningTrap(B_0, V_0, d, pl); //Make a trap
      trap.add_particle(singly_charged_Calcium); //Add the particle to the trap

      //Finding the analytical solution
      arma::vec x(int(t/h[i])+1);
      x(0) = trap.particle_l[0].position[0];
      arma::vec y(int(t/h[i])+1);
      y(0) = trap.particle_l[0].position[1];
      arma::vec z(int(t/h[i])+1);
      z(0) = trap.particle_l[0].position[2];
      analytical_solution(x, y, z, h[i], w_z, w_p, w_m, A_p, A_m);

      std::ofstream myfile;
      dt = h[i]; //The given time step

      //std::vector<double> r_vec = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
      //cout <<  abs(sqrt(pow(trap.particle_l[0].position[0],2)+pow(trap.particle_l[0].position[1],2)+pow(trap.particle_l[0].position[1],2))-(sqrt(pow(x,2)+pow(y,2)+pow(z,2))))/(sqrt(pow(trap.particle_l[0].position[0],2)+pow(trap.particle_l[0].position[1],2)+pow(trap.particle_l[0].position[1],2))) <<'\n';
      // Runge Kutta error
      myfile.open("coordinates_rk" + std::to_string(dt) + ".txt");
      myfile << std::scientific << 0 << " " << 0 << " "<< 0 <<"\n";
      for(int i=0; i<int(t/dt); i++){

        trap.evolve_RK4(dt, particle_interaction1particle); //Evolve the particle using runge-kutta
        double error =  sqrt(pow(trap.particle_l[0].position[0]-x(i+1),2)+pow(trap.particle_l[0].position[1]-y(i+1),2)+pow(trap.particle_l[0].position[2]-z(i+1),2)) / sqrt(pow(x(i+1),2)+pow(y(i+1),2)+pow(z(i+1),2));

        myfile << std::scientific << dt*(i+1) << " " << std::scientific << error << " " << sqrt(pow(trap.particle_l[0].position[0]-x(i+1),2)+pow(trap.particle_l[0].position[1]-y(i+1),2)+pow(trap.particle_l[0].position[2]-z(i+1),2)) << '\n';
      }
    }
    //End RK4 solution part

    //Start euler solution part
    vector<Particle> pl2; //A new vector for storing particles
    Particle singly_charged_Calcium2 = Particle(q, m, r, v); //A new singly charged calcium ion
    PenningTrap trap2 = PenningTrap(B_0, V_0, d, pl2); //A new Penning trap
    trap2.add_particle(singly_charged_Calcium2); //Add the particle

    std::ofstream myfile2;
    myfile2.open("coordinates_euler.txt");
    myfile2 << std::scientific << 0 << " " << std::scientific << trap2.particle_l[0].position[0] << " " << std::scientific << trap2.particle_l[0].position[1] <<" " << std::scientific << trap2.particle_l[0].position[2] << "\n";
    dt = 0.001;
    for(int i=0; i<int(t/dt); i++){
      trap2.evolve_forward_Euler(dt, particle_interaction1particle); //Evolve the particle using forward euler
      myfile2 << std::scientific << dt*(i+1) << " " << std::scientific << trap2.particle_l[0].position[0] << " " << std::scientific << trap2.particle_l[0].position[1] <<" " << std::scientific << trap2.particle_l[0].position[2] << "\n"; //Write the number of iteration needed for convergence to file
    }

    // Forward Euler error
    for(int i=0; i<5; i++){
      singly_charged_Calcium = Particle(q, m, r, v); //Make a particle
      pl.clear(); //Clear the particle vector
      trap2 = PenningTrap(B_0, V_0, d, pl); //Make a trap
      trap2.add_particle(singly_charged_Calcium); //Add the particle to the trap

      std::ofstream myfile3;
      dt = h[i];

      //Find analytical solution
      arma::vec x(int(t/h[i])+1);
      x(0) = trap2.particle_l[0].position[0];
      arma::vec y(int(t/h[i])+1);
      y(0) = trap2.particle_l[0].position[1];
      arma::vec z(int(t/h[i])+1);
      z(0) = trap2.particle_l[0].position[2];
      analytical_solution(x, y, z, h[i], w_z, w_p, w_m, A_p, A_m);

      myfile3.open("coordinates_eu" + std::to_string(dt) + ".txt");
      myfile3 << std::scientific << 0 << " " << 0 << " "<< 0 << "\n";
      for(int i=0; i<int(t/dt); i++){

        trap2.evolve_forward_Euler(dt, particle_interaction1particle); //Evolve the particle using forward euler
        double error =  sqrt(pow(trap2.particle_l[0].position[0]-x(i+1),2)+pow(trap2.particle_l[0].position[1]-y(i+1),2)+pow(trap2.particle_l[0].position[2]-z(i+1),2)) / sqrt(pow(x(i+1),2)+pow(y(i+1),2)+pow(z(i+1),2));

        myfile3 << std::scientific << dt*(i+1) << " " << std::scientific << error << " " << sqrt(pow(trap2.particle_l[0].position[0]-x(i+1),2)+pow(trap2.particle_l[0].position[1]-y(i+1),2)+pow(trap2.particle_l[0].position[2]-z(i+1),2)) << '\n';

      }
    }
    //End euler solution part


    //Two particles:
    bool particle_interaction2particle_without = false; //Start by setting particle_interaction to false

    arma::vec r2 = {0.5, 2.3, 0.2}; //A new position vector
    arma::vec v2 = {-0.3, 1., 2.2}; //A new velocity vector
    std::vector<Particle> pl2particles; //A new vector for storing particles
    Particle singly_charged_Calcium_nr2 = Particle(q, m, r2, v2); //A new particle
    singly_charged_Calcium = Particle(q, m, r, v); //The same particle as in the one particle case
    PenningTrap trap2particles = PenningTrap(B_0, V_0, d, pl2particles); //A new trap
    trap2particles.add_particle(singly_charged_Calcium); //Add the first particle
    trap2particles.add_particle(singly_charged_Calcium_nr2); //Add the new particle

    std::ofstream myfile2particles;
    myfile2particles.open("rk2particles_without.txt");
    myfile2particles << std::scientific << 0 << " " << std::scientific << trap2particles.particle_l[0].position[0] << " " << std::scientific << trap2particles.particle_l[0].position[1] <<" " << std::scientific << trap2particles.particle_l[0].position[2] << " " << std::scientific << trap2particles.particle_l[0].velocity[0] << " " << std::scientific << trap2particles.particle_l[0].velocity[1] <<" " << std::scientific << trap2particles.particle_l[0].velocity[2] << " " << std::scientific << trap2particles.particle_l[1].position[0] << " " << std::scientific << trap2particles.particle_l[1].position[1] <<" " << std::scientific << trap2particles.particle_l[1].position[2] << " " << std::scientific << trap2particles.particle_l[1].velocity[0] << " " << std::scientific << trap2particles.particle_l[1].velocity[1] <<" " << std::scientific << trap2particles.particle_l[1].velocity[2] <<"\n";
    for(int i=0; i<int(t/dt); i++){
      trap2particles.evolve_RK4(dt, particle_interaction2particle_without); //Evolve both particles in the trap using runge-kutta without interactions
      myfile2particles << std::scientific << (dt*(i+1)) << " " << std::scientific << trap2particles.particle_l[0].position[0] << " " << std::scientific << trap2particles.particle_l[0].position[1] <<" " << std::scientific << trap2particles.particle_l[0].position[2] << " " << std::scientific << trap2particles.particle_l[0].velocity[0] << " " << std::scientific << trap2particles.particle_l[0].velocity[1] <<" " << std::scientific << trap2particles.particle_l[0].velocity[2] << " " << std::scientific << trap2particles.particle_l[1].position[0] << " " << std::scientific << trap2particles.particle_l[1].position[1] <<" " << std::scientific << trap2particles.particle_l[1].position[2] << " " << std::scientific << trap2particles.particle_l[1].velocity[0] << " " << std::scientific << trap2particles.particle_l[1].velocity[1] <<" " << std::scientific << trap2particles.particle_l[1].velocity[2] <<"\n";
    }

    //Two particles:
    bool particle_interaction2particle_with = true; //Set the intaractions to true

    //Make a new trap with equal initial condition as the previous
    std::vector<Particle> pl2particles2;
    singly_charged_Calcium_nr2 = Particle(q, m, r2, v2);
    singly_charged_Calcium = Particle(q, m, r, v);
    PenningTrap trap2particles2 = PenningTrap(B_0, V_0, d, pl2particles2);
    trap2particles2.add_particle(singly_charged_Calcium);
    trap2particles2.add_particle(singly_charged_Calcium_nr2);

    std::ofstream myfile2particles2;
    myfile2particles2.open("rk2particles_with.txt");
    myfile2particles2 << std::scientific << 0 << " " << std::scientific << trap2particles2.particle_l[0].position[0] << " " << std::scientific << trap2particles2.particle_l[0].position[1] <<" " << std::scientific << trap2particles2.particle_l[0].position[2] << " " << std::scientific << trap2particles2.particle_l[0].velocity[0] << " " << std::scientific << trap2particles2.particle_l[0].velocity[1] <<" " << std::scientific << trap2particles2.particle_l[0].velocity[2] << " " << std::scientific << trap2particles2.particle_l[1].position[0] << " " << std::scientific << trap2particles2.particle_l[1].position[1] <<" " << std::scientific << trap2particles2.particle_l[1].position[2] << " " << std::scientific << trap2particles2.particle_l[1].velocity[0] << " " << std::scientific << trap2particles2.particle_l[1].velocity[1] <<" " << std::scientific << trap2particles2.particle_l[1].velocity[2] <<"\n";
    for(int i=0; i<int(t/dt); i++){
      trap2particles2.evolve_RK4(dt, particle_interaction2particle_with); //Evolve both particles in the trap using runge-kutta with interactions
      myfile2particles2 << std::scientific << dt*(i+1) << " " << std::scientific << trap2particles2.particle_l[0].position[0] << " " << std::scientific << trap2particles2.particle_l[0].position[1] <<" " << std::scientific << trap2particles2.particle_l[0].position[2] << " " << std::scientific << trap2particles2.particle_l[0].velocity[0] << " " << std::scientific << trap2particles2.particle_l[0].velocity[1] <<" " << std::scientific << trap2particles2.particle_l[0].velocity[2] << " " << std::scientific << trap2particles2.particle_l[1].position[0] << " " << std::scientific << trap2particles2.particle_l[1].position[1] <<" " << std::scientific << trap2particles2.particle_l[1].position[2] << " " << std::scientific << trap2particles2.particle_l[1].velocity[0] << " " << std::scientific << trap2particles2.particle_l[1].velocity[1] <<" " << std::scientific << trap2particles2.particle_l[1].velocity[2] <<"\n";
    }


    return 0;
}
