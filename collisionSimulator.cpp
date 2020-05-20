/* Created as part of a summer 
   internship at Argonne National 
   Lab.
   
   The following simulates the
   collision of Flourine-18 with
   Deuterium.
   
   Author: Karime Maamari
   Date: 07/2018                */

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

#define PI 3.14159265
#define C 3.0*pow(10,8)

// Converts momentum to kinetic energy
double momToKin(double Ml, double mom){
  return (sqrt(pow(Ml,2) + pow(mom,2)) - Ml);
}

// Converts mass and energy to momentum
double momentum(double mass, double kEnergy){
  return (sqrt((2*mass*kEnergy)+(pow(kEnergy,2))));
}

// Gives the possible roots for momentum via quadratic calculation
double roots(double a, double b, double c){
  double discriminant = pow(b,2)-(4*a*c);

  if(!discriminant)
    return -1;
  double root1 = (-b+sqrt(discriminant))/(2*a);
  return root1;
}

/* Function to simulate the collision. Takes in the masses (in MeV/c^2) and
   converts them to their energy counterparts (MeV). Passes the masses, angle, and
   momentum to a quadratic function which returns the momentum of the proton. This
   momentum is then passed to a function to convert from momentum to kinetic
   energy of the proton. */
double collision(double Kb, double Mb, double Mt, double Mh, double Ml, double angle){
  // Convert masses to energy
  double energyB = Mb + Kb;
  double Pb = momentum(Mb, Kb);

  // Set up coefficients of the quadratic equation
    // Constants
  double massConst = pow(Mh,2)-pow(Mb,2)-pow(Mt,2)-pow(Ml,2)-(2*Mt*energyB);
  double energyConst = 4*(pow(energyB*Ml,2)+(2*energyB*Mt*(pow(Ml,2)))+pow(Mt*Ml,2));
  double constants = energyConst - pow(massConst,2);

    // A in the usual quadratic form (-B +/- sqrt(B^2 - 4AC)/2A)
  double quadCoeff1 = 4*(pow(energyB,2)+(2*energyB*Mt)+pow(Mt,2)-pow(Pb*cos(angle*PI/180),2));

    // B in the usual quadratic form (-B +/- sqrt(B^2 - 4AC)/2A)
  double quadCoeff2 = 4*massConst*Pb*cos((angle*PI)/180);

  // Send over to quadratic function to find roots
  double mom = roots(quadCoeff1, quadCoeff2, constants);

  // Send momentum found to function to convert to KE
  double keL = momToKin(Ml, mom);
  return keL;
}

int main(){
  double Mb, Mt, Mh, Mh2, Mh3, Ml, angle, Kb;
  ofstream output("output.txt");
  ofstream output2("output2.txt");
  ofstream output3("output3.txt");

  // Masses in amu (optional functionality to input masses) and convert to MeV/c^2
  // cout << "Please input the following in amu:" << endl << "Mass of beam: ";
  // cin >> Mb;
  Mb = 18.009380;
  Mb *= 931.49432;

  // cout << "Mass of target: ";
  // cin >> Mt;
  Mt = 2.0141017;
  Mt *= 931.49432;

   // cout << "Mass of resultant: ";
   // cin >> Mh;
  Mh = 18.998403*931.49432;
  Mh2 = 19.017*931.49432;
  Mh3 = 19.029*931.49432;

  // cout << "Mass of proton: ";
  // cin >> Ml;
  Ml = 1.0072764;
  Ml *= 931.49432;

  // Beam kinetic energy
  // cout << "Energy beam: ";
  // cin >> Kb;
  Kb = 180;

  // cout << "Would you like to output results to screen? (y/n) ";
  // char answer;
  // cin >> answer;

  // Take in angle of discharge
  for(double i = 0; i < 180; i++){
    angle = i;
    output << angle << "\t";
    output2 << angle << "\t";
    output3 << angle << "\t";

    // Send to calculate the collision
    double kineticEnergy1 = collision(Kb, Mb, Mt, Mh, Ml, angle);
    double kineticEnergy2 = collision(Kb, Mb, Mt, Mh2, Ml, angle);
    double kineticEnergy3 = collision(Kb, Mb, Mt, Mh3, Ml, angle);

    // if(answer == 'y'){
    //   // cout << "The kinetic energy of the proton is " <<  kineticEnergy << " MeV at an angle of " << angle << " degrees." << endl;
    //   cout << "**GROUND STATE**" << endl << angle << " " << kineticEnergy1 << endl;
    //   cout << endl << "**10 MeV Excitation**" << endl << angle << " " << kineticEnergy2 << endl;
    //   cout << endl << "**20 MeV Excitation**" << endl << angle << " " << kineticEnergy3 << endl;
    // }

    output << kineticEnergy1 << endl;
    output2 << kineticEnergy2 << endl;
    output3 << kineticEnergy3 << endl;
  }
  return 0;
}
