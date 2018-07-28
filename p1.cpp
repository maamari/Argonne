#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

#define PI 3.14159265
#define C 3.0*pow(10,8)

//Converts momentum to kinetic energy
double momToKin(double Ml, double mom){
  return (sqrt(pow(Ml,2) + pow(mom,2)) - Ml);
}

//Converts mass and energy to momentum
double momentum(double mass, double kEnergy){
  return (sqrt((2*mass*kEnergy)+(pow(kEnergy,2))));
}

//Gives the possible roots for momentum via quadratic calculaion
double roots(double a, double b, double c){
  double discriminant = pow(b,2)-(4*a*c);

  if(!discriminant)
    return -1;
  //else {
    double root1 = (-b+sqrt(discriminant))/(2*a);
    return root1;
    //cout << "Root 1 is " << root1 << endl;

    //double root2 = (-b + sqrt(discriminant))/2*a;
    //cout << "Root 2 is " << root2 << endl;

    //int choice;
    //cout << "Which root would you like, 1 or 2? ";
    //cin >> choice;

    //if(choice == 1)
      //return root1;
    //else if(choice == 2)
      //return root2;
  //}
  //cout << "I'm sorry that was an invalid entry." << endl;
  //return -1;
}

/* Function to simulate the collision. Takes in the masses (in MeV/c^2) and
   converts them to their energy counterparts (MeV). Passes the masses, angle, and
   momentum to a quadratic function which returns the momentum of the proton. This
   momentum is then passed to a function to convert from momentum to kinetic
   energy of the proton. */
double collision(double Kb, double Mb, double Mt, double Mh, double Ml, double angle){
  //Convert masses to energy
  double energyB = Mb + Kb;
  double Pb = momentum(Mb, Kb);
  // cout << Pb << endl;

  //Set up coefficients of the quadratic equation
    //constants
  double a0 = pow(Mh,2)-pow(Mb,2)-pow(Mt,2)-pow(Ml,2)-(2*Mt*energyB);
  // double a1 = 4*(pow(energyB*Ml,2)+(2*energyB*Mt*(pow(Ml,2)))+pow(Mt*Ml,2));
  double constants = 4*pow(energyB+Mt,2)-pow(a0,2);
    //A
  double coefficient1 = 4*pow(energyB+Mt,2)-pow(Pb*cos(angle*PI/180),2);
    //B
  double coefficient2 = 4*a0*Pb*cos((angle*PI)/180);

  //Send over to quadratic function to find roots
  double mom = roots(coefficient1, coefficient2, constants);

  // cout << mom << endl;
  //Send momentum found to function to convert to KE
  double keL = momToKin(Ml, mom);
  return keL;
}

int main(){
  double Mb, Mt, Mh, Ml, angle, Kb;
  ofstream output("output.txt");

  //Take in masses in amu and convert to MeV/c^2
  // cout << "Please input the following in amu:" << endl << "Mass of beam: ";
  // cin >> Mb;
  Mb = 18;
  Mb *= 931.49432;

  // cout << "Mass of target: ";
  // cin >> Mt;
  Mt = 2;
  Mt *= 931.49432;

  // cout << "Mass of resultant: ";
  // cin >> Mh;
  Mh = 19;
  Mh *= 931.49432;

  // cout << "Mass of proton: ";
  // cin >> Ml;
  Ml = 1;
  Ml *= 931.49432;

  //ask for beam kinetic (180MeV)
  // cout << "Energy beam: ";
  Kb = 180;
  // cin >> Kb;

  cout << "Would you like to output results to screen? (y/n) ";
  char answer;
  cin >> answer;

  //Take in angle of discharge
  for(double i = 0; i < 180; i++){
    angle = i;
    output << angle << "\t";

    //Send to calculate the collision
    double kineticEnergy = collision(Kb, Mb, Mt, Mh, Ml, angle);
    if(answer == 'y')
      // cout << "The kinetic energy of the proton is " <<  kineticEnergy << " MeV at an angle of " << angle << " degrees." << endl;
      cout << angle << " " << kineticEnergy << endl;
    output << kineticEnergy << endl;
  }

  return 0;
}
