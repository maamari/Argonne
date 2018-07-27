#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

#define PI 3.14159265
#define C 3.0e8
double momToKin(double massL, double mom);

//Converts mass and energy to momentum
double momentum(double mass, double kEnergy){
  return (sqrt(2*mass*kEnergy*(1+kEnergy/2/mass)));
}

//Gives the possible roots for momentum via quadratic calculaion
double roots(double a, double b, double c){
  double discriminant = pow(b,2) - 4*a*c;
  
  if(!discriminant)
    return -1;
  //else {
    double root1 = (-b - sqrt(discriminant))/2*a;
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
double collision(double Kb, double massB, double massT, double massH, double massL, double angle){
  //Convert masses to energy
  double energyB = massB + Kb;
  double energyT = massT;
  double energyH;
  double energyL;
  double momentumB = momentum(massB, 180);

  //Set up coefficients of the quadratic equation
  //constants
  double constants = pow(massH,2) - pow(massB,2) - pow(massT,2) - pow(massL,2) + 2.0*Kb*massL + 2.0*massL*massT; // double check formulas
  //A
  double coefficient1 = 2*(Kb/massL + massT/massL);
  //B
  double coefficient2 = 2*(Kb)*cos((angle*PI)/180);
  //Send over to quadratic function to find roots
  double mom = roots(coefficient1, coefficient2, constants);

  //Send momentum found to function to convert to KE
  double keL = momToKin(massL, mom);

  return keL;
}

//Converts momentum to kinetic energy
double momToKin(double massL, double mom){
  //cout << mom << endl;
  double kineticEnergy = sqrt(pow(massL,2) + pow(mom,2)) - massL;
  //cout << kineticEnergy << endl;

  return kineticEnergy;
}

int main(){
  double massB, massT, massH, massL, angle, Kb;
  ofstream output("output.txt");

  //Take in masses in amu and convert to MeV/c^2
  cout << "Please input the following in amu:" << endl << "Mass of beam: ";
  cin >> massB;
  massB *= 931.5;

  cout << "Mass of target: ";
  cin >> massT;
  massT *= 931.5;

  cout << "Mass of resultant: ";
  cin >> massH;
  massH *= 931.5;

  cout << "Mass of proton: ";
  cin >> massL;
  massL *= 931.5;

  //ask for beam kinetic (180MeV)
  cout << "Energy beam: ";
  cin >> Kb;

  cout << "Would you like to output results to screen? (y/n) ";
  char answer;
  cin >> answer;
 
  //Take in angle of discharge
  for(double i = 20; i < 61; i++){
    angle = i;
    output << angle << "\t";

    //Send to calculate the collision
    double kineticEnergy = collision(Kb, massB, massT, massH, massL, angle);
    if(answer == 'y')
      cout << "The kinetic energy of the proton is " <<  kineticEnergy << " at an angle of " << angle << " degrees." << endl;
    output << kineticEnergy << endl;
  }

  return 0;
}
