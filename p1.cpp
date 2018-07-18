#include <iostream>
#include <cmath>
using namespace std;

#define PI 3.14159265

//Converts mass and energy to momentum
double momentum(double mass, double energy){
  return (pow(energy,2)-pow(mass,2));
}

//Gives the possible roots for momentum via quadratic calculaion
double roots(double a, double b, double c){
  double discriminant = pow(b,2) + 4*a*c;

  if(!discriminant)
    return -1;
  else {
    double root1 = (-b - sqrt(discriminant))/2*a;
    cout << "Root 1 is " << root1 << endl;

    double root2 = (-b + sqrt(discriminant))/2*a;
    cout << "Root 2 is " << root2 << endl;

    int choice;
    cout << "Which root would you like, 1 or 2? ";
    cin >> choice;

    if(choice == 1)
      return root1;
    else if(choice == 2)
      return root2;
  }
  cout << "I'm sorry that was an invalid entry." << endl;
  return -1;
}

/* Function to simulate the collision. Takes in the masses (in MeV/c^2) and
converts them to their energy counterparts (MeV). Passes the masses, angle, and
momentum to a quadratic function which returns the momentum of the proton. This
momentum is then passed to a function to convert from momentum to kinetic
energy of the proton. */
double collision(double massB, double massT, double massH, double massL, double angle){
  double energyB = massB*(3.0*pow(10,8));
  double energyT = massT*(3.0*pow(10,8));
  double energyH = massH*(3.0*pow(10,8));
  double energyL = massL*(3.0*pow(10,8));
  double momentumB = momentum(massB, energyB);

  double constants = pow(massH,2) - pow(massB,2) - pow(massT,2) - pow(massL,2) + 2.0*energyB*massL + 2.0*massL*massT;
  double coefficient1 = 2*(energyB/massL + massT/massL);
  double coefficient2 = 2*(momentumB)*cos((angle*PI)/180);
  double momentumL = roots(coefficient1, coefficient2, constants);
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
  double massB, massT, massH, massL, angle;

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

  cout << endl << "Please input the angle of the proton discharged in degrees: ";
  cin >> angle;

  double mom = collision(massB, massT, massH, massL, angle);

  return 0;
}
