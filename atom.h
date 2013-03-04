#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <string>
#include "../../Desktop/FYS4460/include/armadillo"
using namespace std;
using namespace arma;

class Atom
{
public:
  //  Atom();
    Atom(vec r);
    vec3 getPosition();
    vec3 getVelocity();
    vec3 getForce();
    void setPosition(const vec &inPosition);
    void setVelocity(const vec &inVelocity);
    void setForce(const vec &inForce);
    void addForce(const vec &force);

    double getPotential();
    double getKinetic();
    double getPressure();
    void setPotential(double inPotential);
    void setKinetic(double inKinetic);
    void setPressure(double inPressure);

    vec3 getDisplacement();
    void setDisplacement(const vec &inDisplacement);



protected:
    vec3 r;
    vec3 v;
    vec3 f;
    vec3 d;

    double potential;
    double kinetic;
    double pressure;
};

#endif // ATOM_H
