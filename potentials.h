#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <iostream>
#include <string>
#include "atom.h"
#include "lib.h"
#include "../../Desktop/FYS4460/include/armadillo"

using namespace std;
using namespace arma;

class Potentials
{
public:
    Potentials();
    vec Zero_potential(const vec3 &dr);
    vec Lennard_Jones_potential(const vec3 &dr);
protected:
    vec3 f;
    vec3 zeroVector;
    double innerProduct;
};

#endif // POTENTIALS_H
