#include "potentials.h"

Potentials::Potentials()
{
    zeroVector = zeros(3);
}

vec Potentials::Zero_potential(const vec3 &dr)
{
    return zeroVector;
}

vec Potentials::Lennard_Jones_potential(const vec3 &dr)
{
    innerProduct = dot(dr,dr);
    f = 24.0 * (2.0 * pow(innerProduct,-3.0) - 1.0 ) * dr * pow(innerProduct,-4.0);
    return f;
}
