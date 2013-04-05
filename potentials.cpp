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
    double tempProduct = 1/(innerProduct*innerProduct*innerProduct);
    double tempProduct2 = tempProduct / innerProduct;
    f = 24.0 * (2.0 * tempProduct - 1.0 ) * dr * tempProduct2;
    return f;
}
