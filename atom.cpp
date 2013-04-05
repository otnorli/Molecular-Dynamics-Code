#include "atom.h"

Atom::Atom(vec ri)
{
    r = ri;
}

vec3 Atom::getPosition()
{
    return r;
}

vec3 Atom::getVelocity()
{
    return v;
}

vec3 Atom::getForce()
{
    return f;
}

void Atom::setPosition(const vec &inPosition)
{
    r = inPosition;
}

void Atom::setVelocity(const vec &inVelocity)
{
    v = inVelocity;
}

void Atom::setForce(const vec &inForce)
{
    f += inForce;
}

void Atom::addForce(const vec &force)
{
    f += force;
}

double Atom::getPotential()
{
    return potential;
}

double Atom::getPressure()
{
    return pressure;
}

void Atom::setPotential(double inPotential)
{
    potential += inPotential;
}

double Atom::getKinetic()
{
    return kinetic;
}

void Atom::setPressure(double inPressure)
{
    pressure = inPressure;
}

vec3 Atom::getDisplacement()
{
    return d;
}

void Atom::setKinetic(double inKinetic)
{
    kinetic = inKinetic;
}

void Atom::setDisplacement(const vec &inDisplacement)
{
    d += inDisplacement;
}
