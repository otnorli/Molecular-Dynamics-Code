#include "pressurecells.h"

PressureCells::PressureCells()
{
}

void PressureCells::addAtoms(Atom *A)
{
    atoms.push_back(A);
}


const vector<Atom *> &PressureCells::getAtoms()
{
    return atoms;
}
