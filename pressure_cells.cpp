#include "pressure_cells.h"

Pressure_Cells::Pressure_Cells()
{
}

void Pressure_Cells::addAtoms(Atom *A)
{
    atoms.push_back(A);
}

const vector<Atom *> &Pressure_Cells::getAtoms()
{
    return atoms;
}
