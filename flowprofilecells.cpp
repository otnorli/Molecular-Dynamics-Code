#include "flowprofilecells.h"

FlowProfileCells::FlowProfileCells()
{
}

void FlowProfileCells::addAtoms(Atom *A)
{
     atoms.push_back(A);
}

const vector<Atom *> &FlowProfileCells::getAtoms()
{
     return atoms;
}
