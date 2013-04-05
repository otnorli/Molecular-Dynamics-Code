#ifndef PRESSURE_CELLS_H
#define PRESSURE_CELLS_H

#include <iostream>
#include <string>
#include "atom.h"

using namespace std;

class Pressure_Cells
{
public:
    Pressure_Cells();
    void addAtoms(Atom *A);
    const vector<Atom*>& getAtoms();
    vector<Atom*> atoms;
    vec PcellIndices;
};

#endif // PRESSURE_CELLS_H
