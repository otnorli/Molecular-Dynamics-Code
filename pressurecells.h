#ifndef PRESSURECELLS_H
#define PRESSURECELLS_H

#include <vector>
#include <iostream>
#include <string>
#include "atom.h"
#include "../../Desktop/FYS4460/include/armadillo"

using namespace std;
using namespace arma;

class PressureCells
{
public:
    PressureCells();
    void addAtoms(Atom *A);
    const vector<Atom*>& getAtoms();
    vec cell_index;
    vector<Atom*> atoms;
};

#endif // PRESSURECELLS_H
