#ifndef FLOWPROFILECELLS_H
#define FLOWPROFILECELLS_H
#include <vector>
#include <iostream>
#include <string>
#include "../../Desktop/FYS4460/include/armadillo"
#include "atom.h"
#include "lib.h"

using namespace std;
using namespace arma;

class FlowProfileCells
{
public:
    FlowProfileCells();
    void addAtoms(Atom *A);
    const vector<Atom*>& getAtoms();
    vec cell_Index;
    vector<Atom*> atoms;
};

#endif // FLOWPROFILECELLS_H
