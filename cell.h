#ifndef CELL_H
#define CELL_H

#include <vector>
#include <iostream>
#include <string>
#include "atom.h"
#include "../../Desktop/FYS4460/include/armadillo"

using namespace std;
using namespace arma;

class Cell
{
public:
    Cell();
    void addAtoms(Atom *A);
    void addNeighbor(Cell *cell);
    const vector<Atom*>& getAtoms();
    const vector<Cell*>& getNeighbor();
    vec cellIndices;
    bool hasCalculated;
    vector<Atom*> atoms;
    vector<Cell*> neighbor;

};

#endif // CELL_H
