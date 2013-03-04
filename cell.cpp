#include "cell.h"

Cell::Cell()
{
}

void Cell::addAtoms(Atom *A)
{
    atoms.push_back(A);
}

void Cell::addNeighbor(Cell *cell)
{
    neighbor.push_back(cell);
}

const vector<Atom *> &Cell::getAtoms()
{
    return atoms;
}

const vector<Cell *> &Cell::getNeighbor()
{
    return neighbor;
}
