#ifndef GENERATEQUANTITIES_H
#define GENERATEQUANTITIES_H

#include <iostream>
#include <string>
#include "atom.h"
#include "cell.h"
#include "lib.h"
#include "potentials.h"
#include "../../Desktop/FYS4460/include/armadillo"



using namespace std;
using namespace arma;

class GenerateQuantities
{
public:
    //Krefter og hastigheter
    GenerateQuantities(string Thermostat);
    double NormalDist(long *idum);
    void PrintToFile();
    void printVelocity();
    void generatePosition();
    void generateVelocity();
    void generateForce();
    void generateCells();
    void acceleration(Atom *atom1, Atom *atom2);
    void ClearCells();

    //Fysiske størrelser
    double radialDF(Atom *atom1, Atom *atom2);

    //Termostater
    double BerendsenThermostat(double Temperature);
    void AndersenThermostat(Atom *atom);



protected:
    //Krefter og posisjoner og hastigheter
    Cell* cells;
    vector<Atom*> atoms;
    Potentials* potential;
    double temporary;
    int time_steps;
    int N, N_tot, Lc, nCells;
    int i,j,k,n,p,t;
    double Temp;
    double sigma, tau, k_b, T0, T_bath, b;
    double m, eV, epsilon, standev, L, dt, d;
    double cellSize, r_cut;
    long idum;
    int BinIndex;
    string atom_name;
    double dr;
    int termostat_limit;
    vec r1;
    vec r2;
    vec drvec;
    vec force;

    int B;

    double density, volume;


    //Termostat
    string thermostat;
    double relaxation_time;
    double andersen_std;
};

#endif // GENERATEQUANTITIES_H
