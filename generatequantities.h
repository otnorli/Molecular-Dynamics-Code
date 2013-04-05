#ifndef GENERATEQUANTITIES_H
#define GENERATEQUANTITIES_H

#include <iostream>
#include <string>
#include "atom.h"
#include "cell.h"
#include "lib.h"
#include "potentials.h"
#include "../../Desktop/FYS4460/include/armadillo"
#include <cmath>
#include "pressurecells.h"
#include "flowprofilecells.h"

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
    void printPressure(const vec &pressss);
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

    //Flow
    vec3 FlowF;
    bool flow_on_off;
    bool Sylindric_or_circle;
    int Permeability_Counter_pluss, Permeability_Counter_minus;
    int Patrykt_Kraft_Retning;



protected:
    //Krefter og posisjoner og hastigheter
    Cell* cells;
    PressureCells* pressCells;
    vector<Atom*> atoms;
    Potentials* potential;
    double temporary;
    double pi;
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
    double radius_grense, Radius_Reskalering;
    int B;
    mat Por_index;
    bool who_moves;
    double Porosity;
    bool half_density;
    double density, volume;
    vec3 Reset_Force1;
    bool who_moves_not;
    bool PORER_INNI_TITSLOOPEN;
    int number_of_pores;
    bool senter_pore;

    //Termostat
    string thermostat;
    double relaxation_time;
    double andersen_std;

    //Pressure cellene
    int nLpc;
    double Lpc;
    int nPressCells;
    double volumePressCells;

    //Flow profile
    FlowProfileCells* flowCells;
    int nFlowCells;
    int nflow;
    double Lfc;
};

#endif // GENERATEQUANTITIES_H
