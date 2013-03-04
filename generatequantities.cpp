#include "generatequantities.h"

GenerateQuantities::GenerateQuantities(string Termostat)
{
    potential = new Potentials();

    //Deklarasjoner
    N = 8;
    N_tot = N*N*N*4;
    atom_name = "Ar";
    m = 39.948; //amu
    k_b = 1.3806503 * pow(10.0, -23.0);
    //b = 5.260;
    b = 5.72;
    sigma = 3.405;
    b = b/sigma;
    eV = 1.60217646 * pow(10.0, -19.0);
    epsilon = 1.0318*pow(10.0, -2.0) * eV;
    T0 = 119.74;
    L = b*N;
    volume = L * L * L;
    density = volume/N_tot;
    r1 = zeros(3);
    r2 = zeros(3);
    drvec = zeros(3);
    force = zeros(3);

    //Input
    Temp = 100.0;
    T_bath = 0.851*T0;
    r_cut = 3;
    time_steps = 1000;
    dt = 0.005;
    relaxation_time = 15*dt;
    termostat_limit = 500;

    //Reskalering;
    Temp /= T0;
    T_bath /= T0;
    standev = sqrt(Temp);
    Lc = (int) (L/r_cut);
    r_cut = L/Lc + 0.0000001;
    nCells = Lc * Lc * Lc;
    idum = -1;
    thermostat = Termostat;
    andersen_std = sqrt(T_bath);
}

double GenerateQuantities::NormalDist(long *idum)
{
    double u = ran0(idum) * 2 - 1;
    double v = ran0(idum) * 2 - 1;
    double r = u*u + v*v;
    if (r==0 || r > 1){
        return NormalDist(idum);}
    else{
        double c = sqrt(-2 * log(r)/r);
        return u*c;
    }
}

void GenerateQuantities::PrintToFile()
{
    ofstream myfile ("Argon6.xyz", ios_base::app);
    if (myfile.is_open()){
        myfile << N_tot << endl;
        myfile << "Argon atomer beveger seg" << endl;

        for (i = 0; i < atoms.size(); i++){
            Atom *A = atoms[i];
            myfile << atom_name << " " << A->getPosition()(0) << " " <<  A->getPosition()(1) << " " << A->getPosition()(2) << endl;
        }
        myfile.close();
    }
    else{
        cout << "open file dude" << endl;
    }
}

void GenerateQuantities::printVelocity()
{
    ofstream minfil;
    minfil.open("velocity.txt");
     if (minfil.is_open()){
         for (i=0; i<atoms.size(); i++){
             Atom *A = atoms[i];
             minfil << A->getVelocity()(0) << " " << A->getVelocity()(1) << " " << A->getVelocity()(2) << endl;
         }
         minfil.close();
     }
     else{
         cout << "open file dude" << endl;
     }
 }


void GenerateQuantities::generatePosition()
{
    Atom *A;
    vec R(3);
    vec r(3);

    //Lager krystallstrukturen
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            for (k=0;k<N;k++){
                //Lager de 4 atomene i hver av krystallene våre

                R(0) = k*b;
                R(1) = j*b;
                R(2) = i*b;

                A = new Atom(R);
                atoms.push_back(A);

                r(0) = R(0) + 0.5*b;
                r(1) = R(1) + 0.5*b;
                r(2) = R(2);

                A = new Atom(r);
                atoms.push_back(A);

                r(0) = R(0);
                r(1) = R(1) + 0.5*b;
                r(2) = R(2) + 0.5*b;

                A = new Atom(r);
                atoms.push_back(A);

                r(0) = R(0) + 0.5*b;
                r(1) = R(1);
                r(2) = R(2) + 0.5*b;
                A = new Atom(r);
                atoms.push_back(A);
            }
        }
    }
}

void GenerateQuantities::generateVelocity()
{
    vec velocity(3);
    vec3 velocitySum = zeros(3);

    //Gir alle atomene en normalfordelt hastighet std = 1;
    for (i = 0; i < atoms.size(); i++){

        //Normalfordeling
        velocity(0) = NormalDist(&idum);
        velocity(1) = NormalDist(&idum);
        velocity(2) = NormalDist(&idum);

        //Uniform fordeling
        //velocity(0) = randu();
        //velocity(1) = randu();
        //velocity(2) = randu();
        //velocity *= sqrt(12.0);
        //////////////////////////////////

        velocity *= standev;

        Atom* A = atoms.at(i);
        A->setVelocity(velocity);
        velocitySum += A->getVelocity();
    }

    velocitySum /= N_tot;

    for (i = 0; i<atoms.size(); i++){
        Atom* A = atoms.at(i);
        velocity = A->getVelocity() - velocitySum;
        A->setVelocity(velocity);
    }
}

void GenerateQuantities::generateForce()
{
    double norm;
    double Epot;
    double dWork;
    vec force(3);
    vec DeltaR_ij(3);

    for (i=0;i<atoms.size();i++){
        for (j=i+1;j<atoms.size();j++){

            //Finner minste avstand mellom atomene
            DeltaR_ij = atoms[i]->getPosition() - atoms[j]->getPosition();
            for (k=0; k<3; k++){
                temporary = DeltaR_ij(k) - L;
                if (abs(temporary) < abs(DeltaR_ij(k))){
                    DeltaR_ij(k) = temporary;
                }

                temporary = ( DeltaR_ij(k) + L);
                if (abs(temporary) < abs(DeltaR_ij(k))){
                    DeltaR_ij(k) = temporary;
                }
            }

            //Beregner kraften
            norm = dot(DeltaR_ij, DeltaR_ij);
            force = (2.0 * pow(norm, -3.0) - 1) * DeltaR_ij * pow(norm, -4.0) * 24;
            Epot = (pow( norm, -6.0) - pow(norm, -3.0));
            dWork = dot(force, DeltaR_ij);

            //Bruker newtons tredje lov
            atoms[i]->setForce(force);
            atoms[j]->setForce(-force);

            atoms[i]->setPotential(Epot);
            atoms[j]->setPotential(Epot);

            atoms[i]->setPressure(dWork);
        }
    }
}

void GenerateQuantities::generateCells()
{
    int ix, iy, iz;
    int nBins = 10;
    double Epot, kinetic, press;

    cells = new Cell[nCells];

    vec r(3);
    vec d(3);
    vec cellIndices(3);
    vec indexes(3);
    vec3 newVelocity;
    vec3 oldPosition;
    vec3 newPosition;
    vec3 tempVelocity;
    vec3 force;
    vec kinetic_energy(time_steps);
    vec potential_energy(time_steps);
    vec total_energy(time_steps);
    vec temperature(time_steps);
    vec pressure(time_steps);
    vec displacement(time_steps);
    vec bins(nBins);
    vec DConstant(time_steps);

    mat Index_Matrix;
    Index_Matrix = zeros<mat>(3,27);

    //Lager cellene. Hver celle får et nummer og har 3 indekser som gir posisjonen til cella.
    int cellNum = 0;
    for (i = 0; i < Lc; i++){
        for (j = 0; j < Lc; j++){
            for (k = 0; k < Lc; k++){
                cells[cellNum].cellIndices << i << j << k;
                cellNum++;
            }
        }
    }

    //Konstruerer en index matrise for å finne naboene senere
    int l=0;
    for (i = -1; i < 2; i++){
        for (j = -1; j < 2; j++){
            for (k = -1; k < 2; k++){
                Index_Matrix(0,l) = i;
                Index_Matrix(1,l) = j;
                Index_Matrix(2,l) = k;
                l += 1;
            }
        }
    }

    //Finner nabocellene og legger de til i matrisa
    for (i = 0; i < nCells; i++){
        for (j = 0; j<Index_Matrix.n_cols; j++){
            indexes = Index_Matrix.col(j) + cells[i].cellIndices;
            for (k = 0; k < 3; k++){
                if (indexes(k) > Lc-1){
                    indexes(k) -= Lc;
                }
                else if (indexes(k) < 0){
                    indexes(k) += Lc;
                }
            }

            for (k = 0; k < nCells; k++){
                if (cells[k].cellIndices(0) == indexes(0) &&
                        cells[k].cellIndices(1) == indexes(1) &&
                        cells[k].cellIndices(2) == indexes(2) &&
                        i != k){
                    cells[i].addNeighbor(&cells[k]);
                }
            }
        }
    }

    //Sett atomer i cellene
    for (i = 0; i < atoms.size(); i++){
        r = atoms[i]->getPosition();

        ix = (int) (r(0)/r_cut);
        iy = (int) (r(1)/r_cut);
        iz = (int) (r(2)/r_cut);

        cellIndices << ix << iy << iz;

        for (j = 0; j < nCells; j++){
            if (cells[j].cellIndices(0) == cellIndices(0) &&
                    cells[j].cellIndices(1) == cellIndices(1) &&
                    cells[j].cellIndices(2) == cellIndices(2)){
                cells[j].addAtoms(atoms[i]);
            }
        }
    }

    //Legger dem til i listene

    generateForce();
    PrintToFile();

    //Initialiserer g(r)
    for (i = 0; i<nBins; i++){
        bins(i) = 0;
    }

    //Finner fysiske størrelser på t=0
    kinetic_energy(0) = 0;
    potential_energy(0) = 0;
    displacement(0) = 0;
    pressure(0) = 0;
    DConstant(0) = 0;

    for (i = 0; i < atoms.size(); i++){
        Atom *atom = atoms[i];
        kinetic_energy(0) += dot(atom->getVelocity(), atom->getVelocity());
        potential_energy(0) += atom->getPotential();
        pressure(0) += atom->getPressure();
    }

    temperature(0) = (kinetic_energy(0)) / (3.0 * N_tot);
    kinetic_energy(0) *= 0.5;
    potential_energy(0) *= 4;
    total_energy(0) = kinetic_energy(0) + potential_energy(0);
    pressure(0) = density * temperature(0) + pressure(0) / (3.0*volume);

    //***************************************
    //Tidsutvikling
    //***************************************

    for (t = 1; t < time_steps; t++){
        cout << "Timestep number:" << t << endl;

        //Kalkulerer temp velocity og nye posisjoner for hvert atom
        for (cellNum = 0; cellNum < nCells; cellNum++){
            for (j = 0; j < (cells[cellNum].atoms).size(); j++){
                force = cells[cellNum].atoms[j]->getForce();
                tempVelocity = cells[cellNum].atoms[j]->getVelocity() + 0.5 * force * dt;
                oldPosition = cells[cellNum].atoms[j]->getPosition();
                newPosition = oldPosition + tempVelocity * dt;
                d = newPosition - oldPosition;

                for (k=0; k<3; k++){
                    if (newPosition(k) > L){
                        newPosition(k) -= L;
                    }
                    if (newPosition(k) < 0){
                        newPosition(k) += L;
                    }
                }

                cells[cellNum].atoms[j]->setPosition(newPosition);
                cells[cellNum].atoms[j]->setVelocity(tempVelocity);

                cells[cellNum].atoms[j]->setDisplacement(d);
            }
        }

    //Klarerer cellene
    for (cellNum = 0; cellNum < nCells; cellNum++){
        cells[cellNum].atoms.clear();
    }

    //Setter atomer inn i cellene
    for (i=0; i<atoms.size(); i++){
        r = atoms[i]->getPosition();

        ix = (int) (r(0)/r_cut);
        iy = (int) (r(1)/r_cut);
        iz = (int) (r(2)/r_cut);

        cellIndices << ix << iy << iz;

        for (j=0; j<nCells; j++){
            if (    cells[j].cellIndices(0) == cellIndices(0) &&
                    cells[j].cellIndices(1) == cellIndices(1) &&
                    cells[j].cellIndices(2) == cellIndices(2)){
                cells[j].addAtoms(atoms[i]);
            }
        }
    }

    //Resetter kreftene og fysiske størrelser
    for (cellNum = 0; cellNum < nCells; cellNum++){
        cells[cellNum].hasCalculated = false;
        for (j = 0; j < (cells[cellNum].atoms).size(); j++){
            Atom *atom = cells[cellNum].atoms[j];
            force = atom->getForce();
            cells[cellNum].atoms[j]->setForce(-force);

            Epot = atom->getPotential();
            kinetic = 0.0;
            press = 0.0;
            cells[cellNum].atoms[j]->setPotential(-Epot);
            cells[cellNum].atoms[j]->setKinetic(kinetic);
            cells[cellNum].atoms[j]->setPressure(press);
        }
    }

    //Kalkulerer kreftene mellom atomene i samme celle
    for (cellNum = 0; cellNum < nCells; cellNum++){
        for (i = 0; i < (cells[cellNum].atoms).size(); i++){
            for (j = i+1; j < (cells[cellNum].atoms).size(); j++){
                Atom *atom1 = cells[cellNum].atoms[i];
                Atom *atom2 = cells[cellNum].atoms[j];
                acceleration(atom1, atom2);
            }
        }

        //Kalkulerer kreftene mellom atomene i nabocellene også

        for (p = 0; p < (cells[cellNum].neighbor).size(); p++){
            Cell* neighbor = cells[cellNum].neighbor[p];
            if (!neighbor->hasCalculated == false){
                for (i = 0; i < (cells[cellNum].atoms).size();i++){
                    for (j = 0; j < (neighbor->atoms).size(); j++){
                        Atom *atom1 = cells[cellNum].atoms[i];
                        Atom *atom2 = neighbor->atoms[j];
                        acceleration(atom1, atom2);
                    }
                }
            }
        }
        cells[cellNum].hasCalculated = true;
    }

    //Kalkulerer nye velocities
    for (cellNum = 0; cellNum < nCells; cellNum++){
        for (j = 0; j < (cells[cellNum].atoms).size(); j++){
            Atom *atom = cells[cellNum].atoms[j];

            newVelocity = atom->getVelocity() + 0.5*atom->getForce()* dt;
            cells[cellNum].atoms[j]->setVelocity(newVelocity);

            cells[cellNum].atoms[j]->setKinetic(dot(newVelocity, newVelocity));
        }
    }

    //Finner g(r)
    for (i= 0; i< atoms.size(); i++){
        for (j = 0; j < atoms.size(); j++){
            Atom *atom1 = atoms[i];
            Atom *atom2 = atoms[j];
            dr = abs(radialDF(atom1, atom2));
            BinIndex = nBins * (int) 2*dr/L;
            if (BinIndex >= nBins){}
            else {
                bins(BinIndex) += 1;
            }
        }
    }

    //Kalkulerer fysiske størrelser
    kinetic_energy(t) = 0.0;
    potential_energy(t) = 0.0;
    pressure(t)         = 0.0;
    displacement(t)     = 0.0;

    for(int i = 0; i < atoms.size(); i++)
    {
        kinetic_energy(t) += atoms[i]->getKinetic();
        potential_energy(t) += atoms[i]->getPotential();
        pressure(t) += atoms[i]->getPressure();
        displacement(t) += dot(atoms[i]->getDisplacement(), atoms[i]->getDisplacement());
    }

    potential_energy(t) *= 4;
    temperature(t) = kinetic_energy(t)/(3.0 * N_tot);
    kinetic_energy(t) *= 0.5;
    total_energy(t) = kinetic_energy(t) + potential_energy(t);
    pressure(t) = density * temperature(t) + pressure(t) / (3.0 * volume);
    displacement(t) /= N_tot;
    DConstant(t) = displacement(t)/(6*t);

    //Bruker termostaten
    if (thermostat == "noThermostat"){}

    else if (thermostat == "Berendsen"){
        if (t == termostat_limit){
            thermostat = "noThermostat";
        }
        for (cellNum = 0; cellNum < nCells; cellNum++)
                {
                    for(j = 0; j < (cells[cellNum].atoms).size(); j++)
                    {
                        Atom *atom = cells[cellNum].atoms[j];
                        newVelocity = atom->getVelocity();
                        newVelocity *= BerendsenThermostat(temperature(t));
                        cells[cellNum].atoms[j]->setVelocity(newVelocity);
                    }
        }
    }

    else if(thermostat == "Andersen"){
        if (t == termostat_limit){
            thermostat = "noThermostat";
        }
        for (cellNum = 0; cellNum < nCells; cellNum++){
            for (j=0; j < (cells[cellNum].atoms).size(); j++){
                Atom *atom = cells[cellNum].atoms[j];
                AndersenThermostat(atom);
            }
        }
    }


    PrintToFile();
    }

    //****************************************
    //Slutt på tidsloopen
    //****************************************


    //g(r) normalisering
    for (i = 0; i < nBins; i++){
        bins(i) /= (time_steps*N_tot*N_tot);
    }

    //Skriver ut fysiske egenskaper


    // Write dynamical properties to file
    ofstream myfilz ("Energy.txt");
    if (myfilz.is_open())
    {
        for(int t = 0; t < time_steps; t++)
        {
            myfilz << t << " " << total_energy(t) << " " << temperature(t) << " "<< kinetic_energy(t)
                   << " " << potential_energy(t) << " " << pressure(t)<< " " << displacement(t)*sigma*sigma << " " << DConstant(t) << endl;
        }
        myfilz.close();
    }
    else cout << "Unable to open file";

    ofstream myfil ("radial.txt");
    if (myfil.is_open())
    {
        for(int t = 0; t < nBins; t++)
        {
            myfil << bins(t) << endl;
        }
        myfil.close();
    }
    else cout << "Unable to open file";

    printVelocity();
}

void GenerateQuantities::acceleration(Atom *atom1, Atom *atom2)
{
    double Ep, p;

    r1 = atom1->getPosition();
    r2 = atom2->getPosition();

    drvec = r1 - r2;

    //Finner minste avstand mellom atomene
    for (k=0; k<3; k++){
        if( abs(drvec(k)-L) < abs(drvec(k))) {
            drvec(k) = drvec(k) - L;
        }

        if( abs(drvec(k)+L) < abs(drvec(k))) {
            drvec(k) = drvec(k) + L;
        }
    }

    Ep = (pow(norm(drvec,2), -3.0) - 1)*pow(norm(drvec,2), -3.0);
    p = dot(force, drvec);
    force = potential->Lennard_Jones_potential(drvec);

    atom1->setForce(force);
    atom2->setForce(-force);

    atom1->setPotential(Ep);
    atom2->setPotential(Ep);

    atom1->setPressure(p);
    //atom2->setPressure(p);

}

void GenerateQuantities::ClearCells()
{
    for (int cellNum = 0; cellNum < nCells; cellNum++){
        cells[cellNum].atoms.clear();
    }
}

double GenerateQuantities::                                                                                                                                                                                                                                                         F(Atom *atom1, Atom *atom2)
{
    vec r1(3);
    vec r2(3);
    vec dr(3);
    r1 = atom1->getPosition();
    r2 = atom2->getPosition();

    dr = r1 - r2;
    for(int k = 0; k<3; k++)
    {
        if(abs(dr(k)-L) < abs(dr(k)))
        {
             dr(k) = dr(k) - L;
        }
        if(abs(dr(k)+L) < abs(dr(k)))
        {
            dr(k) = dr(k) + L;
        }
    }
    return norm(dr,2);
}

double GenerateQuantities::BerendsenThermostat(double Temperature)
{
    double gamma;
    double temporary_variable = (T_bath / Temperature) - 1.0;
    temporary_variable *= (dt/relaxation_time);
    gamma = sqrt(1.0 + temporary_variable);
    return gamma;
}

void GenerateQuantities::AndersenThermostat(Atom *atom)
{
    double random_nymber = (double) rand()/RAND_MAX;
    vec velocity(3);
    if (random_nymber < dt/relaxation_time){
        velocity(0) = NormalDist(&idum);
        velocity(1) = NormalDist(&idum);
        velocity(2) = NormalDist(&idum);
        velocity *= andersen_std;
        atom->setVelocity(velocity);
    }
}
