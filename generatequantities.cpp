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
    //b = 5.260; //Prosjekt 1
    b = 5.72; //Prosjekt 2

    //Input half density eller ikke
    half_density = true;

    if (half_density == true){
        b = b * pow(2.0, 1.0/3.0);}
    sigma = 3.405;
    b = b/sigma;
    eV = 1.60217646 * pow(10.0, -19.0);
    epsilon = 1.0318*pow(10.0, -2.0) * eV;
    T0 = 119.74;
    L = b*N;
    volume = L * L * L;
    r1 = zeros(3);
    r2 = zeros(3);
    drvec = zeros(3);
    force = zeros(3);
    radius_grense = 2;
    pi = acos(-1.0);
    Permeability_Counter_pluss = 0;
    Permeability_Counter_minus = 0;
    FlowF << 0 << 0 << 0.1;
    r_cut = 3;
    dt = 0.005;

    //Input

    //Temperatur
    Temp = 1.5*T0;
    T_bath = 1.5*T0;
    time_steps = 2000;
    relaxation_time = 15*dt;
    termostat_limit = 1000; //Må settes til 0 hvis ingen termostat !!!

    //Porer
    who_moves = true; // True = inni porene, False = utenfor porene
    who_moves_not = false; //False = inni porene, True = utenfor porene
    Sylindric_or_circle = true; //True gir sylinder-porer langs z retning, false gir kule-porer
    number_of_pores = 1;
    senter_pore = true;
    PORER_INNI_TITSLOOPEN = false;
    Radius_Reskalering = 2; //Verdi 1 her gir riktig størrelse på porer, større verdi gir mindre porer

    //Flow
    nflow = 35;
    Lfc = L / nflow;
    nFlowCells = nflow * nflow * nflow;
    FlowF *= 50; //Øke kraften litt
    flow_on_off = true; //True gir flow, false gir ikke flow.
    Patrykt_Kraft_Retning = 2; // 0 = x retning, 1 = y retning, 2 = z retning

    //Reskalering;
    density = volume/N_tot;
    Temp /= T0;
    T_bath /= T0;
    standev = sqrt(Temp);
    Lc = (int) (L/r_cut);
    r_cut = L/Lc + 0.0000001;
    nCells = Lc * Lc * Lc;
    idum = -1;
    thermostat = Termostat;
    andersen_std = sqrt(T_bath);
    if (senter_pore == true){
        number_of_pores = 1;}
    if (flow_on_off == false){
        FlowF(0) = 0;
        FlowF(1) = 0;
        FlowF(2) = 0;
    }

    //Pressure cells
    nLpc = 10;//(int) b*N;
    Lpc = L / nLpc;
    nPressCells = nLpc * nLpc * nLpc;
    volumePressCells = Lpc * Lpc * Lpc;
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
        myfile << atoms.size() << endl;
        myfile << "Argon atomer beveger seg" << endl;

        for (i = 0; i < atoms.size(); i++){
            Atom *A = atoms[i];
            if (A->ismoving == false){
                atom_name = "O";
            }
            if (A->ismoving == true){
                atom_name = "Ar";
            }
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

void GenerateQuantities::printPressure(const vec &pressss)
{
    ofstream pressfil ("Pressure.txt", ios_base::app);
    if (pressfil.is_open()){
        for (i=0; i < nPressCells; i++){
            pressfil << pressss(i) << " ";
        }
    pressfil << endl;
    }
    else{
        cout << "open file dude" << endl;
    }
    pressfil.close();
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

    velocitySum /= atoms.size();

    for (i = 0; i< atoms.size(); i++){
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
            atoms[j]->setPressure(dWork);
        }
    }
}

void GenerateQuantities::generateCells()
{
    int ix, iy, iz;
    int nBins = 10;
    double Epot, kinetic, press;

    ofstream Utinfo ("FlowProfile.txt");
    ofstream Utinfo2("FlowProfile2.txt");
    pressCells = new PressureCells[nPressCells];
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
    int test;
    double RADIUS;
    vec Pore_Radiusses(number_of_pores);

    //Pressurecellslister
    vec local_kinetic_energy(nPressCells);
    vec local_pressure(nPressCells);
    vec cellPressure(nPressCells);

    //Andre ting
    mat Index_Matrix;
    Index_Matrix = zeros<mat>(3,27);

    //Flow celler
    int rawNumber = nFlowCells - nflow;
    flowCells = new FlowProfileCells[nFlowCells];
    vec flowVelocityx = zeros(nflow*nflow);
    vec flowVelocityy = zeros(nflow*nflow);
    vec flowVelocityz = zeros(nflow*nflow);
    vec counter = zeros(nflow*nflow);

    //Initialiserer alle atomene i porer

    if (PORER_INNI_TITSLOOPEN == true){
        for (i = 0; i < atoms.size(); i++){
            Atom *atom = atoms[i];
            atom->ismoving = who_moves;
        }
    }


    if (PORER_INNI_TITSLOOPEN == false){
        for (i = 0; i < atoms.size(); i++){
            Atom *atom = atoms[i];
            atom->ismoving = who_moves_not;
        }


    //For hver pore
    Porosity = 0;
    for(j = 0; j < number_of_pores; j++){
        test = 1;
        vec rr2(3);
        while (test != 0){
            test = 0;

            //Lager sentrum til poren

            if (senter_pore == false){
            for (k = 0; k<3; k++){
                rr2(k) = L*ran0(&idum);}
            cout << "Sentrum i poren: " << endl << rr2;
            }

            if (senter_pore == true){
            for (k=0; k<3;k++){
                rr2(k) = L/2;}
            cout << "Sentrum i poren: " << endl << rr2;
            }

            //Lager radius til poren
            RADIUS = ran0(&idum) + radius_grense;
            RADIUS *= 10;
            RADIUS /= sigma;
            RADIUS /= Radius_Reskalering;
            cout << "Radius = " << RADIUS << endl;
            Pore_Radiusses(j) = RADIUS;

            //Sjekker at poren ikke går inn i noen andre porer
            for (i=0; i < atoms.size(); i++){

                //Hvert atom

                //For kule-porer
                if (Sylindric_or_circle == false){
                Atom *atom = atoms[i];
                vec rr1(3);
                vec drr(3);

                rr1 = atom->getPosition();
                drr = rr1 - rr2;

                    for(k = 0; k<3; k++)
                    {
                        if(abs(drr(k)-L) < abs(drr(k)))
                        {
                             drr(k) = drr(k) - L;
                        }
                        if(abs(drr(k)+L) < abs(drr(k)))
                        {
                           drr(k) = drr(k) + L;
                        }
                    }

                    //Sjekker at vi ikke beveger oss
                    if (norm(drr, 2) <= RADIUS){
                        if (atom->ismoving == who_moves){
                            test = 1;
                        }
                    }
                }

                //For sylinderporer
                else if(Sylindric_or_circle == true){
                    Atom *atom = atoms[i];

                    vec rr1(3);
                    vec ddr(3);
                    rr1 = atom->getPosition();
                    ddr = rr1 - rr2;

                    vec rr12d(2);
                    vec ddr2d(2);

                    rr12d << rr1(0) << rr1(0);
                    ddr2d << ddr(0) << ddr(1);

                    for (k = 0; k < 2; k++){
                        if (abs(ddr2d(k) - L) < abs(ddr2d(k))){
                            ddr2d(k) = ddr2d(k) - L;
                        }

                        if (abs(ddr2d(k) + L) < abs(ddr2d(k))){
                            ddr2d(k) = ddr2d(k) + L;
                        }
                    }

                    //Sjekker at vi ikke beveger oss
                    if (norm(ddr2d,2) <= RADIUS){
                        if (atom->ismoving == who_moves){
                            test = 1;
                        }
                    }
                }
            }
        }

        cout << "Pore nummer " << j+1 << " is completed!" << endl;

        if (Sylindric_or_circle == false){
            Porosity += pow(RADIUS, 3.0);
        }

        else if (Sylindric_or_circle == true){
            Porosity += pow(RADIUS, 2.0);
        }

        //Atomene i denne poren beveger seg ikke
        for (i=0; i < atoms.size(); i++){
            Atom *atom = atoms[i];

            if (Sylindric_or_circle == false){
                vec rr1(3);
                rr1 = atom->getPosition();
                vec drr(3);
                drr = rr1 - rr2;

                for(int k = 0; k<3; k++)
                {
                    if(abs(drr(k)-L) < abs(drr(k)))
                    {
                         drr(k) = drr(k) - L;
                    }
                    if(abs(drr(k)+L) < abs(drr(k)))
                    {
                        drr(k) = drr(k) + L;
                    }
                }

                if (norm(drr, 2) <= RADIUS){
                    atom->ismoving = who_moves;
                }
            }

            else if (Sylindric_or_circle == true){
                vec rr1(3);
                rr1 = atom->getPosition();
                vec drr(3);
                drr = rr1-rr2;

                vec rr12d(2);
                vec ddr2d(2);

                rr12d << rr1(0) << rr1(1);
                ddr2d << drr(0) << drr(1);

                for (k = 0; k < 2; k++){
                    if (abs(ddr2d(k) - L) < abs(ddr2d(k))){
                        ddr2d(k) = ddr2d(k) - L;
                    }

                    if (abs(ddr2d(k) + L) < abs(ddr2d(k))){
                        ddr2d(k) = ddr2d(k) + L;
                    }
                }

                if (norm(ddr2d,2) <= RADIUS){
                    atom->ismoving = who_moves;
                }
            }
        }
    }

    if (Sylindric_or_circle == false){
        Porosity = Porosity * 4 / 3 * pi;
    }

    else if (Sylindric_or_circle == true){
        Porosity = Porosity * pi * L;
    }
    Porosity /= volume;

    cout << "Porosity = " << Porosity << endl;
    }

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

    //Generate flowprofilecells
    int flowCellNumber = 0;
            for(int i = 0; i < nflow; i++)
            {
                for(int j = 0; j < nflow; j++)
                {
                    for(int k = 0; k < nflow; k++)
                    {
                        flowCells[flowCellNumber].cell_Index << i << j << k;
                        flowCellNumber++;
                    }
                }
            }

    // Generate the pressure cells
    int pressCellNumber = 0;
    for(int i = 0; i < nLpc; i++)
    {
        for(int j = 0; j < nLpc; j++)
        {
            for(int k = 0; k < nLpc; k++)
            {
                pressCells[pressCellNumber].cell_index << i << j << k;
                pressCellNumber++;
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

    //Setter atomene i pressure cells
    for(i = 0; i < atoms.size(); i++)
    {
        r = atoms[i]->getPosition();
        ix = int (r(0)/Lpc);
        iy = int (r(1)/Lpc);
        iz = int (r(2)/Lpc);
        cellIndices << ix << iy << iz;
        for(j = 0; j < nPressCells; j++)
        {
            if(pressCells[j].cell_index(0)==cellIndices(0) &&
                    pressCells[j].cell_index(1)==cellIndices(1) &&
                    pressCells[j].cell_index(2)==cellIndices(2))
            {
                pressCells[j].addAtoms(atoms[i]);
            }
        }
    }

    //Legger dem til i listene

    generateForce();
    if (PORER_INNI_TITSLOOPEN==false){
        PrintToFile();
    }

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

    temperature(0) = (kinetic_energy(0)) / (3.0 * atoms.size());
    kinetic_energy(0) *= 0.5;
    potential_energy(0) *= 4;
    total_energy(0) = kinetic_energy(0) + potential_energy(0);
    pressure(0) = density * temperature(0) + pressure(0) / (2.0*3.0*volume);


    //Regner trykket internt i pressurecellene
    for(int cellNumber = 0; cellNumber < nPressCells; cellNumber++)
    {
        local_kinetic_energy(cellNumber) = 0.0;
        local_pressure(cellNumber) = 0.0;
        for(int i = 0; i < (pressCells[cellNumber].atoms).size(); i++)
        {
            local_kinetic_energy(cellNumber) += dot(pressCells[cellNumber].atoms[i]->getVelocity(), pressCells[cellNumber].atoms[i]->getVelocity());
            local_pressure(cellNumber) += pressCells[cellNumber].atoms[i]->getPressure();
        }
        cellPressure(cellNumber) = ( local_kinetic_energy(cellNumber) + local_pressure(cellNumber));
        cellPressure(cellNumber) /= (3 * volumePressCells);
    }
    printPressure(cellPressure);


    //***************************************
    //Tidsutvikling
    //***************************************

    for (t = 1; t < time_steps; t++){
        cout << "Timestep number:" << t << endl;

        //Skrur av flowen etter en gitt tid.
        if (t == termostat_limit)
        {
            FlowF(2) = 0;
        }

        //Kalkulerer temp velocity og nye posisjoner for hvert atom
        for (cellNum = 0; cellNum < nCells; cellNum++){
            for (j = 0; j < (cells[cellNum].atoms).size(); j++){
                if (cells[cellNum].atoms[j]->ismoving == who_moves){
                force = cells[cellNum].atoms[j]->getForce();
                tempVelocity = cells[cellNum].atoms[j]->getVelocity() + 0.5 * force * dt;

                oldPosition = cells[cellNum].atoms[j]->getPosition();
                newPosition = oldPosition + tempVelocity * dt;
                d = newPosition - oldPosition;
                    for (k=0; k<3; k++){
                        if (newPosition(k) > L){
                            newPosition(k) -= L;
                            if (k == Patrykt_Kraft_Retning){
                                Permeability_Counter_minus += 1;
                                if (t > termostat_limit){
                                    Utinfo << cells[cellNum].atoms[j]->getPosition()(0) << " " << cells[cellNum].atoms[j]->getPosition()(1) << " " << cells[cellNum].atoms[j]->getPosition()(2) << " " <<  cells[cellNum].atoms[j]->getVelocity()(0) << " "  <<  cells[cellNum].atoms[j]->getVelocity()(1) << " "  <<  cells[cellNum].atoms[j]->getVelocity()(2) << " " << endl;
                                }
                            }
                        }
                        if (newPosition(k) < 0){
                            newPosition(k) += L;
                            if (k == Patrykt_Kraft_Retning){
                                Permeability_Counter_pluss += 1;
                                if (t > termostat_limit){
                                    Utinfo2 << cells[cellNum].atoms[j]->getPosition()(0) << " " << cells[cellNum].atoms[j]->getPosition()(1) << " " << cells[cellNum].atoms[j]->getPosition()(2) << " " <<  cells[cellNum].atoms[j]->getVelocity()(0) << " "  <<  cells[cellNum].atoms[j]->getVelocity()(1) << " "  <<  cells[cellNum].atoms[j]->getVelocity()(2) << " " << endl;
                                }
                            }
                        }
                    }

                    cells[cellNum].atoms[j]->setDisplacement(d);
                    cells[cellNum].atoms[j]->setPosition(newPosition);
                    cells[cellNum].atoms[j]->setVelocity(tempVelocity);
                }
            }
        }

    //Klarerer cellene
    for (cellNum = 0; cellNum < nCells; cellNum++){
        cells[cellNum].atoms.clear();
    }

    //Klarerer flowcellene
    for(int cellNumber = 0; cellNumber < nFlowCells; cellNumber++)
    {
        flowCells[cellNumber].atoms.clear();
    }

    // Set atoms in flow profile cells
        for(int i = 0; i < atoms.size(); i++)
        {
            r = atoms[i]->getPosition();
            ix = int (r(0)/Lfc);
            iy = int (r(1)/Lfc);
            iz = int (r(2)/Lfc);
            cellIndices << ix << iy << iz;
            for(int j = 0; j < nFlowCells; j++)
            {
                if(flowCells[j].cell_Index(0)==cellIndices(0) &&
                        flowCells[j].cell_Index(1)==cellIndices(1) &&
                        flowCells[j].cell_Index(2)==cellIndices(2))
                {
                    flowCells[j].addAtoms(atoms[i]);
                }
            }
        }

    //Klarerer pressurecellene
    for (int cellNumber = 0; cellNumber < nPressCells; cellNumber++){
        pressCells[cellNumber].atoms.clear();
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

    //Resetter kreftene og fysiske størrelser og setter på flow
    kinetic = 0.0;
    for (cellNum = 0; cellNum < nCells; cellNum++){
        cells[cellNum].hasCalculated = false;
        for (j = 0; j < (cells[cellNum].atoms).size(); j++){
            Atom *atom = cells[cellNum].atoms[j];
            if (atom->ismoving == who_moves){
                force = atom->getForce();
                cells[cellNum].atoms[j]->setForce(-force);
                cells[cellNum].atoms[j]->setForce(FlowF);
                Epot = atom->getPotential();
                press = 0.0;
                cells[cellNum].atoms[j]->setPotential(-Epot);
                cells[cellNum].atoms[j]->setPressure(press);
            }
            cells[cellNum].atoms[j]->setKinetic(kinetic);
        }
    }

    //Kalkulerer kreftene mellom atomene i samme celle
    for (cellNum = 0; cellNum < nCells; cellNum++){
        for (i = 0; i < (cells[cellNum].atoms).size(); i++){
            for (j = i+1; j < (cells[cellNum].atoms).size(); j++){
                Atom *atom1 = cells[cellNum].atoms[i];
                Atom *atom2 = cells[cellNum].atoms[j];
                if (atom1->ismoving == who_moves || atom2->ismoving == who_moves){
                    acceleration(atom1, atom2);}
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
                        if (atom1->ismoving == who_moves || atom2->ismoving == who_moves){
                            acceleration(atom1, atom2);
                        }
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

            if (atom->ismoving == who_moves){
                newVelocity = atom->getVelocity() + 0.5*atom->getForce()* dt;
                cells[cellNum].atoms[j]->setVelocity(newVelocity);
                cells[cellNum].atoms[j]->setKinetic(dot(newVelocity, newVelocity));
            }
        }
    }


    // Calculate the flow profile in flow profile cells
        if(t >= termostat_limit)
        {
            int l = 0;
            for(int raw = 0; raw < rawNumber; raw+=nflow)
            {
                for(int j = raw; j < (nflow + raw); j++)
                {
                    for(int k = 0; k < (flowCells[j].atoms).size(); k++)
                    {
                        if(flowCells[j].atoms[k]->ismoving == who_moves)
                        {
                            flowVelocityx[l] += flowCells[j].atoms[k]->getVelocity()(0);
                            flowVelocityy[l] += flowCells[j].atoms[k]->getVelocity()(1);
                            flowVelocityz[l] += flowCells[j].atoms[k]->getVelocity()(2);
                            counter[l] +=1;
                        }
                    }
                }
                l +=1;
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

    //Resetter kreftene og pressure for atomer som IKKE beveger seg
    for (i=0; i < atoms.size(); i++){
        Atom *A = atoms[i];
        if (A->ismoving != who_moves){
            Reset_Force1 = A->getForce();
            A->setPressure(0);
            A->setForce(-Reset_Force1);
        }
    }

    //Finner fysiske størrelser inni tidsloopen
    for(i = 0; i < atoms.size(); i++)
    {
        kinetic_energy(t) += atoms[i]->getKinetic();
        potential_energy(t) += atoms[i]->getPotential();
        pressure(t) += atoms[i]->getPressure();
        displacement(t) += dot(atoms[i]->getDisplacement(), atoms[i]->getDisplacement());
    }

    potential_energy(t) *= 4;
    temperature(t) = kinetic_energy(t)/(3.0 * atoms.size());
    kinetic_energy(t) *= 0.5;
    total_energy(t) = kinetic_energy(t) + potential_energy(t);

    pressure(t) = density * temperature(t) + pressure(t) / (6.0 * volume);
    displacement(t) /= atoms.size();
    DConstant(t) = displacement(t)/(6*t);


    //Setter atomene i pressure cells
    for(int i = 0; i < atoms.size(); i++)
    {
        r = atoms[i]->getPosition();
        ix = int (r(0)/Lpc);
        iy = int (r(1)/Lpc);
        iz = int (r(2)/Lpc);
        cellIndices << ix << iy << iz;
        for(int j = 0; j < nPressCells; j++)
        {
            if(pressCells[j].cell_index(0)==cellIndices(0) &&
                    pressCells[j].cell_index(1)==cellIndices(1) &&
                    pressCells[j].cell_index(2)==cellIndices(2))
            {
                pressCells[j].addAtoms(atoms[i]);
            }
        }
    }

    //Regner trykket internt i pressurecellene
    for(int cellNumber = 0; cellNumber < nPressCells; cellNumber++)
    {
        local_kinetic_energy(cellNumber) = 0.0;
        local_pressure(cellNumber) = 0.0;
        for(int i = 0; i < (pressCells[cellNumber].atoms).size(); i++)
        {
            if (pressCells[cellNumber].atoms[i]->ismoving == who_moves){
                local_kinetic_energy(cellNumber) += pressCells[cellNumber].atoms[i]->getKinetic();
                local_pressure(cellNumber) += pressCells[cellNumber].atoms[i]->getPressure();
            }
        }
        cellPressure(cellNumber) = ( local_kinetic_energy(cellNumber) + local_pressure(cellNumber));
        cellPressure(cellNumber) /= (3 * volumePressCells);
    }

    printPressure(cellPressure);
    /////////////////////////////////////////////////////////

    //Hvis vi setter porene inni tidsloopen
    if (PORER_INNI_TITSLOOPEN == true){
    if (t == termostat_limit){

        for (i = 0; i < atoms.size(); i++){
            Atom *atom = atoms[i];
            atom->ismoving = who_moves_not;
        }


    //For hver pore
    Porosity = 0;
    for(j = 0; j < number_of_pores; j++){
        test = 1;
        vec rr2(3);
        while (test != 0){
            test = 0;

            //Lager sentrum til poren
            for (k = 0; k<3; k++){
                rr2(k) = L*ran0(&idum);}

            //Lager radius til poren
            RADIUS = ran0(&idum) + radius_grense;
            RADIUS *= 10;
            RADIUS /= sigma;
            RADIUS /= Radius_Reskalering;
            Pore_Radiusses(j) = RADIUS;

            //Sjekker at poren ikke går inn i noen andre porer
            for (i=0; i < atoms.size(); i++){

                //Hvert atom
                //For kule-porer
                if (Sylindric_or_circle == false){
                Atom *atom = atoms[i];
                vec rr1(3);
                vec drr(3);

                rr1 = atom->getPosition();
                drr = rr1 - rr2;

                    for(k = 0; k<3; k++)
                    {
                        if(abs(drr(k)-L) < abs(drr(k)))
                        {
                             drr(k) = drr(k) - L;
                        }
                        if(abs(drr(k)+L) < abs(drr(k)))
                        {
                           drr(k) = drr(k) + L;
                        }
                    }

                    //Sjekker om vi beveger oss
                    if (norm(drr, 2) <= RADIUS){
                        if (atom->ismoving == who_moves){
                            test = 1;
                        }
                    }
                }

                //For sylinderporer
                else if(Sylindric_or_circle == true){
                    Atom *atom = atoms[i];

                    vec rr1(3);
                    vec ddr(3);
                    rr1 = atom->getPosition();
                    ddr = rr1 - rr2;

                    vec rr12d(2);
                    vec ddr2d(2);

                    rr12d << rr1(0) << rr1(0);
                    ddr2d << ddr(0) << ddr(1);

                    for (k = 0; k < 2; k++){
                        if (abs(ddr2d(k) - L) < abs(ddr2d(k))){
                            ddr2d(k) = ddr2d(k) - L;
                        }

                        if (abs(ddr2d(k) + L) < abs(ddr2d(k))){
                            ddr2d(k) = ddr2d(k) + L;
                        }
                    }

                    //Sjekker om vi beveger oss
                    if (norm(ddr2d,2) <= RADIUS){
                        if (atom->ismoving == who_moves){
                            test = 1;
                        }
                    }
                }
            }
        }

        cout << "Pore nummer " << j+1 << " is completed!" << endl;

        if (Sylindric_or_circle == false){
            Porosity += pow(RADIUS, 3.0);
        }

        else if (Sylindric_or_circle == true){
            Porosity += pow(RADIUS, 2.0);
        }

        //Atomene i denne poren beveger seg ikke
        for (i=0; i < atoms.size(); i++){
            Atom *atom = atoms[i];

            if (Sylindric_or_circle == false){
                vec rr1(3);
                rr1 = atom->getPosition();
                vec drr(3);
                drr = rr1 - rr2;

                for(int k = 0; k<3; k++)
                {
                    if(abs(drr(k)-L) < abs(drr(k)))
                    {
                         drr(k) = drr(k) - L;
                    }
                    if(abs(drr(k)+L) < abs(drr(k)))
                    {
                        drr(k) = drr(k) + L;
                    }
                }

                if (norm(drr, 2) <= RADIUS){
                    atom->ismoving = who_moves;
                }
            }

            else if (Sylindric_or_circle == true){
                vec rr1(3);
                rr1 = atom->getPosition();
                vec drr(3);
                drr = rr1-rr2;

                vec rr12d(2);
                vec ddr2d(2);

                rr12d << rr1(0) << rr1(1);
                ddr2d << drr(0) << drr(1);

                for (k = 0; k < 2; k++){
                    if (abs(ddr2d(k) - L) < abs(ddr2d(k))){
                        ddr2d(k) = ddr2d(k) - L;
                    }

                    if (abs(ddr2d(k) + L) < abs(ddr2d(k))){
                        ddr2d(k) = ddr2d(k) + L;
                    }
                }

                if (norm(ddr2d,2) <= RADIUS){
                    atom->ismoving = who_moves;
                }
            }
        }
    }

    if (Sylindric_or_circle == false){
        Porosity = Porosity * 4 / 3 * pi;
    }

    else if (Sylindric_or_circle == true){
        Porosity = Porosity * pi * L;
    }
    Porosity /= volume;

    cout << "Porosity = " << Porosity << endl;
    }}
    //////////////////////////

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

        if (PORER_INNI_TITSLOOPEN == true){
            if (t >= termostat_limit){
                PrintToFile();
            }
        }


        else{
            PrintToFile();
        }
    }


    //****************************************
    //Slutt på tidsloopen
    //****************************************

    Utinfo.close();
    Utinfo2.close();

    double Permeability;
    double r_sum_squared=0;
    for (j = 0; j < number_of_pores; j++){
        r_sum_squared += Pore_Radiusses(j) * Pore_Radiusses(j);
    }
    Permeability = Porosity * r_sum_squared / 8;
    cout << "Permeability = " << Permeability*sigma*sigma << endl;


    //g(r) normalisering
    for (i = 0; i < nBins; i++){
        bins(i) /= (time_steps*atoms.size()*atoms.size());
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
    cout << "Permeability counter pluss kom til: " << Permeability_Counter_pluss << endl;
    cout << "Permeability counter minus kom til: " << Permeability_Counter_minus << endl;
    cout << "Porøsitet var: " << Porosity << endl;

    //Printe flow cellene
    ofstream ofil ("FlowProfile.txt");
    if (ofil.is_open())
    {
        for(int l = 0; l < nflow*nflow; l++)
        {
            if(counter[l] == 0)
            {
                ofil << flowVelocityx[l] <<" " << flowVelocityy[l] <<" "<< flowVelocityz[l]<<endl;
            }else{
                ofil << flowVelocityx[l]/counter[l] <<" " << flowVelocityy[l]/counter[l] <<" "<< flowVelocityz[l]/counter[l] <<endl;
            }
        }
        ofil.close();
    }
    else cout << "Unable to open file";

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

    double normtemp = norm(drvec,2.0);

    normtemp = 1/(normtemp*normtemp*normtemp);

    Ep = (normtemp - 1)*normtemp;

    force = potential->Lennard_Jones_potential(drvec);

    p = dot(force, drvec);

    atom1->setForce(force);
    atom2->setForce(-force);

    atom1->setPotential(Ep);
    atom2->setPotential(Ep);

    atom1->setPressure(p);
    atom2->setPressure(p);
}

void GenerateQuantities::ClearCells()
{
    for (int cellNum = 0; cellNum < nCells; cellNum++){
        cells[cellNum].atoms.clear();
    }
}

double GenerateQuantities::radialDF(Atom *atom1, Atom *atom2)
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
