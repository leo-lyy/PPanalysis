# include "../include/rg_calculation.h"

using namespace std;

void dumpIO_Rg(System& system, ifstream& dumpFilein, int tarchain, int endchain)
{
    // molMass [molID], define a vector to store the mass of each chain
    vector<double> molMass(system.molnumMAX+1,0);
    // msrg [chain number]
    vector<double> msrg(system.molnumMAX+1,0);
    // define a vector to store the center of mass in x, y, z direction
    vector<vector<double>> com(system.molnumMAX+1, vector<double>(4,0));
    string line;
    istringstream ss(line);
    ofstream msrgout("msrg.txt");
    ofstream comXout("comX.txt");
    ofstream comYout("comY.txt");
    ofstream comZout("comZ.txt");
    // for better readability, we can add headers to the output files, ignore for better performance in MATLAB
    // msrgout << "Step        ";
    // for (int ii = tarchain; ii <= endchain; ii++) msrgout << "Mol_ID" << ii << "      ";
    // msrgout << endl;
    // comXout << "Step        ";
    // for (int ii = tarchain; ii <= endchain; ii++) comXout << "Mol_ID" << ii << "      ";
    // comXout << endl;
    // comYout << "Step        ";
    // for (int ii = tarchain; ii <= endchain; ii++) comYout << "Mol_ID" << ii << "      ";
    // comYout << endl;
    // comZout << "Step        ";
    // for (int ii = tarchain; ii <= endchain; ii++) comZout << "Mol_ID" << ii << "      ";
    // comZout << endl;


    long int f = 0;
    while (f < system.frames)
    {
        // long int iderr = 0;
        // auto fstart = std::chrono::high_resolution_clock::now();
        f++;
        long int id, flable, ix, iy, iz;
        int mol, type;
        double x, y, z, vx, vy, vz;
        for (long int k = 0; k < 9; k++)
        {
            getline(dumpFilein, line);

            if (line == "ITEM: TIMESTEP")
            {
                getline(dumpFilein, line);
                ss >> flable;
                k++;
            }
            else if (line == "ITEM: BOX BOUNDS pp pp pp")
            {
                system.box_type = 0; // orthogonal box
                getline(dumpFilein, line);
                ss >> system.xlo >> system.xhi;
                getline(dumpFilein, line);
                ss >> system.ylo >> system.yhi;
                getline(dumpFilein, line);
                ss >> system.zlo >> system.zhi;
                system.bx = system.xhi - system.xlo;
                system.by = system.yhi - system.ylo;
                system.bz = system.zhi - system.zlo;
                k += 3;
            }
            else if (line == "ITEM: BOX BOUNDS xy xz yz pp pp pp")
            {
                system.box_type = 1; // triclinic box
                getline(dumpFilein, line);
                ss >> system.xlo >> system.xhi >> system.xtilt;
                getline(dumpFilein, line);
                ss >> system.ylo >> system.yhi >> system.ytilt;
                getline(dumpFilein, line);
                ss >> system.zlo >> system.zhi >> system.ztilt;
                system.bx = system.xhi - system.xlo;
                system.by = system.yhi - system.ylo;
                system.bz = system.zhi - system.zlo;
                k += 3;
            }
            else if (line == "ITEM: ATOMS id mol type x y z ix iy iz")
            {
                for (long int i = 0; i < system.num_atoms; i++)
                {
                    dumpFilein >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz;
                    id--;  // the atom id in dump file starts from 1
                    system.atoms[i].mol = mol;
                    system.atoms[i].type = type;
                    system.atoms[i].x = x;
                    system.atoms[i].y = y;
                    system.atoms[i].z = z;
                    system.atoms[i].ix = ix;
                    system.atoms[i].iy = iy;
                    system.atoms[i].iz = iz;
                    if (f == 1) molMass[mol] += system.mass[type];
                }
                // if (iderr != 0) cout << "ID error: " << iderr <<", frame:"<<f<< endl;
            }
            else if (line == "ITEM: ATOMS id mol type x y z ix iy iz vx vy vz")
            {
                for (long int i = 0; i < system.num_atoms; i++)
                {
                    dumpFilein >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz >> vx >> vy >> vz;
                    id--;  // the atom id in dump file starts from 1
                    system.atoms[i].mol = mol;
                    system.atoms[i].type = type;
                    system.atoms[i].x = x;
                    system.atoms[i].y = y;
                    system.atoms[i].z = z;
                    system.atoms[i].ix = ix;
                    system.atoms[i].iy = iy;
                    system.atoms[i].iz = iz;
                    if (f == 1) molMass[mol] += system.mass[type];
                }
            }

        }
        unwrap(system);
        // cout << "Reading frame " << f << endl;

        // Calculate the radius of gyration
        // initialize the com and msrg
        for (int i = tarchain; i <= endchain; i++)
        {
            com[i][1] = 0;
            com[i][2] = 0;
            com[i][3] = 0;
            msrg[i] = 0;
        }
        for (long int i = 0; i < system.num_atoms; i++)
        {
            mol = system.atoms[i].mol;
            if (mol >= tarchain && mol <= endchain)
            {
                com[mol][1] += system.atoms[i].x * system.mass[system.atoms[i].type];
                com[mol][2] += system.atoms[i].y * system.mass[system.atoms[i].type];
                com[mol][3] += system.atoms[i].z * system.mass[system.atoms[i].type];
            }
        }
        for (int i = tarchain; i <= endchain; i++)
        {
            com[i][1] /= molMass[i];
            com[i][2] /= molMass[i];
            com[i][3] /= molMass[i];
        }
        for (long int j = 0; j < system.num_atoms; j++)
        {
            int k = system.atoms[j].mol;
            if (k >= tarchain && k <= endchain)
            {
                msrg[k] += system.mass[system.atoms[j].type] * ((system.atoms[j].x - com[k][1]) * (system.atoms[j].x - com[k][1]) + 
                                                               (system.atoms[j].y - com[k][2]) * (system.atoms[j].y - com[k][2]) +
                                                               (system.atoms[j].z - com[k][3]) * (system.atoms[j].z - com[k][3]));
            }
        }
        // wrap the com into the box make sure the com fall in the box lo and hi
        for (int i = tarchain; i <= endchain; i++)
        {
            while (com[i][1] < system.xlo) com[i][1] += system.bx;
            while (com[i][1] >= system.xhi) com[i][1] -= system.bx;
            while (com[i][2] < system.ylo) com[i][2] += system.by;
            while (com[i][2] >= system.yhi) com[i][2] -= system.by;
            while (com[i][3] < system.zlo) com[i][3] += system.bz;
            while (com[i][3] >= system.zhi) com[i][3] -= system.bz;
        }
        for (int i = tarchain; i <= endchain; i++) msrg[i] = sqrt(msrg[i] / molMass[i]);
        msrgout << fixed << setprecision(0) << f << "        ";
        comXout << fixed << setprecision(0) << f << "        ";
        comYout << fixed << setprecision(0) << f << "        "; 
        comZout << fixed << setprecision(0) << f << "        ";
        
        msrgout << fixed << setprecision(5);
        comXout << fixed << setprecision(5);
        comYout << fixed << setprecision(5);
        comZout << fixed << setprecision(5);
        
        for (int i = tarchain; i <= endchain; i++)
        {
            msrgout << msrg[i] << "      ";
            comXout << com[i][1] << "      ";
            comYout << com[i][2] << "      ";
            comZout << com[i][3] << "      ";
        }
        msrgout << endl;
        comXout << endl;
        comYout << endl;
        comZout << endl;
        printProgressBar(f, system.frames);
    }
    msrgout.close();
    comXout.close();
    comYout.close();
    comZout.close();
}
