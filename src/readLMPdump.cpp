
#include"../include/readLMPdump.h"

using namespace std;


long int countDumpFrame(const string& dumpFileName)
{
    ifstream dumpFile(dumpFileName);
    string line;
    long int frameCount = 0;
    while (getline(dumpFile, line))
    {
        if (line == "ITEM: TIMESTEP")
        {
            frameCount++;
        }
    }
    dumpFile.close();
    return frameCount;
}

void dumpIO_p2(System& system, ifstream& dumpFilein)
{
    string line;
    istringstream ss(line);
    ofstream p2out("P2.txt");
    long int f = 0;
    while (f < system.frames)
    {
        // long int iderr = 0;
        auto fstart = std::chrono::high_resolution_clock::now();
        f++;
        cout << endl << "Calculating frame " << f << "/" << system.frames << " ..." << endl;
        long int id, flable;
        int mol, type;
        double x, y, z, vx, vy, vz;
        long int ix, iy, iz;
        // update the atom position and velocity
        while (line != "ITEM: ATOMS id mol type x y z ix iy iz vx vy vz")
        {
            getline(dumpFilein, line);
            // cout << line << endl;
        }
        for (long int i = 0; i < system.num_atoms; i++)
        {
            getline(dumpFilein, line);
            istringstream ss(line);
            ss >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz >> vx >> vy >> vz;
            id--;  // the atom id in dump file starts from 1
            system.atoms[id].x = x;
            system.atoms[id].y = y;
            system.atoms[id].z = z;
            system.atoms[id].mol = mol;
            system.atoms[id].type = type;
            system.atoms[id].ix = ix;
            system.atoms[id].iy = iy;
            system.atoms[id].iz = iz;
            system.atoms[id].vx = vx;
            system.atoms[id].vy = vy;
            system.atoms[id].vz = vz;
            // if (system.atoms[id].id != id + 1) iderr++;
        }
        // if (iderr != 0) cout << "ID error: " << iderr <<", frame:"<<f<< endl;
        unwrap(system);
        // calculate the P2
        double p2x = p2calculation(system);
        p2out << f << " " << p2x << endl;
        cout << "P2 = " << p2x << endl;
        p2out.flush();
    }
}

void Profile_velocity(System& system, ifstream& dumpFilein, double zlol, double zhil, int nslice)
{
    string line;
    istringstream ss(line);
    ofstream p2out("P2.txt");
    long int f = 0;
    vector <double> zslicev(nslice, 0);
    vector <int> sn(nslice, 0);
    int k;
    p2out << "TIMESTEP" << "      ";
    for(int i=0; i<nslice; i++) p2out << "SLICE_" << i << "      ";
    p2out << endl;
    while (f < system.frames)
    {
        // long int iderr = 0;
        auto fstart = std::chrono::high_resolution_clock::now();
        f++;
        cout << endl << "Calculating frame " << f << "/" << system.frames << " ..." << endl;
        long int id, flable;
        int mol, type;
        double x, y, z, vx, vy, vz;
        long int ix, iy, iz;
        double slicelen;  // the length of each slice
        slicelen = (double)(zhil - zlol) / nslice;

        // update the atom position and velocity
        while (line != "ITEM: ATOMS id mol type x y z ix iy iz vx vy vz")
        {
            getline(dumpFilein, line);
            // cout << line << endl;
        }
        for (long int i = 0; i < system.num_atoms; i++)
        {
            getline(dumpFilein, line);
            istringstream ss(line);
            ss >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz >> vx >> vy >> vz;
            id--;  // the atom id in dump file starts from 1
            // system.atoms[id].x = x;
            // system.atoms[id].y = y;
            // system.atoms[id].z = z;
            // system.atoms[id].mol = mol;
            // system.atoms[id].type = type;
            // system.atoms[id].ix = ix;
            // system.atoms[id].iy = iy;
            // system.atoms[id].iz = iz;
            // system.atoms[id].vx = vx;
            // system.atoms[id].vy = vy;
            // system.atoms[id].vz = vz;
            if (z >= zlol && z <= zhil)
            {
                k = floor((z - zlol) / slicelen);
                zslicev[k] += vy;
                sn[k]++;
            }
            
        }
        for(int c=0; c < nslice; c++)
        {
            zslicev[c] /= sn[c];
        }
        p2out << f <<"      ";
        for(int c=0; c < nslice; c++)
        {
            p2out << zslicev[c] <<"      ";
        }
        p2out << endl;   
        p2out.flush();
    }
}
