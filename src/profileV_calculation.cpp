
#include"../include/profileV_calculation.h"

using namespace std;

void Profile_velocity(System& system, ifstream& dumpFilein, double zlol, double zhil, int nslice)
{
    string line;
    istringstream ss(line);
    ofstream pvout("Profile_velocity.txt");
    long int f = 0;
    vector <double> zslicev(nslice, 0);
    vector <int> sn(nslice, 0);
    double slicelen;  // the length of each slice
    slicelen = (double)(zhil - zlol) / nslice;
    int k;
    pvout << "TIMESTEP" << "      ";
    for(int i=0; i<nslice; i++) pvout << "SLICE_" << i << "      ";
    pvout << endl;
    while (f < system.frames)
    {
        // long int iderr = 0;
        // auto fstart = std::chrono::high_resolution_clock::now();
        f++;
        printProgressBar(f, system.frames);
        long int id, flable;
        int mol, type;
        double x, y, z, vx, vy, vz;
        long int ix, iy, iz;
        fill(zslicev.begin(), zslicev.end(), 0);
        fill(sn.begin(), sn.end(), 0);

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
        pvout << f <<"      ";
        for(int c=0; c < nslice; c++)
        {
            pvout << zslicev[c] <<"      ";
        }
        pvout << endl;   
        pvout.flush();
    }
}
