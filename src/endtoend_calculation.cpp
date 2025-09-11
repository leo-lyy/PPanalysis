# include "../include/endtoend_calculation.h"

using namespace std;

void reeCalculation(System& system, long int f, ofstream& Reeout)
{
    // define the Ree as a vector of the end to end distance
    axis Ree[system.molnumMAX + 1];
    axis Ree_sum;
    Ree_sum.x = 0;
    Ree_sum.y = 0;
    Ree_sum.z = 0;
    long int head_id = 0;                       // for 0-based index
    long int tail_id = system.molen[1] - 3;
    for (long int i = 1; i < system.molnumMAX + 1; i++)
    {
        // Reeout << f << " " << head_id + 1 << " " << tail_id + 1 << endl;
        // calculate the mass center of the end monomer
        axis Ree_head, Ree_tail;
        double mass_head = system.mass[system.atoms[head_id].type - 1] + system.mass[system.atoms[head_id + 1].type -1] + system.mass[system.atoms[head_id + 2].type -1];
        double mass_tail = system.mass[system.atoms[tail_id].type - 1] + system.mass[system.atoms[tail_id + 1].type -1] + system.mass[system.atoms[tail_id + 2].type -1];
        // Reeout << mass_head << " " << mass_tail << endl;
        Ree_head.x = (system.atoms[head_id].x * system.mass[system.atoms[head_id].type - 1] + system.atoms[head_id + 1].x * system.mass[system.atoms[head_id + 1].type -1]  + system.atoms[head_id + 2].x * system.mass[system.atoms[head_id + 2].type -1] )/ mass_head;
        Ree_head.y = (system.atoms[head_id].y * system.mass[system.atoms[head_id].type - 1] + system.atoms[head_id + 1].y * system.mass[system.atoms[head_id + 1].type -1]  + system.atoms[head_id + 2].y * system.mass[system.atoms[head_id + 2].type -1] )/ mass_head;
        Ree_head.z = (system.atoms[head_id].z * system.mass[system.atoms[head_id].type - 1] + system.atoms[head_id + 1].z * system.mass[system.atoms[head_id + 1].type -1]  + system.atoms[head_id + 2].z * system.mass[system.atoms[head_id + 2].type -1] )/ mass_head;
        Ree_tail.x = (system.atoms[tail_id].x * system.mass[system.atoms[tail_id].type - 1] + system.atoms[tail_id + 1].x * system.mass[system.atoms[tail_id + 1].type -1]  + system.atoms[tail_id + 2].x * system.mass[system.atoms[tail_id + 2].type -1] )/ mass_tail;
        Ree_tail.y = (system.atoms[tail_id].y * system.mass[system.atoms[tail_id].type - 1] + system.atoms[tail_id + 1].y * system.mass[system.atoms[tail_id + 1].type -1]  + system.atoms[tail_id + 2].y * system.mass[system.atoms[tail_id + 2].type -1] )/ mass_tail;
        Ree_tail.z = (system.atoms[tail_id].z * system.mass[system.atoms[tail_id].type - 1] + system.atoms[tail_id + 1].z * system.mass[system.atoms[tail_id + 1].type -1]  + system.atoms[tail_id + 2].z * system.mass[system.atoms[tail_id + 2].type -1] )/ mass_tail;
        Ree[i].x = Ree_tail.x - Ree_head.x;
        Ree[i].y = Ree_tail.y - Ree_head.y;
        Ree[i].z = Ree_tail.z - Ree_head.z;
        Ree[i].id = i;
        Ree_sum.x += Ree[i].x;
        Ree_sum.y += Ree[i].y;
        Ree_sum.z += Ree[i].z;
        if (i + 1 < system.molnumMAX + 1)
        {
            head_id += system.molen[i];
            tail_id = head_id + system.molen[i+1] - 3;
        }
    }
    Ree_sum.x /= system.molnumMAX;
    Ree_sum.y /= system.molnumMAX;
    Ree_sum.z /= system.molnumMAX;
    Reeout << f << " " << Ree_sum.x << " " << Ree_sum.y << " " << Ree_sum.z << endl;
}

void dumpIO_Ree(System& system, ifstream& dumpFilein)
{
    // dump file format: ITEM: ATOMS id mol type x y z ix iy iz
    string line;
    istringstream ss(line);
    ofstream Reeout("Ree.txt");
    long int f = 0;
    while (f < system.frames)
    {
        // long int iderr = 0;
        // auto fstart = std::chrono::high_resolution_clock::now();
        printProgressBar(f+1, system.frames);
        f++;
        long int id, flable;
        int mol, type;
        double x, y, z, vx, vy, vz;
        long int ix, iy, iz;
        // update the atom position and velocity
        while (getline(dumpFilein, line))
        {
            if (line.find("ITEM: ATOMS") != string::npos)
            {
                break;
            }
        }
        bool hasVel = line.find("vx") != string::npos;
        for (long int i = 0; i < system.num_atoms; i++)
        {
            getline(dumpFilein, line);
            istringstream ss(line);
            if (hasVel)
            {
                ss >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz >> vx >> vy >> vz;
            }
            else
            {
                ss >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz;
            }
            id--;  // the atom id in dump file starts from 1
            system.atoms[id].x = x;
            system.atoms[id].y = y;
            system.atoms[id].z = z;
            system.atoms[id].mol = mol;
            system.atoms[id].type = type;
            system.atoms[id].ix = ix;
            system.atoms[id].iy = iy;
            system.atoms[id].iz = iz;
            // if (system.atoms[id].id != id + 1) iderr++;
        }
        unwrap(system);

        // calculate the Ree
        reeCalculation(system, f, Reeout);
    }
}
