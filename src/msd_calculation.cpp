#include "../include/msd_calculation.h"

using namespace std;

namespace {
bool readDumpFrame(System& system, ifstream& dumpFilein, long long& timestep)
{
    string line;
    bool gotTimestep = false;
    bool gotBox = false;
    bool gotAtoms = false;

    while (getline(dumpFilein, line))
    {
        if (line.rfind("ITEM: TIMESTEP", 0) == 0)
        {
            if (!getline(dumpFilein, line)) return false;
            istringstream ss(line);
            ss >> timestep;
            gotTimestep = true;
        }
        else if (line.rfind("ITEM: NUMBER OF ATOMS", 0) == 0)
        {
            // skip one line containing atom count
            if (!getline(dumpFilein, line)) return false;
        }
        else if (line.rfind("ITEM: BOX BOUNDS", 0) == 0)
        {
            bool triclinic = (line.find("xy") != string::npos);
            system.box_type = triclinic ? 1 : 0;

            if (!getline(dumpFilein, line)) return false;
            {
                istringstream ss(line);
                if (triclinic) ss >> system.xlo >> system.xhi >> system.xtilt;
                else ss >> system.xlo >> system.xhi;
            }
            if (!getline(dumpFilein, line)) return false;
            {
                istringstream ss(line);
                if (triclinic) ss >> system.ylo >> system.yhi >> system.ytilt;
                else ss >> system.ylo >> system.yhi;
            }
            if (!getline(dumpFilein, line)) return false;
            {
                istringstream ss(line);
                if (triclinic) ss >> system.zlo >> system.zhi >> system.ztilt;
                else ss >> system.zlo >> system.zhi;
            }

            system.bx = system.xhi - system.xlo;
            system.by = system.yhi - system.ylo;
            system.bz = system.zhi - system.zlo;
            gotBox = true;
        }
        else if (line.rfind("ITEM: ATOMS", 0) == 0)
        {
            bool hasVel = (line.find("vx") != string::npos);
            for (long int i = 0; i < system.num_atoms; i++)
            {
                long int id;
                int mol, type;
                long int ix, iy, iz;
                double x, y, z, vx = 0.0, vy = 0.0, vz = 0.0;

                if (hasVel)
                {
                    if (!(dumpFilein >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz >> vx >> vy >> vz)) return false;
                }
                else
                {
                    if (!(dumpFilein >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz)) return false;
                }

                if (id <= 0 || id > system.num_atoms) return false;
                size_t idx = static_cast<size_t>(id - 1);
                system.atoms[idx].id = id;
                system.atoms[idx].mol = mol;
                system.atoms[idx].type = type;
                system.atoms[idx].x = x;
                system.atoms[idx].y = y;
                system.atoms[idx].z = z;
                system.atoms[idx].ix = ix;
                system.atoms[idx].iy = iy;
                system.atoms[idx].iz = iz;
                system.atoms[idx].vx = vx;
                system.atoms[idx].vy = vy;
                system.atoms[idx].vz = vz;
            }

            // clear the trailing newline left by operator>>
            getline(dumpFilein, line);

            gotAtoms = true;
            break;
        }
    }

    if (!(gotTimestep && gotBox && gotAtoms)) return false;
    if (system.bx <= 0.0 || system.by <= 0.0 || system.bz <= 0.0) return false;

    unwrap(system);
    return true;
}
} // namespace

void dumpIO_MSD(System& system, ifstream& dumpFilein)
{
    ofstream msdout("msd.txt");
    msdout << "# frame timestep msd_x msd_y msd_z msd_total" << endl;

    vector<axis> ref(system.num_atoms);

    long int f = 0;
    while (f < system.frames)
    {
        long long timestep = 0;
        if (!readDumpFrame(system, dumpFilein, timestep)) break;
        f++;

        if (f == 1)
        {
            for (long int i = 0; i < system.num_atoms; i++)
            {
                ref[i].x = system.atoms[i].x;
                ref[i].y = system.atoms[i].y;
                ref[i].z = system.atoms[i].z;
            }
            msdout << f << " " << timestep << " "
                   << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
            printProgressBar(f, system.frames);
            continue;
        }

        double msdX = 0.0, msdY = 0.0, msdZ = 0.0;
        for (long int i = 0; i < system.num_atoms; i++)
        {
            double dx = system.atoms[i].x - ref[i].x;
            double dy = system.atoms[i].y - ref[i].y;
            double dz = system.atoms[i].z - ref[i].z;
            msdX += dx * dx;
            msdY += dy * dy;
            msdZ += dz * dz;
        }

        msdX /= system.num_atoms;
        msdY /= system.num_atoms;
        msdZ /= system.num_atoms;
        double msd = msdX + msdY + msdZ;

        msdout << f << " " << timestep << " "
               << msdX << " " << msdY << " " << msdZ << " " << msd << endl;
        printProgressBar(f, system.frames);
    }

    msdout.close();
}
