// find the position vector of two atoms seperated by 6 bonds
// define the orientation vector ovec = (0, 1, 0)
// calculate the angle between vector and ovec
// calculate the order parameter by Second-Order Legendre Polynomial
// return the order parameter
# include "../include/p2_calculation.h"
using namespace std;

void printProgressP2(int k, int kmax, double p2) {
    int barWidth = 50;  // Width of the progress bar
    float progress = static_cast<float>(k) / kmax;
    int pos = static_cast<int>(barWidth * progress);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " % (" << k << "/" << kmax << ", P2 = " << p2 << ")\r";
    std::cout.flush();
}

double p2calculation(System& system, long int f)
{
    axis r1;
    long int p2num = 0;
    double p2sum = 0;
    for (long int i = 0; i < system.num_atoms - 9; i++)
    {
        if(system.atoms[i].mol == system.atoms[i + 9].mol && ((system.atoms[i].type == 2 || system.atoms[i].type == 1) && (system.atoms[i + 9].type == 2 || system.atoms[i + 9].type == 1)))
        {
            // cout << i+1 << "  ";
            r1.x = system.atoms[i + 9].x - system.atoms[i].x;
            r1.y = system.atoms[i + 9].y - system.atoms[i].y;
            r1.z = system.atoms[i + 9].z - system.atoms[i].z;
            axis ovec;
            ovec.x = 0;
            ovec.y = 1;
            ovec.z = 0;
            double cos_theta = dot(r1, ovec) / norm(r1);
            // cout << "cos_theta: " << cos_theta << endl;
            p2sum += (0.5 * (3 * cos_theta * cos_theta - 1));
            p2num++;
        }
    }
    double p2 = p2sum / p2num;
    printProgressP2(f, system.frames, p2);
    return p2;
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
        double p2x = p2calculation(system, f);
        p2out << f << " " << p2x << endl;
        p2out.flush();
    }
}

