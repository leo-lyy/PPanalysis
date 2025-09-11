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
            axis ovec{0, 1, 0, 0};
            double cos_theta = dot(r1, ovec) / norm(r1);
            // cout << "cos_theta: " << cos_theta << endl;
            p2sum += (0.5 * (3 * cos_theta * cos_theta - 1));
            p2num++;
        }
    }
    double p2 = p2num > 0 ? p2sum / p2num : 0.0;
    printProgressP2(f, system.frames, p2);
    return p2;
}

void dumpIO_p2(System& system, ifstream& dumpFilein)
{
    string line;
    ofstream p2out("P2.txt");
    long int f = 0;
    while (f < system.frames)
    {
        // long int iderr = 0;
        auto fstart = std::chrono::high_resolution_clock::now();
        f++;
        long long timestep = -1;
        bool gotAtoms = false;
        bool gotBox = false;
        bool gotTimestep = false;
        string header;
        while (std::getline(dumpFilein, header))
        {
            if (header.rfind("ITEM: TIMESTEP", 0) == 0) {
                if (!std::getline(dumpFilein, line)) return;
                { istringstream iss(line); iss >> timestep; }
                gotTimestep = true;
            }
            else if (header.rfind("ITEM: NUMBER OF ATOMS", 0) == 0) {
                if (!std::getline(dumpFilein, line)) return;
            }
            else if (header.rfind("ITEM: BOX BOUNDS", 0) == 0) {
                bool triclinic = header.find("xy") != string::npos;
                system.box_type = triclinic ? 1 : 0;
                if (!getline(dumpFilein, line)) return;
                if (!getline(dumpFilein, line)) return;
                if (!getline(dumpFilein, line)) return;
                gotBox = true;
            }
            else if (header.rfind("ITEM: ATOMS", 0) == 0) {
                bool hasVel = header.find("vx") != string::npos;
                if (!hasVel) {
                    cerr << "Error: P2 requires velocities in dump." << endl;
                    return;
                }
                for (long int i = 0; i < system.num_atoms; i++) {
                    long long id;
                    int mol, type;
                    double x, y, z, vx, vy, vz;
                    int ix, iy, iz;
                    if (!(dumpFilein >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz >> vx >> vy >> vz)) return;
                    if (id <= 0 || id > system.num_atoms) return;
                    size_t idx = size_t(id - 1);
                    auto& a = system.atoms[idx];
                    a.mol = mol;
                    a.type = type;
                    a.x = x;
                    a.y = y;
                    a.z = z;
                    a.ix = ix;
                    a.iy = iy;
                    a.iz = iz;
                    a.vx = vx;
                    a.vy = vy;
                    a.vz = vz;
                }
                gotAtoms = true;
                std::getline(dumpFilein, line);
                break;
            }
            if (gotAtoms) break;
        }
        if (!(gotAtoms && gotTimestep)) break;
        unwrap(system);
        double p2x = p2calculation(system, f);
        p2out << f << " " << p2x << "\n";
        p2out.flush();
    }
}

