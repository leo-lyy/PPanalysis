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

void dumpIO_p2_slice(System& system, ifstream& dumpFilein, double zlo, double zhi, int nslice)
{
    string line;
    ofstream p2out_slice("P2_slice.txt");
    ofstream p2out_global("P2.txt");
    long int f = 0;
    double slice_width = (zhi - zlo) / nslice;
    
    while (f < system.frames)
    {
        f++;
        // Initialize slice statistics for this frame
        vector<double> p2sum_slice(nslice, 0.0);
        vector<long int> p2num_slice(nslice, 0);
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
        
        // Calculate P2 for each slice
        axis r1;
        for (long int i = 0; i < system.num_atoms - 9; i++)
        {
            if(system.atoms[i].mol == system.atoms[i + 9].mol && 
               ((system.atoms[i].type == 2 || system.atoms[i].type == 1) && 
                (system.atoms[i + 9].type == 2 || system.atoms[i + 9].type == 1)))
            {
                r1.x = system.atoms[i + 9].x - system.atoms[i].x;
                r1.y = system.atoms[i + 9].y - system.atoms[i].y;
                r1.z = system.atoms[i + 9].z - system.atoms[i].z;
                
                // Calculate midpoint z coordinate
                double z_mid = (system.atoms[i].z + system.atoms[i + 9].z) / 2.0;
                
                // Determine which slice this vector belongs to
                int slice_idx = static_cast<int>((z_mid - zlo) / slice_width);
                
                // Check if within bounds
                if (slice_idx >= 0 && slice_idx < nslice)
                {
                    axis ovec{0, 1, 0, 0};
                    double cos_theta = dot(r1, ovec) / norm(r1);
                    double p2_value = 0.5 * (3 * cos_theta * cos_theta - 1);
                    p2sum_slice[slice_idx] += p2_value;
                    p2num_slice[slice_idx]++;
                }
            }
        }
        
        // Output results for this frame
        if (f == 1) {
            // First line: z coordinates
            for (int s = 0; s < nslice; s++)
            {
                double z_center = zlo + (s + 0.5) * slice_width;
                p2out_slice << z_center;
                if (s < nslice - 1) p2out_slice << " ";
            }
            p2out_slice << "\n";
        }
        // Each frame: P2 values
        for (int s = 0; s < nslice; s++)
        {
            double p2_avg = p2num_slice[s] > 0 ? p2sum_slice[s] / p2num_slice[s] : 0.0;
            p2out_slice << p2_avg;
            if (s < nslice - 1) p2out_slice << " ";
        }
        p2out_slice << "\n";
        p2out_slice.flush();
        
        // Calculate and output global P2 for this frame
        double avg_p2 = 0.0;
        long int total_num = 0;
        for (int s = 0; s < nslice; s++) {
            avg_p2 += p2sum_slice[s];
            total_num += p2num_slice[s];
        }
        avg_p2 = total_num > 0 ? avg_p2 / total_num : 0.0;
        p2out_global << f << " " << avg_p2 << "\n";
        p2out_global.flush();
        
        // Print progress
        printProgressP2(f, system.frames, avg_p2);
    }
    std::cout << std::endl;  // New line after progress bar completes
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

