#include "../include/dihedral_calculation.h"
#include <algorithm>

using namespace std;

namespace {

struct BackboneDihedral {
    long int dihedralId;
    long int atom1;
    long int atom2;
    long int atom3;
    long int atom4;
};

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

bool isBackbonePattern(int t1, int t2, int t3, int t4)
{
    // target patterns for iPP backbone:
    // begin : 3-1-2-1
    // middle: 1-2-1-2 or 2-1-2-1
    // end   : 2-1-2-2
    // NOTE:
    // Keep begin/end in the exact direction above. Including reversed begin/end
    // (1-2-1-3, 2-2-1-2) overcounts by also capturing non-backbone dihedrals.
    if (t1 == 3 && t2 == 1 && t3 == 2 && t4 == 1) return true;
    if ((t1 == 1 && t2 == 2 && t3 == 1 && t4 == 2) || (t1 == 2 && t2 == 1 && t3 == 2 && t4 == 1)) return true;
    if (t1 == 2 && t2 == 1 && t3 == 2 && t4 == 2) return true;
    return false;
}

bool isStrictIncreasingAtomId(long int a1, long int a2, long int a3, long int a4)
{
    // User rule: atom ids along iPP backbone are strictly incremental.
    return (a1 < a2) && (a2 < a3) && (a3 < a4);
}

double calcDihedralDeg(const axis& r1, const axis& r2, const axis& r3, const axis& r4)
{
    axis b1{r2.x - r1.x, r2.y - r1.y, r2.z - r1.z, 0};
    axis b2{r3.x - r2.x, r3.y - r2.y, r3.z - r2.z, 0};
    axis b3{r4.x - r3.x, r4.y - r3.y, r4.z - r3.z, 0};

    axis n1, n2;
    cross(n1, b1, b2);
    cross(n2, b2, b3);

    double n1n = norm(n1);
    double n2n = norm(n2);
    double b2n = norm(b2);
    if (n1n < 1e-14 || n2n < 1e-14 || b2n < 1e-14) return 0.0;

    axis b2u{b2.x / b2n, b2.y / b2n, b2.z / b2n, 0};
    axis m1;
    cross(m1, n1, b2u);

    double x = dot(n1, n2);
    double y = dot(m1, n2);
    return atan2(y, x) * 180.0 / M_PI;
}

char classifyRIS(double deg)
{
    // T: (-pi, -5/6*pi] U [5/6*pi, pi]  => (-180, -150] U [150, 180]
    // g: [-pi/2, -pi/6]                  => [-90, -30]
    // G: [ pi/6,  pi/2]                  => [ 30,  90]
    if ((deg > -180.0 && deg <= -150.0) || (deg >= 150.0 && deg <= 180.0)) return 'T';
    if (deg >= -90.0 && deg <= -30.0) return 'g';
    if (deg >= 30.0 && deg <= 90.0) return 'G';
    return 'X'; // other
}

} // namespace

void dumpIO_backboneDihedral(System& system, ifstream& dumpFilein)
{
    // Set to false to turn off detailed per-chain outputs.
    const bool outputBackboneDetails = false;

    map<int, vector<BackboneDihedral>> backboneByMol;

    for (const auto& dih : system.dihedrals)
    {
        if (dih.atom1 <= 0 || dih.atom2 <= 0 || dih.atom3 <= 0 || dih.atom4 <= 0) continue;
        if (dih.atom1 > system.num_atoms || dih.atom2 > system.num_atoms || dih.atom3 > system.num_atoms || dih.atom4 > system.num_atoms) continue;

        const Atom& a1 = system.atoms[dih.atom1 - 1];
        const Atom& a2 = system.atoms[dih.atom2 - 1];
        const Atom& a3 = system.atoms[dih.atom3 - 1];
        const Atom& a4 = system.atoms[dih.atom4 - 1];

        if (!(a1.mol == a2.mol && a2.mol == a3.mol && a3.mol == a4.mol)) continue;
        if (!isStrictIncreasingAtomId(dih.atom1, dih.atom2, dih.atom3, dih.atom4)) continue;
        if (!isBackbonePattern(a1.type, a2.type, a3.type, a4.type)) continue;

        BackboneDihedral item{dih.id, dih.atom1, dih.atom2, dih.atom3, dih.atom4};
        backboneByMol[a1.mol].push_back(item);
    }

    for (auto& kv : backboneByMol)
    {
        sort(kv.second.begin(), kv.second.end(), [](const BackboneDihedral& lhs, const BackboneDihedral& rhs) {
            return lhs.dihedralId < rhs.dihedralId;
        });
    }

    ofstream dout;
    ofstream risout;
    if (outputBackboneDetails)
    {
        dout.open("backbone_dihedral_degree.txt");
        dout << "# frame timestep mol dihedral_count dihedral_deg_1 ... dihedral_deg_n" << endl;
        risout.open("backbone_dihedral_RIS.txt");
        risout << "# frame timestep mol dihedral_count RIS_1 ... RIS_n (T/g/G/X)" << endl;
    }
    ofstream nturnout("Nturns.txt");
    nturnout << "# frame timestep <Nturns>_InterpretationA" << endl;

    long int f = 0;
    while (f < system.frames)
    {
        long long timestep = 0;
        if (!readDumpFrame(system, dumpFilein, timestep)) break;
        f++;

        double sumChainNturns = 0.0;
        long int chainCount = 0;

        for (const auto& kv : backboneByMol)
        {
            int mol = kv.first;
            const auto& list = kv.second;
            vector<char> risList;
            risList.reserve(list.size());

            if (outputBackboneDetails)
            {
                dout << f << " " << timestep << " " << mol << " " << list.size();
                risout << f << " " << timestep << " " << mol << " " << list.size();
            }
            for (const auto& d : list)
            {
                const Atom& a1 = system.atoms[d.atom1 - 1];
                const Atom& a2 = system.atoms[d.atom2 - 1];
                const Atom& a3 = system.atoms[d.atom3 - 1];
                const Atom& a4 = system.atoms[d.atom4 - 1];
                axis r1{a1.x, a1.y, a1.z, 0};
                axis r2{a2.x, a2.y, a2.z, 0};
                axis r3{a3.x, a3.y, a3.z, 0};
                axis r4{a4.x, a4.y, a4.z, 0};
                double deg = calcDihedralDeg(r1, r2, r3, r4);
                if (outputBackboneDetails) dout << " " << deg;
                char ris = classifyRIS(deg);
                if (outputBackboneDetails) risout << " " << ris;
                risList.push_back(ris);
            }
            if (outputBackboneDetails)
            {
                dout << endl;
                risout << endl;
            }

            // <Nturns> (Interpretation A) for one chain in this frame.
            // Algorithm notes:
            // 1) Backbone dihedrals per chain are 2N-3, so N = (nDihed+3)/2.
            // 2) Build NON-overlapping adjacent-dihedral pairs:
            //      (d1,d2), (d3,d4), ...
            //    This gives (2N-4)/2 = N-2 pairs, matching the definition.
            // 3) Tag pair as RH if "Tg", LH if "TG", otherwise non-helical.
            // 4) Group consecutive same-tag pairs into helical groups.
            //    If group lengths are m_k, then:
            //      f(N) = N / (3*(N-2))
            //      Nturns(chain) = [f(N)/nGrp] * sum_k(m_k)
            // 5) Frame output is average over chains.
            long int nDihed = static_cast<long int>(risList.size());
            long int Nmono = (nDihed + 3) / 2;
            double fconv = 0.0;
            if (Nmono > 2)
            {
                fconv = static_cast<double>(Nmono) / (3.0 * static_cast<double>(Nmono - 2));
            }

            // Build NON-overlapping adjacent-dihedral pairs:
            // (d1,d2), (d3,d4), ..., total pairs = (2N-4)/2 = N-2
            long int nPair = (nDihed - 1) / 2;
            vector<int> pairType(nPair, 0); // 0: other, 1: Tg, 2: TG
            for (long int p = 0; p < nPair; p++)
            {
                long int i0 = 2 * p;
                long int i1 = i0 + 1;
                if (risList[i0] == 'T' && risList[i1] == 'g') pairType[p] = 1;
                else if (risList[i0] == 'T' && risList[i1] == 'G') pairType[p] = 2;
            }

            // Group consecutive Tg or TG pairs (RH/LH separately), then combine all groups.
            long int nGrp = 0;
            long int sumM = 0;
            long int p = 0;
            while (p < nPair)
            {
                if (pairType[p] == 0)
                {
                    p++;
                    continue;
                }

                int t = pairType[p];
                long int m = 1;
                p++;
                while (p < nPair && pairType[p] == t)
                {
                    m++;
                    p++;
                }
                nGrp++;
                sumM += m;
            }

            double chainNturns = 0.0;
            if (nGrp > 0) chainNturns = fconv * static_cast<double>(sumM) / static_cast<double>(nGrp);
            sumChainNturns += chainNturns;
            chainCount++;
        }

        double avgNturns = 0.0;
        if (chainCount > 0) avgNturns = sumChainNturns / static_cast<double>(chainCount);
        nturnout << f << " " << timestep << " " << avgNturns << endl;

        printProgressBar(f, system.frames);
    }

    if (outputBackboneDetails)
    {
        dout.close();
        risout.close();
    }
    nturnout.close();
}
