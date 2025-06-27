#include"../include/readLMPdata.h"

using namespace std;
bool compareAtomID(const Atom& a, const Atom& b){return a.id < b.id;}
void readGeneralInfo(ifstream& file, System& system)
{
    string line;
    getline(file, line);
    getline(file, line);
    while (getline(file, line))
    {
        istringstream ss(line);
        if (line.find("atoms") != string::npos)         ss >> system.num_atoms;
        // cout <<"!!!!!"<< system.num_atoms << endl;
        if (line.find("atom types") != string::npos)    ss >> system.num_atom_types;
        if (line.find("bonds") != string::npos)         ss >> system.num_bonds;
        if (line.find("bond types") != string::npos)    ss >> system.num_bond_types;
        if (line.find("angles") != string::npos)        ss >> system.num_angles;
        if (line.find("angle types") != string::npos)   ss >> system.num_angle_types;
        if (line.find("dihedrals") != string::npos)     ss >> system.num_dihedrals;
        if (line.find("dihedral types") != string::npos)ss >> system.num_dihedral_types;
        if (line.find("impropers") != string::npos)     ss >> system.num_impropers;
        if (line.find("improper types") != string::npos)ss >> system.num_improper_types;
        if (line.empty()) break;
    }
}
void readBox(ifstream& file, System& system)
{
    string line;
    while (getline(file, line))
    {
        if (line.find("xlo") != string::npos)
        {
            istringstream ss(line);
            ss >> system.xlo >> system.xhi;
            system.bx = system.xhi - system.xlo;
        }
        else if (line.find("ylo") != string::npos)
        {
            istringstream ss(line);
            ss >> system.ylo >> system.yhi;
            system.by = system.yhi - system.ylo;
        }
        else if (line.find("zlo") != string::npos)
        {
            istringstream ss(line);
            ss >> system.zlo >> system.zhi;
            system.bz = system.zhi - system.zlo;
            break;
        }
    }
}
void readMasses(ifstream& file, System& system)
{
    string line;
    while (getline(file, line))
    {
        if (line == "Masses")
        {
            getline(file, line);
            for (int i = 0; i < system.num_atom_types; i++)
            {
                getline(file, line);
                istringstream ss(line);
                double m;
                int id;
                ss >> id >> m;
                system.mass.push_back(m);
            }
            break;
        }
    }
}
void readAtoms(ifstream& file, System& system)
{
    string line;
    while (getline(file, line))
    {
        if (line == "Atoms # molecular")
        {
            getline(file, line);
            for (long int i = 0; i < system.num_atoms; i++)
            {
                getline(file, line);
                istringstream ss(line);
                Atom atom;
                ss >> atom.id >> atom.mol >> atom.type >> atom.x >> atom.y >> atom.z >> atom.ix >> atom.iy >> atom.iz;
                system.atoms.push_back(atom);
                if (system.molnumMAX < atom.mol) system.molnumMAX = atom.mol;
            }
            break;
        }
    }
}
void readVelocities(ifstream& file, System& system)
{
    string line;
    long int vid;
    while (getline(file, line))
    {
        if (line == "Velocities")
        {
            getline(file, line);
            for (long int i = 0; i < system.num_atoms; i++)
            {
                getline(file, line);
                istringstream ss(line);
                ss >> vid;
                vid--;
                ss >> system.atoms[vid].vx >> system.atoms[vid].vy >> system.atoms[vid].vz;
            }
            break;
        }
    }
}
void readBonds(ifstream& file, System& system)
{
    string line;
    
    while (getline(file, line))
    {
        if (line == "Bonds")
        {
            getline(file, line);
            for (long int i = 0; i < system.num_bonds; i++)
            {
                getline(file, line);
                istringstream ss(line);
                Bond bond;
                ss >> bond.id >> bond.type >> bond.atom1 >> bond.atom2;
                system.bonds.push_back(bond);
            }
            break;
        }
    }
}
void readAngles(ifstream& file, System& system)
{
    string line;
    while (getline(file, line))
    {
        if (line == "Angles")
        {
            getline(file, line);
            for (long int i = 0; i < system.num_angles; i++)
            {
                getline(file, line);
                istringstream ss(line);
                Angle angle;
                ss >> angle.id >> angle.type >> angle.atom1 >> angle.atom2 >> angle.atom3;
                system.angles.push_back(angle);
            }
            break;
        }
    }

}
void readDihedrals(ifstream& file, System& system)
{
    string line;
    while (getline(file, line))
    {
        if (line == "Dihedrals")
        {
            getline(file, line);
            for (long int i = 0; i < system.num_dihedrals; i++)
            {
                getline(file, line);
                istringstream ss(line);
                Dihedral dihedral;
                ss >> dihedral.id >> dihedral.type >> dihedral.atom1 >> dihedral.atom2 >> dihedral.atom3 >> dihedral.atom4;
                system.dihedrals.push_back(dihedral);
            }
            break;
        }
    }
}
void readImpropers(ifstream& file, System& system)
{
    string line;
    while (getline(file, line))
    {
        if (line == "Impropers")
        {
            getline(file, line);
            for (long int i = 0; i < system.num_impropers; i++)
            {
                getline(file, line);
                istringstream ss(line);
                Improper improper;
                ss >> improper.id >> improper.type >> improper.atom1 >> improper.atom2 >> improper.atom3 >> improper.atom4;
                system.impropers.push_back(improper);
            }
            break;
        }
    }
}

void countmolen(System& system)
{
    system.molen.resize(system.molnumMAX+1);
    fill(system.molen.begin(), system.molen.end(), 0);
    for (long int i = 0; i < system.num_atoms; i++)system.molen[system.atoms[i].mol]++;
}
bool fileExists(const std::string& filename) {
    return std::filesystem::exists(filename);
}
void readLammpsData(const string& dataFileName, System& system)
{
    // check if the file exists
    if (!fileExists(dataFileName)) 
    {
        cout << "data file does not exist!" << endl;
    }

    auto start = std::chrono::high_resolution_clock::now();// 记录读取开始的时间
    ifstream file(dataFileName);
    string line;

    // read the head part
    readGeneralInfo(file, system);

    // read the box info
    readBox(file, system);
    
    // read the mass info
    readMasses(file, system);

    // read the atom info
    cout << "Reading atom info ..." << endl;
    readAtoms(file, system);
    cout << "Done! Total atom number: " << system.atoms.size() << endl;
    // sort the system.atoms by system.atoms.id
    sort(system.atoms.begin(), system.atoms.end(), compareAtomID);

    // read the velocity
    cout << "Reading velocity info ..." << endl;
    readVelocities(file, system);
    cout << "Done! Toteal atom velocity number: "<< system.atoms.size() << endl;
    // read the bonds
    cout << "Reading bond info ..." << endl;
    readBonds(file, system);
    cout << "Done! Total bond number: " << system.bonds.size() << endl;
    // read the angles
    cout << "Reading angle info ..." << endl;
    readAngles(file, system);
    cout << "Done! Total angle number: " << system.angles.size() << endl;
    // read the dihedrals
    cout << "Reaning dihedral info ..." << endl;
    readDihedrals(file, system);
    cout << "Done! Total dihedral number: " << system.dihedrals.size() << endl;
    // read the impropers
    cout << "Reading improper info ..." << endl;
    readImpropers(file, system);
    cout << "Done! Total improper number: " << system.impropers.size() << endl;

    file.close();
    auto end = std::chrono::high_resolution_clock::now();// 记录读取结束的时间
    std::chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time for reading Data file: " << elapsed.count() << " s" << endl;
    
    unwrap(system);    // unwrap the box for data file
    countmolen(system);

}