#pragma once
#include <vector>
#include "data_types.h"
using namespace std;
struct System{
    int num_atom_types, num_bond_types, num_angle_types,num_dihedral_types, num_improper_types;
    long int num_atoms = 0, num_bonds = 0, num_angles = 0, num_dihedrals = 0, num_impropers = 0;
    vector<Atom> atoms;
    vector<Bond> bonds;
    vector<Angle> angles;
    vector<Dihedral> dihedrals;
    vector<Improper> impropers;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double bx, by, bz;
    int box_type = 0; // 0: orthogonal, 1: triclinic
    double xtilt, ytilt, ztilt; // triclinic box tilt factors
    vector<double> mass;
    long int molnumMAX = -1;
    vector<long int> molen;
    long int frames = 0;
};