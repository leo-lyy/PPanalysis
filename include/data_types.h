#pragma once
#include <vector>
struct axis{
    double x;
    double y;
    double z;
    long int id;
};
struct Atom{
    long int id;
    int mol;
    int type;
    double x, y, z;
    int ix, iy, iz;
    double vx, vy, vz;
    bool crystal;
    double mass;
    
    Atom()
    : id(0), mol(0), type(0),
     x(0), y(0), z(0),
     ix(0), iy(0), iz(0),
     vx(0), vy(0), vz(0),
     crystal(false),
     mass(0)
    {}

};
struct Bond{
    long int id;
    int type;
    long int atom1, atom2;

    Bond()
    : id(0), type(0),
     atom1(0), atom2(0)
    {}
};
struct Angle{
    long int id;
    int type;
    long int atom1, atom2, atom3;

    Angle()
    : id(0), type(0),
     atom1(0), atom2(0), atom3(0)
    {}
};
struct Dihedral{
    long int id;
    int type;
    long int atom1, atom2, atom3, atom4;

    Dihedral()
    : id(0), type(0),
     atom1(0), atom2(0), atom3(0), atom4(0)
    {}
};
struct Improper{
    long int id;
    int type;
    long int atom1, atom2, atom3, atom4;

    Improper()
    : id(0), type(0),
     atom1(0), atom2(0), atom3(0), atom4(0)
    {}
};