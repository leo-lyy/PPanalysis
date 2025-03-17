#include "../include/functions.h"

using namespace std;

void unwrap(System& system)
{
    for( long int i = 0; i < system.num_atoms; i++)
    {
        system.atoms[i].x = system.atoms[i].x + system.atoms[i].ix * system.bx;
        system.atoms[i].y = system.atoms[i].y + system.atoms[i].iy * system.by;
        system.atoms[i].z = system.atoms[i].z + system.atoms[i].iz * system.bz;
    }
}
void wrap(System& system)
{
    for( long int i = 0; i < system.num_atoms; i++)
    {
        system.atoms[i].x = system.atoms[i].x - system.atoms[i].ix * system.bx;
        system.atoms[i].y = system.atoms[i].y - system.atoms[i].iy * system.by;
        system.atoms[i].z = system.atoms[i].z - system.atoms[i].iz * system.bz;
    }
}

double norm(const axis& rr)
{
    return sqrt(rr.x * rr.x + rr.y * rr.y + rr.z * rr.z);
}
void cross(axis& zz, axis& xx, axis yy)  // zz = xx cross yy
{
    zz.x = xx.y * yy.z - yy.y * xx.z;
    zz.y = yy.x * xx.z - yy.z * xx.x;
    zz.z = xx.x * yy.y - yy.x * xx.y;
}
double dot(axis& A, axis& B) // C = A · B
{
    return A.x * B.x + A.y * B.y + A.z * B.z;
}
axis unitVec(axis& vec)
{
    double n = norm(vec);
    return {
        vec.x / n,
        vec.y / n,
        vec.z / n
    };

}