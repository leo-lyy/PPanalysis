// find the position vector of two atoms seperated by 6 bonds
// define the orientation vector ovec = (0, 1, 0)
// calculate the angle between vector and ovec
// calculate the order parameter by Second-Order Legendre Polynomial
// return the order parameter
# include "../include/p2calculation.h"
using namespace std;


double p2calculation(System& system)
{
    axis r1;
    long int p2num = 0;
    double p2sum = 0;
    // cout << system.num_atoms <<endl;
    // cout << system.atoms.size() <<endl;
    for (long int i = 0; i < system.num_atoms - 9; i++)
    {
        // cout << system.atoms[i].id << "  "<<system.atoms[i].mol << "  "<<system.atoms[i].type << "  "<<system.atoms[i].x << "  "<<system.atoms[i].y << "  "<<system.atoms[i].z << "  "<< endl;
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
    return p2sum / p2num;
}