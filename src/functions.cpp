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
void vecInit(vector<vector<axis>>& v, long int molNum, const long int maxhelix)
{
    v.resize(molNum + 1);
    for (long int i = 1; i<= molNum; i++)
    {
        v[i].resize(maxhelix);
    }
}
double degVec(axis& a, axis& b)
{

    double norm_a = norm(a);
    double norm_b = norm(b);

    if (norm_a == 0 || norm_b == 0) {

        throw std::invalid_argument("One of the vectors has zero length");
    }

    double dot_product = a.x * b.x + a.y * b.y + a.z * b.z;
    double cos_theta = dot_product / (norm_a * norm_b);

    // Ensure the value is within the valid range for acos
    if (cos_theta < -1.0) cos_theta = -1.0;
    if (cos_theta > 1.0) cos_theta = 1.0;

    double angle_radians = acos(cos_theta);
    double angle_degrees = angle_radians * (180.0 / M_PI);

    return angle_degrees;

}
void printProgressBar(int k, int kmax) {
    int barWidth = 50;  // Width of the progress bar
    float progress = static_cast<float>(k) / kmax;
    int pos = static_cast<int>(barWidth * progress);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " % (" << k << "/" << kmax << ")\r";
    std::cout.flush();
}