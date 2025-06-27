#include "../include/functions.h"

using namespace std;

// // Function to compute the inverse of a 3x3 matrix
// std::array<std::array<double, 3>, 3> inverse_3x3(const std::array<std::array<double, 3>, 3>& matrix) {
//     // Initialize output matrix
//     std::array<std::array<double, 3>, 3> inv_matrix{};

//     // Compute determinant
//     double det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
//                - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
//                + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

//     // Check for singular matrix
//     if (std::abs(det) < 1e-10) {
//         throw std::runtime_error("Matrix is singular (determinant is zero or near-zero)");
//     }

//     // Compute inverse: adjugate matrix divided by determinant
//     double inv_det = 1.0 / det;

//     // Cofactor matrix elements
//     inv_matrix[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * inv_det;
//     inv_matrix[0][1] = -(matrix[0][1] * matrix[2][2] - matrix[0][2] * matrix[2][1]) * inv_det;
//     inv_matrix[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * inv_det;
    
//     inv_matrix[1][0] = -(matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) * inv_det;
//     inv_matrix[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * inv_det;
//     inv_matrix[1][2] = -(matrix[0][0] * matrix[1][2] - matrix[0][2] * matrix[1][0]) * inv_det;
    
//     inv_matrix[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * inv_det;
//     inv_matrix[2][1] = -(matrix[0][0] * matrix[2][1] - matrix[0][1] * matrix[2][0]) * inv_det;
//     inv_matrix[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * inv_det;

//     return inv_matrix;
// }

void unwrap(System& system)
{
    if (system.box_type == 0) // orthogonal box
    {
        for( long int i = 0; i < system.num_atoms; i++)
        {
            system.atoms[i].x = system.atoms[i].x + system.atoms[i].ix * system.bx;
            system.atoms[i].y = system.atoms[i].y + system.atoms[i].iy * system.by;
            system.atoms[i].z = system.atoms[i].z + system.atoms[i].iz * system.bz;
        }
    }
    else if (system.box_type == 1) // triclinic box
    {
        for( long int i = 0; i < system.num_atoms; i++)
        {
            double H[3][3];
            H[0][0] = system.bx;
            H[0][1] = system.xtilt;
            H[0][2] = system.ytilt;
            H[1][0] = 0.0;
            H[1][1] = system.by;
            H[1][2] = system.ztilt;
            H[2][0] = 0.0;
            H[2][1] = 0.0;
            H[2][2] = system.bz;
            system.atoms[i].x = system.atoms[i].x + system.atoms[i].ix * H[0][0] + system.atoms[i].iy * H[0][1] + system.atoms[i].iz * H[0][2];
            system.atoms[i].y = system.atoms[i].y + system.atoms[i].ix * H[1][0] + system.atoms[i].iy * H[1][1] + system.atoms[i].iz * H[1][2];
            system.atoms[i].z = system.atoms[i].z + system.atoms[i].ix * H[2][0] + system.atoms[i].iy * H[2][1] + system.atoms[i].iz * H[2][2];
        }
    }
}
void wrap(System& system)
{
    if (system.box_type == 0) // orthogonal box
    {
        for( long int i = 0; i < system.num_atoms; i++)
        {
            system.atoms[i].x = fmod(system.atoms[i].x, system.bx);
            system.atoms[i].y = fmod(system.atoms[i].y, system.by);
            system.atoms[i].z = fmod(system.atoms[i].z, system.bz);
        }
    }
    else if (system.box_type == 1) // triclinic box
    {
        for( long int i = 0; i < system.num_atoms; i++)
        {
            double H[3][3];
            H[0][0] = system.bx;
            H[0][1] = system.xtilt;
            H[0][2] = system.ytilt;
            H[1][0] = 0.0;
            H[1][1] = system.by;
            H[1][2] = system.ztilt;
            H[2][0] = 0.0;
            H[2][1] = 0.0;
            H[2][2] = system.bz;
            // double invH[3][3];
            // invH = inverse_3x3(H);
            system.atoms[i].x = system.atoms[i].x - H[0][0] * system.atoms[i].ix;
            system.atoms[i].y = system.atoms[i].y - H[1][1] * system.atoms[i].iy;
            system.atoms[i].z = system.atoms[i].z - H[2][2] * system.atoms[i].iz;
        }
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