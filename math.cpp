#include "math.h"
using namespace std;

vec_3d::vec_3d()
{
    for (int i = 0; i < 3; i++)
        v[0] = 0;
};

vec_3d::vec_3d(double v0, double v1, double v2)
{
    v[0] = v0; v[1] = v1; v[2] = v2;
};

vec_3d::vec_3d(const vec_3d &v1)
{
    for (int i = 0; i < 3; i++)
        v[i] = v1.v[i];
};

vec_3d& vec_3d::operator=(const vec_3d &v1)
{
    for (int i = 0; i < 3; i++)
        v[i] = v1.v[i];
    return *this;
};

double& vec_3d::operator[](int idx)
{
    return v[idx];
};

vec_3d operator+(const vec_3d &v1, const vec_3d &v2)
{
    return vec_3d(v1.v[0] + v2.v[0], v1.v[1] + v2.v[1], v1.v[2] + v2.v[2]);
};

vec_3d operator-(const vec_3d &v1, const vec_3d &v2)
{
    return vec_3d(v1.v[0] - v2.v[0], v1.v[1] - v2.v[1], v1.v[2] - v2.v[2]);
};

double vec_3d::length2()
{
    return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
};
