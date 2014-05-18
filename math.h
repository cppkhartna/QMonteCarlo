#pragma once
//#include "const.h"
#include <cmath>

struct vec_3d
{
    double v[3];
public:
    vec_3d();
    vec_3d(double v0, double v1, double v2);
    vec_3d(const vec_3d& v1);
    ~vec_3d(){};
    vec_3d& operator=(const vec_3d &v1);
    double& operator[](int idx);
    double length2();
    inline void add_random(int k, double& change, double& x_max, double& x_min)
    { 
        double v_k = v[k];
        v_k += change;
        if (v_k > x_max)
            v_k = x_max;
        if (v_k < x_min)
            v_k = x_min;
        v[k] = v_k; 
    };
    friend vec_3d operator+(const vec_3d &v1, const vec_3d &v2);
    friend vec_3d operator-(const vec_3d &v1, const vec_3d &v2);
    inline friend double sub_len(const vec_3d& v1, const vec_3d& v2)
    {
        double a0 = v1.v[0] - v2.v[0];
        double a1 = v1.v[1] - v2.v[1];
        double a2 = v1.v[2] - v2.v[2];
        return sqrt(a0*a0 + a1*a1 + a2*a2);
    };
};
