#ifndef __ENGINE_HPP__
#define __ENGINE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
//
#include "shape.hpp"

// Definitions
class engine_t{
    private:
    public:
        engine_t();
        ~engine_t();
};

struct projection_1d_t{
    vector_t<real_t> p0, P0_u;
    real_t l_m, l_p;
    real_t P0, P_m, P_p;
};

struct projection_2d_t{
    projection_1d_t para_1d[3];
    vector_t<real_t> u[3];
    real_t R0[3], R_m[3], R_p[3], d;
    vector_t<real_t> n;
};

struct projection_3d_t{
    projection_2d_t para_2d[4];
    vector_t<real_t> n[4];
};

// Functions

#endif