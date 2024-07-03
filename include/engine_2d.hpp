#ifndef __ENGINE_2D_HPP__
#define __ENGINE_2D_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
//
#include "shape.hpp"
#include "engine.hpp"
#include "projection.hpp"

// Definitions


// Functions
complex_t get_phi_mn_2d(const complex_t k, 
    const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    quadl_domain_t quadl, int &flag);
void R_mn_2d(const real_t alpha_m, const real_t beta_m, const real_t alpha_n, const real_t beta_n, 
    const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    real_t &R_mn_mm, real_t &R_mn_mp, real_t &R_mn_pm, real_t &R_mn_pp);
void g_mn_2d(const real_t alpha_m, const real_t beta_m, const real_t alpha_n, const real_t beta_n, 
    const complex_t k, const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    complex_t &g_mn_mm, complex_t &g_mn_mp, complex_t &g_mn_pm, complex_t &g_mn_pp);
complex_t phi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    quadl_domain_t quadl, int &flag);


#endif