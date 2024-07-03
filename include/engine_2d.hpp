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
class engine_2d_t{
    private:
    public:
        quadl_domain_t quadl;
        shape_t shape;
        matrix_t<complex_t> Z_mn;
        engine_2d_t();
        ~engine_2d_t();
        void compute_Z_mn();
        void save_Z_mn(const char *filename);
        
};

// Functions
complex_t psi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    quadl_domain_t quadl, int &flag);
complex_t phi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    quadl_domain_t quadl, int &flag);
complex_t Z_mn_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const complex_t eta, quadl_domain_t quadl);

#endif