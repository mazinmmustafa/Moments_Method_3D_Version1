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
struct RCS_2d{
    real_t sigma_theta=0.0;
    real_t sigma_phi=0.0;
};

class engine_2d_t{
    private:
    public:
        quadl_domain_t quadl;
        shape_t shape;
        matrix_t<complex_t> Z_mn, V_m, I_n;
        size_t N=0;
        engine_2d_t();
        ~engine_2d_t();
        void compute_Z_mn();
        void save_Z_mn(const char *filename);
        void load_Z_mn(const char *filename);
        void compute_V_m_plane_wave(const complex_t E_TM, const complex_t E_TE,
            const real_t theta_i, const real_t phi);
        RCS_2d RCS_plane_wave_2d(const real_t theta_s, const real_t phi_s);
};

// Functions
complex_t psi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    quadl_domain_t quadl, int &flag);
complex_t phi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    quadl_domain_t quadl, int &flag);
complex_t Z_mn_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const complex_t eta, quadl_domain_t quadl);


#endif