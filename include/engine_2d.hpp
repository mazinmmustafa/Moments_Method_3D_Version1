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
#include "projection.hpp"

// Definitions
struct RCS_2d_t{
    real_t sigma_theta=0.0;
    real_t sigma_phi=0.0;
};

struct field_2d_t{
    vector_t<complex_t> E=vector_t<complex_t>(0.0, 0.0, 0.0);
    vector_t<complex_t> H=vector_t<complex_t>(0.0, 0.0, 0.0);
};

class engine_2d_t{
    private:
        quadl_domain_t quadl;
        const size_t k_max=40;
        const real_t tol=1.0E-3;
        int is_Z_mn_available=false;
        int is_V_m_available=false;
        int is_I_n_available=false;
        int is_mesh_obtained=false;
    public:
        shape_t shape;
        size_t N=0;
        complex_t k=0.0, eta=0.0;
        real_t freq=0.0, lambda=0.0;
        matrix_t<complex_t> Z_mn, V_m, I_n;
        engine_2d_t();
        ~engine_2d_t();
        void compute_Z_mn();
        void save_Z_mn(const char *filename);
        void load_Z_mn(const char *filename);
        void compute_V_m_plane_wave(const complex_t E_TM, const complex_t E_TE,
            const real_t theta_i, const real_t phi);
        RCS_2d_t RCS_plane_wave_2d(const real_t theta_s, const real_t phi_s);
        field_2d_t compute_near_field(const vector_t<real_t> p);
        //
        size_t get_N(){return this->N;}
        void set_medium(const complex_t mu, const complex_t eps, const real_t freq);
        void mesh(const char *filename, const real_t clmax);
        void solve_currents();
        void unset();
        shape_info_t get_shape_info(){
            return this->shape.get_shape_info();
        }
        field_2d_t compute_incident_plane_wave_field(const real_t theta_i, const real_t phi_i,
            const complex_t E_TM, const complex_t E_TE, const vector_t<real_t> p);
        void export_currents();
};

// Functions
complex_t psi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const real_t lambda, quadl_domain_t quadl, int &flag);
complex_t phi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const real_t lambda, quadl_domain_t quadl, int &flag);
complex_t Z_mn_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const real_t lambda, const complex_t eta, quadl_domain_t quadl);

#endif