//
#include "engine.hpp"



engine_t::engine_t(){

}

engine_t::~engine_t(){

}

projection_1d_t get_projection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> p){
    projection_1d_t para_1d;
    vector_t<real_t> v21=v2-v1;
    real_t alpha=v21*(p-v1)/(v21*v21);
    vector_t<real_t> p0=v1+alpha*v21;
    vector_t<real_t> l=unit(v2-v1);
    para_1d.p0 = p0;
    para_1d.P0 = mag(p0-p);
    para_1d.P0_u = unit(p0-p);
    para_1d.l_m = (v1-p0)*l;
    para_1d.l_p = (v2-p0)*l;
    para_1d.P_m = mag(p-v1);
    para_1d.P_p = mag(p-v2);
    return para_1d;
}

void get_projection_2d_edge(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, projection_1d_t &para_1d,
    real_t &d, real_t &R0, real_t &R_m, real_t &R_p, vector_t<real_t> &u, vector_t<real_t> &n){
    vector_t<real_t> v21=v2-v1;
    vector_t<real_t> v31=v3-v1;
    vector_t<real_t> v32=v3-v2;
    real_t alpha=((v21*v31)*v31-(v31*v31)*v21)*(v1-p)/((v21*v21)*(v31*v31)-(v21*v31)*(v21*v31));
    real_t beta =((v21*v31)*v21-(v21*v21)*v31)*(v2-p)/((v21*v21)*(v31*v31)-(v21*v31)*(v21*v31));
    vector_t<real_t> p0=v1+alpha*v21+beta*v31;
    para_1d = get_projection_1d(v1, v2, p0);
    d = mag(p-p0);
    R0 = mag(p-para_1d.p0);
    R_m = mag(p-v1);
    R_p = mag(p-v2);
    n = unit(v21^v32);
    u = unit(v21^(v21^v32));
}

projection_2d_t get_projection_2d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p){
    projection_2d_t para_2d;
    get_projection_2d_edge(v2, v3, v1, p, para_2d.para_1d[0],
        para_2d.d, para_2d.R0[0], para_2d.R_m[0], 
        para_2d.R_p[0], para_2d.u[0], para_2d.n);
    get_projection_2d_edge(v3, v1, v2, p, para_2d.para_1d[1],
        para_2d.d, para_2d.R0[1], para_2d.R_m[1], 
        para_2d.R_p[1], para_2d.u[1], para_2d.n);
    get_projection_2d_edge(v1, v2, v3, p, para_2d.para_1d[2],
        para_2d.d, para_2d.R0[2], para_2d.R_m[2], 
        para_2d.R_p[2], para_2d.u[2], para_2d.n);
    return para_2d;
}

void get_projection_3d_triangle(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> p, vector_t<real_t> &n,
    projection_2d_t &para_2d){
    para_2d = get_projection_2d(v1, v2, v3, p);
    n = para_2d.n;
}

projection_3d_t get_projection_3d(const vector_t<real_t> v1, const vector_t<real_t> v2,
    const vector_t<real_t> v3, const vector_t<real_t> v4, const vector_t<real_t> p){
    projection_3d_t para_3d;
    get_projection_3d_triangle(v2, v3, v4, p, para_3d.n[0], para_3d.para_2d[0]);
    get_projection_3d_triangle(v3, v1, v4, p, para_3d.n[1], para_3d.para_2d[1]);
    get_projection_3d_triangle(v1, v2, v4, p, para_3d.n[2], para_3d.para_2d[2]);
    get_projection_3d_triangle(v3, v2, v1, p, para_3d.n[3], para_3d.para_2d[3]);
    return para_3d;
}

//

void R_mn_2d(const real_t alpha_m, const real_t beta_m, const real_t alpha_n, const real_t beta_n, 
    const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    real_t &R_mn_mm, real_t &R_mn_mp, real_t &R_mn_pm, real_t &R_mn_pp){
    vector_t<real_t> R_m_m = +1.0*(basis_m.r_m+alpha_m*basis_m.L_m1+beta_m*basis_m.L_m2);
    vector_t<real_t> R_m_p = -1.0*(basis_m.r_p+alpha_m*basis_m.L_p2+beta_m*basis_m.L_p1);
    vector_t<real_t> R_n_m = +1.0*(basis_n.r_m+alpha_n*basis_n.L_m1+beta_n*basis_n.L_m2);
    vector_t<real_t> R_n_p = -1.0*(basis_n.r_p+alpha_n*basis_n.L_p2+beta_n*basis_n.L_p1);
    R_mn_mm = mag(R_m_m-R_n_m);
    R_mn_mp = mag(R_m_m-R_n_p);
    R_mn_pm = mag(R_m_p-R_n_m);
    R_mn_pp = mag(R_m_p-R_n_p);
}

void g_mn_2d(const real_t alpha_m, const real_t beta_m, const real_t alpha_n, const real_t beta_n, 
    const complex_t k, const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    complex_t &g_mn_mm, complex_t &g_mn_mp, complex_t &g_mn_pm, complex_t &g_mn_pp){
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mn_mm, R_mn_mp, R_mn_pm, R_mn_pp;
    R_mn_2d(alpha_m, beta_m, alpha_n, beta_n, basis_m, basis_n, 
        R_mn_mm, R_mn_mp, R_mn_pm, R_mn_pp);
    g_mn_mm = exp(-j*k*R_mn_mm)/(4.0*pi*R_mn_mm);
    g_mn_mp = exp(-j*k*R_mn_mp)/(4.0*pi*R_mn_mp);
    g_mn_pm = exp(-j*k*R_mn_pm)/(4.0*pi*R_mn_pm);
    g_mn_pp = exp(-j*k*R_mn_pp)/(4.0*pi*R_mn_pp);
}

struct integrand_2d{
    basis_2d_t basis_m, basis_n;
    complex_t k;
};

const real_t eps_sinc=1.0E-15;

complex_t sinc(const complex_t x){
    return abs(x)<eps_sinc ? 1.0 : sin(x)/x;
}

real_t sinc(const real_t x){
    return abs(x)<eps_sinc ? 1.0 : sin(x)/x;
}

complex_t phi_mn_integrand_1_2d(real_t alpha_n, real_t beta_n, void *args_){
    integrand_2d *args=(integrand_2d*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    complex_t k=args->k;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mn_mm, R_mn_mp, R_mn_pm, R_mn_pp;
    R_mn_2d(1.0/3.0, 1.0/3.0, alpha_n, beta_n, basis_m, basis_n, 
        R_mn_mm, R_mn_mp, R_mn_pm, R_mn_pp);
    complex_t I1, I2, I3, I4;
    I1 = exp(-j*k*R_mn_mm/2.0)*sinc(k*R_mn_mm/2.0);
    I2 = exp(-j*k*R_mn_mp/2.0)*sinc(k*R_mn_mp/2.0);
    I3 = exp(-j*k*R_mn_pm/2.0)*sinc(k*R_mn_pm/2.0);
    I4 = exp(-j*k*R_mn_pp/2.0)*sinc(k*R_mn_pp/2.0);
    return (4.0*basis_m.L*basis_n.L)*(-j*k*(I1+I2+I3+I4)/(4.0*pi));
}

complex_t phi_mn_integrand_2_2d(void *args_){
    integrand_2d *args=(integrand_2d*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    complex_t k=args->k;
    vector_t<real_t> rho_m_m=
        vector_t<real_t>(+1.0*(basis_m.r_m+(1.0/3.0)*basis_m.e_1+(1.0/3.0)*basis_m.e_2));
    vector_t<real_t> rho_m_p=
        vector_t<real_t>(-1.0*(basis_m.r_p+(1.0/3.0)*basis_m.e_2+(1.0/3.0)*basis_m.e_1));
    projection_2d_t para_mm=get_projection_2d(basis_n.r_m, basis_n.e_1,
        basis_n.e_2, rho_m_m);
    projection_2d_t para_mp=get_projection_2d(basis_n.r_m, basis_n.e_1,
        basis_n.e_2, rho_m_p);
    projection_2d_t para_pm=get_projection_2d(basis_n.r_p, basis_n.e_2,
        basis_n.e_1, rho_m_m);
    projection_2d_t para_pp=get_projection_2d(basis_n.r_p, basis_n.e_2,
        basis_n.e_1, rho_m_p);
    complex_t I1=0.0, I2=0.0, I3=0.0, I4=0.0;
    for (size_t i=0; i<3; i++){
        I1+=(para_mm.para_1d[i].P0_u*para_mm.u[i])*(
            para_mm.para_1d[i].P0*log((para_mm.R_p[i]+para_mm.para_1d[i].l_p)/(para_mm.R_m[i]+para_mm.para_1d[i].l_m))-para_mm.d)*(
            atan2(para_mm.para_1d[i].P0*para_mm.para_1d[i].l_p, pow(para_mm.R0[i], 2.0)+para_mm.d*para_mm.R_p[i])-
            atan2(para_mm.para_1d[i].P0*para_mm.para_1d[i].l_m, pow(para_mm.R0[i], 2.0)+para_mm.d*para_mm.R_m[i]));
        I2+=(para_mp.para_1d[i].P0_u*para_mp.u[i])*(
            para_mp.para_1d[i].P0*log((para_mp.R_p[i]+para_mp.para_1d[i].l_p)/(para_mp.R_m[i]+para_mp.para_1d[i].l_m))-para_mp.d)*(
            atan2(para_mp.para_1d[i].P0*para_mp.para_1d[i].l_p, pow(para_mp.R0[i], 2.0)+para_mp.d*para_mp.R_p[i])-
            atan2(para_mp.para_1d[i].P0*para_mp.para_1d[i].l_m, pow(para_mp.R0[i], 2.0)+para_mp.d*para_mp.R_m[i]));
        I3+=(para_pm.para_1d[i].P0_u*para_pm.u[i])*(
            para_pm.para_1d[i].P0*log((para_pm.R_p[i]+para_pm.para_1d[i].l_p)/(para_pm.R_m[i]+para_pm.para_1d[i].l_m))-para_pm.d)*(
            atan2(para_pm.para_1d[i].P0*para_pm.para_1d[i].l_p, pow(para_pm.R0[i], 2.0)+para_pm.d*para_pm.R_p[i])-
            atan2(para_pm.para_1d[i].P0*para_pm.para_1d[i].l_m, pow(para_pm.R0[i], 2.0)+para_pm.d*para_pm.R_m[i]));
        I4+=(para_pp.para_1d[i].P0_u*para_pp.u[i])*(
            para_pp.para_1d[i].P0*log((para_pp.R_p[i]+para_pp.para_1d[i].l_p)/(para_pp.R_m[i]+para_pp.para_1d[i].l_m))-para_pp.d)*(
            atan2(para_pp.para_1d[i].P0*para_pp.para_1d[i].l_p, pow(para_pp.R0[i], 2.0)+para_pp.d*para_pp.R_p[i])-
            atan2(para_pp.para_1d[i].P0*para_pp.para_1d[i].l_m, pow(para_pp.R0[i], 2.0)+para_pp.d*para_pp.R_m[i]));
    }
    return (4.0*basis_m.L*basis_n.L)*(I1/basis_n.A_m+I2/basis_n.A_p+I3/basis_n.A_m+I4/basis_n.A_p)/(8.0*pi);
}

complex_t phi_mn_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k,
    quadl_domain_t &quadl){
    integrand_2d args={basis_m, basis_n, k};
    phi_mn_integrand_1_2d
    return quadl.integral_2d(, &args, triangle_domain_t()) 
}