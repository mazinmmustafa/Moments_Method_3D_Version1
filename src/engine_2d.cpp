//
#include "engine_2d.hpp"

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

struct integrand_2d_args{
    basis_2d_t basis_m, basis_n;
    complex_t k;
    real_t alpha, beta;
};

complex_t integrand_2d_phi(const complex_t alpha_, const complex_t beta_, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    complex_t k=args->k;
    const complex_t j=complex_t(0.0, 1.0);
    real_t R_mm, R_mp, R_pm, R_pp;
    real_t alpha, beta;
    alpha = beta = 1.0/3.0;
    R_mn_2d(alpha, beta, real(alpha_), real(beta_), 
        basis_m, basis_n, R_mm, R_mp, R_pm, R_pp);
    real_t Lm=basis_m.L;
    real_t Ln=basis_n.L;
    real_t factor=(Lm*Ln)/pi;
    complex_t I_mm=-j*k*exp(-j*k*R_mm/2.0)*sinc(k*R_mm/2.0);
    complex_t I_mp=-j*k*exp(-j*k*R_mp/2.0)*sinc(k*R_mp/2.0);
    complex_t I_pm=-j*k*exp(-j*k*R_pm/2.0)*sinc(k*R_pm/2.0);
    complex_t I_pp=-j*k*exp(-j*k*R_pp/2.0)*sinc(k*R_pp/2.0);
    return factor*(I_mm-I_mp-I_pm+I_pp);
}

complex_t term_2d_phi(const basis_2d_t &basis_m, const basis_2d_t &basis_n){
    real_t alpha, beta;
    alpha = beta = 1.0/3.0;
    projection_2d_t para_m_m, para_m_p, para_p_m, para_p_p;
    vector_t<real_t> p_m_m=+1.0*(basis_m.r_m+alpha*basis_m.L_m1+beta*basis_m.L_m2);
    vector_t<real_t> p_m_p=+1.0*(basis_m.r_p-alpha*basis_m.L_p1-beta*basis_m.L_p2);
    real_t I_m_m, I_m_p, I_p_m, I_p_p;
    para_m_m = get_projection_2d(basis_n.r_m, basis_n.e_1, basis_m.e_2, p_m_m);    
    para_m_p = get_projection_2d(basis_n.r_m, basis_n.e_1, basis_m.e_2, p_m_p);    
    para_p_m = get_projection_2d(basis_n.r_p, basis_n.e_2, basis_m.e_1, p_m_m);  
    para_p_p = get_projection_2d(basis_n.r_p, basis_n.e_2, basis_m.e_1, p_m_p);  
    real_t R_p, R_m, R0;
    real_t l_p, l_m, P0;
    real_t d;
    real_t term_1, term_2, term_3, term_4;
    vector_t<real_t> u, P0_u;
    for (size_t i=0; i<3; i++){
        // term mm
        R0 = para_m_m.R0[i];
        R_m = para_m_m.R_m[i];
        R_p = para_m_m.R_p[i];
        d = para_m_m.d;
        P0 = para_m_m.para_1d[i].P0;
        l_m = para_m_m.para_1d[i].l_m;
        l_p = para_m_m.para_1d[i].l_p;
        u = para_m_m.u[i];
        P0_u = para_m_m.para_1d[i].P0_u;
        term_1 = P0_u*u;
        term_2 = P0*log((R_p+l_p)/(R_m+l_m));
        term_3 = atan2(P0*l_p, R0*R0+d*R_p);
        term_4 = atan2(P0*l_m, R0*R0+d*R_m);
        I_m_m+=term_1*(term_2-d*(term_3-term_4));
        // term mp
        R0 = para_m_p.R0[i];
        R_m = para_m_p.R_m[i];
        R_p = para_m_p.R_p[i];
        d = para_m_p.d;
        P0 = para_m_p.para_1d[i].P0;
        l_m = para_m_p.para_1d[i].l_m;
        l_p = para_m_p.para_1d[i].l_p;
        u = para_m_p.u[i];
        P0_u = para_m_p.para_1d[i].P0_u;
        term_1 = P0_u*u;
        term_2 = P0*log((R_p+l_p)/(R_m+l_m));
        term_3 = atan2(P0*l_p, R0*R0+d*R_p);
        term_4 = atan2(P0*l_m, R0*R0+d*R_m);
        I_m_p+=term_1*(term_2-d*(term_3-term_4));
        // term pm
        R0 = para_p_m.R0[i];
        R_m = para_p_m.R_m[i];
        R_p = para_p_m.R_p[i];
        d = para_p_m.d;
        P0 = para_p_m.para_1d[i].P0;
        l_m = para_p_m.para_1d[i].l_m;
        l_p = para_p_m.para_1d[i].l_p;
        u = para_p_m.u[i];
        P0_u = para_p_m.para_1d[i].P0_u;
        term_1 = P0_u*u;
        term_2 = P0*log((R_p+l_p)/(R_m+l_m));
        term_3 = atan2(P0*l_p, R0*R0+d*R_p);
        term_4 = atan2(P0*l_m, R0*R0+d*R_m);
        I_p_m+=term_1*(term_2-d*(term_3-term_4));
        // term pp
        R0 = para_p_p.R0[i];
        R_m = para_p_p.R_m[i];
        R_p = para_p_p.R_p[i];
        d = para_p_p.d;
        P0 = para_p_p.para_1d[i].P0;
        l_m = para_p_p.para_1d[i].l_m;
        l_p = para_p_p.para_1d[i].l_p;
        u = para_p_p.u[i];
        P0_u = para_p_p.para_1d[i].P0_u;
        term_1 = P0_u*u;
        term_2 = P0*log((R_p+l_p)/(R_m+l_m));
        term_3 = atan2(P0*l_p, R0*R0+d*R_p);
        term_4 = atan2(P0*l_m, R0*R0+d*R_m);
        I_p_p+=term_1*(term_2-d*(term_3-term_4));
    }
    real_t Lm=basis_m.L;
    real_t Ln=basis_n.L;
    real_t factor=(Lm*Ln)/pi;
    real_t A_n_m=basis_n.A_m;
    real_t A_n_p=basis_n.A_p;
    return factor*(I_m_m/(2.0*A_n_m)-I_m_p/(2.0*A_n_p)-I_p_m/(2.0*A_n_m)+I_p_p/(2.0*A_n_p));
}

complex_t get_phi_mn_2d(const complex_t k, 
    const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    quadl_domain_t quadl, int &flag){
    integrand_2d_args args={basis_m, basis_n, k, 0.0, 0.0};
    triangle_domain_t triangle;
    vector_t<real_t> v1(0.0, 0.0, 0.0);
    vector_t<real_t> v2(1.0, 0.0, 0.0);
    vector_t<real_t> v3(0.0, 1.0, 0.0);
    triangle.v1 = v1;
    triangle.v2 = v2;
    triangle.v3 = v3;
    complex_t I1=quadl.integral_2d(integrand_2d_phi, &args, triangle, flag);
    complex_t I2=term_2d_phi(basis_m, basis_n);
    print(I1);
    print(I2);
    return I1+I2;
}