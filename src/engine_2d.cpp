//
#include "engine_2d.hpp"

const size_t line_max=100;

struct integrand_2d_args{
    basis_2d_t basis_m, basis_n;
    complex_t k=0.0;
    quadl_domain_t quadl;
    complex_t alpha_m=0.0, beta_m=0.0;
    real_t theta_i=0.0, phi_i=0.0;
    real_t theta_s=0.0, phi_s=0.0;
    complex_t E_TM=0.0, E_TE=0.0;
    vector_t<real_t> p=vector_t<real_t>(0.0, 0.0, 0.0); 
    vector_t<real_t> direction=vector_t<real_t>(0.0, 0.0, 0.0);
    real_t lambda=0.0;
    engine_2d_t *engine=null;
};

void R_mn_2d(const real_t alpha_m, const real_t beta_m, const real_t alpha_n, const real_t beta_n, 
    const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    real_t &R_mn_mm, real_t &R_mn_mp, real_t &R_mn_pm, real_t &R_mn_pp){
    vector_t<real_t> R_m_m = +1.0*(basis_m.r_m+alpha_m*basis_m.L_m1+beta_m*basis_m.L_m2);
    vector_t<real_t> R_m_p = +1.0*(basis_m.r_p+alpha_m*basis_m.L_p1+beta_m*basis_m.L_p2);
    vector_t<real_t> R_n_m = +1.0*(basis_n.r_m+alpha_n*basis_n.L_m1+beta_n*basis_n.L_m2);
    vector_t<real_t> R_n_p = +1.0*(basis_n.r_p+alpha_n*basis_n.L_p1+beta_n*basis_n.L_p2);
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

// phi terms

complex_t integrand_phi_2d_inner(const complex_t alpha_n, const complex_t beta_n, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    complex_t k=args->k;
    real_t alpha_m=real(args->alpha_m);
    real_t beta_m=real(args->beta_m);
    real_t R_mm, R_mp, R_pm, R_pp;
    R_mn_2d(alpha_m, beta_m, real(alpha_n), real(beta_n), basis_m, basis_n, R_mm, R_mp, R_pm, R_pp);
    complex_t I_mm, I_mp, I_pm, I_pp;
    const complex_t j=complex_t(0.0, 1.0);
    I_mm = -j*k*exp(-j*k*R_mm/2.0)*sinc(k*R_mm/2.0);
    I_mp = -j*k*exp(-j*k*R_mp/2.0)*sinc(k*R_mp/2.0);
    I_pm = -j*k*exp(-j*k*R_pm/2.0)*sinc(k*R_pm/2.0);
    I_pp = -j*k*exp(-j*k*R_pp/2.0)*sinc(k*R_pp/2.0);
    return (basis_m.L*basis_n.L/pi)*(I_mm-I_mp-I_pm+I_pp);
}

complex_t integrand_phi_2d_outer(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_phi_2d_inner, args, triangle, flag);
}

complex_t integrand_phi_2d_projection(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    real_t lambda=args->lambda;
    real_t l_m, l_p, R_m, R_p, R0, P0, d;
    vector_t<real_t> P0_u, u;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m=basis_m.r_m+real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2;
    vector_t<real_t> rho_p=basis_m.r_p+real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2;
    projection_2d_t para_mm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_m, lambda);
    projection_2d_t para_mp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_m, lambda);
    projection_2d_t para_pm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_p, lambda);
    projection_2d_t para_pp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_p, lambda);
    //
    I_mm = I_mp = I_pm = I_pp = 0.0;
    real_t A, B, C, D;
    for (size_t i=0; i<3; i++){
        // mm
        R0 = para_mm.R0[i];
        R_m = para_mm.R_m[i];
        R_p = para_mm.R_p[i];
        l_m = para_mm.para_1d[i].l_m;
        l_p = para_mm.para_1d[i].l_p;
        P0 = para_mm.para_1d[i].P0;
        P0_u = para_mm.para_1d[i].P0_u;
        u = para_mm.u[i];
        d = para_mm.d;
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_mm+=D*(A-d*(B-C));
        // mp
        R0 = para_mp.R0[i];
        R_m = para_mp.R_m[i];
        R_p = para_mp.R_p[i];
        l_m = para_mp.para_1d[i].l_m;
        l_p = para_mp.para_1d[i].l_p;
        P0 = para_mp.para_1d[i].P0;
        P0_u = para_mp.para_1d[i].P0_u;
        u = para_mp.u[i];
        d = para_mp.d;
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_mp+=D*(A-d*(B-C));
        // pm
        R0 = para_pm.R0[i];
        R_m = para_pm.R_m[i];
        R_p = para_pm.R_p[i];
        l_m = para_pm.para_1d[i].l_m;
        l_p = para_pm.para_1d[i].l_p;
        P0 = para_pm.para_1d[i].P0;
        P0_u = para_pm.para_1d[i].P0_u;
        u = para_pm.u[i];
        d = para_pm.d;
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_pm+=D*(A-d*(B-C));        
        // pp
        R0 = para_pp.R0[i];
        R_m = para_pp.R_m[i];
        R_p = para_pp.R_p[i];
        l_m = para_pp.para_1d[i].l_m;
        l_p = para_pp.para_1d[i].l_p;
        P0 = para_pp.para_1d[i].P0;
        P0_u = para_pp.para_1d[i].P0_u;
        u = para_pp.u[i];
        d = para_pp.d;
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_pp+=D*(A-d*(B-C));
    }   
    return (basis_m.L*basis_n.L/(2.0*pi))*(I_mm/basis_n.A_m-I_mp/basis_n.A_p-I_pm/basis_n.A_m+I_pp/basis_n.A_p);
}

complex_t phi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const real_t lambda, quadl_domain_t quadl, int &flag){
    integrand_2d_args args={basis_m, basis_n, k, quadl, 0, 0};
    args.lambda = lambda;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    return args.quadl.integral_2d(integrand_phi_2d_projection, &args, triangle, flag)+
        args.quadl.integral_2d(integrand_phi_2d_outer, &args, triangle, flag);
}

// psi terms

complex_t integrand_psi_2d_inner(const complex_t alpha_n, const complex_t beta_n, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    complex_t k=args->k;
    real_t alpha_m=real(args->alpha_m);
    real_t beta_m=real(args->beta_m);
    real_t R_mm, R_mp, R_pm, R_pp;
    R_mn_2d(alpha_m, beta_m, real(alpha_n), real(beta_n), basis_m, basis_n, R_mm, R_mp, R_pm, R_pp);
    complex_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m_m, rho_m_p, rho_n_m, rho_n_p;
    rho_m_m = +1.0*(alpha_m*basis_m.L_m1+beta_m*basis_m.L_m2);
    rho_m_p = -1.0*(alpha_m*basis_m.L_p1+beta_m*basis_m.L_p2);
    rho_n_m = +1.0*(real(alpha_n)*basis_n.L_m1+real(beta_n)*basis_n.L_m2);
    rho_n_p = -1.0*(real(alpha_n)*basis_n.L_p1+real(beta_n)*basis_n.L_p2);
    const complex_t j=complex_t(0.0, 1.0);
    I_mm = -j*k*(rho_m_m*rho_n_m)*exp(-j*k*R_mm/2.0)*sinc(k*R_mm/2.0);
    I_mp = -j*k*(rho_m_m*rho_n_p)*exp(-j*k*R_mp/2.0)*sinc(k*R_mp/2.0);
    I_pm = -j*k*(rho_m_p*rho_n_m)*exp(-j*k*R_pm/2.0)*sinc(k*R_pm/2.0);
    I_pp = -j*k*(rho_m_p*rho_n_p)*exp(-j*k*R_pp/2.0)*sinc(k*R_pp/2.0);
    return (basis_m.L*basis_n.L/(4.0*pi))*(I_mm+I_mp+I_pm+I_pp);
}

complex_t integrand_psi_2d_outer(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_psi_2d_inner, args, triangle, flag);
}

complex_t integrand_psi_2d_projection_1(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    real_t lambda=args->lambda;
    real_t l_m, l_p, R_m, R_p, R0, P0, d;
    vector_t<real_t> P0_u, u;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m=basis_m.r_m+real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2;
    vector_t<real_t> rho_p=basis_m.r_p+real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2;
    projection_2d_t para_mm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_m, lambda);
    projection_2d_t para_mp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_m, lambda);
    projection_2d_t para_pm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_p, lambda);
    projection_2d_t para_pp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_p, lambda);
    vector_t<real_t> rho_n_m, rho_n_p;
    vector_t<real_t> rho_m_m=+1.0*(real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2);
    vector_t<real_t> rho_m_p=-1.0*(real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2);
    //
    I_mm = I_mp = I_pm = I_pp = 0.0;
    real_t A, B, C, D;
    for (size_t i=0; i<3; i++){
        // mm
        R0 = para_mm.R0[i];
        R_m = para_mm.R_m[i];
        R_p = para_mm.R_p[i];
        l_m = para_mm.para_1d[i].l_m;
        l_p = para_mm.para_1d[i].l_p;
        P0 = para_mm.para_1d[i].P0;
        P0_u = para_mm.para_1d[i].P0_u;
        u = para_mm.u[i];
        d = para_mm.d;
        rho_n_m = basis_n.r_m-para_mm.p0[i];
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_mm+=(rho_m_m*rho_n_m)*D*(A-d*(B-C));
        // mp
        R0 = para_mp.R0[i];
        R_m = para_mp.R_m[i];
        R_p = para_mp.R_p[i];
        l_m = para_mp.para_1d[i].l_m;
        l_p = para_mp.para_1d[i].l_p;
        P0 = para_mp.para_1d[i].P0;
        P0_u = para_mp.para_1d[i].P0_u;
        u = para_mp.u[i];
        d = para_mp.d;
        rho_n_p = basis_n.r_p-para_mp.p0[i];
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_mp+=(rho_m_m*rho_n_p)*D*(A-d*(B-C));
        // pm
        R0 = para_pm.R0[i];
        R_m = para_pm.R_m[i];
        R_p = para_pm.R_p[i];
        l_m = para_pm.para_1d[i].l_m;
        l_p = para_pm.para_1d[i].l_p;
        P0 = para_pm.para_1d[i].P0;
        P0_u = para_pm.para_1d[i].P0_u;
        u = para_pm.u[i];
        d = para_pm.d;
        rho_n_m = basis_n.r_m-para_pm.p0[i];
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_pm+=(rho_m_p*rho_n_m)*D*(A-d*(B-C));
        // pp
        R0 = para_pp.R0[i];
        R_m = para_pp.R_m[i];
        R_p = para_pp.R_p[i];
        l_m = para_pp.para_1d[i].l_m;
        l_p = para_pp.para_1d[i].l_p;
        P0 = para_pp.para_1d[i].P0;
        P0_u = para_pp.para_1d[i].P0_u;
        u = para_pp.u[i];
        d = para_pp.d;
        rho_n_p = basis_n.r_p-para_pp.p0[i];
        A = P0*log((R_p+l_p)/(R_m+l_m));
        B = atan2(P0*l_p, R0*R0+d*R_p);
        C = atan2(P0*l_m, R0*R0+d*R_m);
        D = P0_u*u;
        I_pp+=(rho_m_p*rho_n_p)*D*(A-d*(B-C));
    }   
    return (basis_m.L*basis_n.L/(8.0*pi))*(-I_mm/basis_n.A_m+I_mp/basis_n.A_p-I_pm/basis_n.A_m+I_pp/basis_n.A_p);
}

complex_t integrand_psi_2d_projection_2(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_m=args->basis_m;
    basis_2d_t basis_n=args->basis_n;
    real_t lambda=args->lambda;
    real_t l_m, l_p, R_m, R_p, R0;
    vector_t<real_t> u;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m=basis_m.r_m+real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2;
    vector_t<real_t> rho_p=basis_m.r_p+real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2;
    projection_2d_t para_mm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_m, lambda);
    projection_2d_t para_mp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_m, lambda);
    projection_2d_t para_pm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_p, lambda);
    projection_2d_t para_pp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_p, lambda);
    vector_t<real_t> rho_m_m=+1.0*(real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2);
    vector_t<real_t> rho_m_p=-1.0*(real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2);
    //
    I_mm = I_mp = I_pm = I_pp = 0.0;
    real_t A, B, C;
    for (size_t i=0; i<3; i++){
        // mm
        R0 = para_mm.R0[i];
        R_m = para_mm.R_m[i];
        R_p = para_mm.R_p[i];
        l_m = para_mm.para_1d[i].l_m;
        l_p = para_mm.para_1d[i].l_p;
        u = para_mm.u[i];
        A = R0*R0*log((R_p+l_p)/(R_m+l_m));
        B = R_p*l_p;
        C = R_m*l_m;
        I_mm+=(0.5*rho_m_m*u)*(A+B-C);
        // mp
        R0 = para_mp.R0[i];
        R_m = para_mp.R_m[i];
        R_p = para_mp.R_p[i];
        l_m = para_mp.para_1d[i].l_m;
        l_p = para_mp.para_1d[i].l_p;
        u = para_mp.u[i];
        A = R0*R0*log((R_p+l_p)/(R_m+l_m));
        B = R_p*l_p;
        C = R_m*l_m;
        I_mp+=(0.5*rho_m_m*u)*(A+B-C);
        // pm
        R0 = para_pm.R0[i];
        R_m = para_pm.R_m[i];
        R_p = para_pm.R_p[i];
        l_m = para_pm.para_1d[i].l_m;
        l_p = para_pm.para_1d[i].l_p;
        u = para_pm.u[i];
        A = R0*R0*log((R_p+l_p)/(R_m+l_m));
        B = R_p*l_p;
        C = R_m*l_m;
        I_pm+=(0.5*rho_m_p*u)*(A+B-C);
        // pp
        R0 = para_pp.R0[i];
        R_m = para_pp.R_m[i];
        R_p = para_pp.R_p[i];
        l_m = para_pp.para_1d[i].l_m;
        l_p = para_pp.para_1d[i].l_p;
        u = para_pp.u[i];
        A = R0*R0*log((R_p+l_p)/(R_m+l_m));
        B = R_p*l_p;
        C = R_m*l_m;
        I_pp+=(0.5*rho_m_p*u)*(A+B-C);
    }   
    return (basis_m.L*basis_n.L/(8.0*pi))*(I_mm/basis_n.A_m-I_mp/basis_n.A_p+I_pm/basis_n.A_m-I_pp/basis_n.A_p);
}

complex_t psi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k,
    const real_t lambda, quadl_domain_t quadl, int &flag){
    integrand_2d_args args={basis_m, basis_n, k, quadl, 0, 0};
    args.lambda = lambda;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    return args.quadl.integral_2d(integrand_psi_2d_projection_1, &args, triangle, flag)+
            args.quadl.integral_2d(integrand_psi_2d_projection_2, &args, triangle, flag)+
            args.quadl.integral_2d(integrand_psi_2d_outer, &args, triangle, flag);
}

complex_t Z_mn_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const real_t lambda, const complex_t eta, quadl_domain_t quadl){
    int flag=false;
    complex_t psi=psi_2d(basis_m, basis_n, k, lambda, quadl, flag); if(flag){flag=false; print("warning: no convergence!\n");}
    complex_t phi=phi_2d(basis_m, basis_n, k, lambda, quadl, flag); if(flag){flag=false; print("warning: no convergence!\n");}
    const complex_t j=complex_t(0.0, 1.0);
    complex_t Z=j*k*eta*psi-j*(eta/k)*phi;
    assert_error(!isnan(abs(Z)), "nan value for Z_mn");
    assert_error(!isinf(abs(Z)), "inf value for Z_mn");
    return Z;
}

complex_t integrand_V_plane_wave_2d(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_m=args->basis_m;
    real_t theta_i=args->theta_i;
    real_t phi_i=args->phi_i;
    complex_t k=args->k;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> rho_m=+1.0*(real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2);
    vector_t<real_t> rho_p=-1.0*(real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2);
    vector_t<real_t> k_i(sin(theta_i)*cos(phi_i), 
                         sin(theta_i)*sin(phi_i), 
                         cos(theta_i));
    complex_t E_TM=args->E_TM;
    complex_t E_TE=args->E_TE;
    complex_t k_i_r_m=k*(k_i*(basis_m.r_m+rho_m));
    complex_t k_i_r_p=k*(k_i*(basis_m.r_p-rho_p));
    vector_t<real_t> theta_i_u(cos(theta_i)*cos(phi_i),
                               cos(theta_i)*sin(phi_i),
                               -sin(theta_i));
    vector_t<real_t> phi_i_u(-sin(phi_i), cos(phi_i), 0.0);
    //
    complex_t xi_m, xi_p;
    xi_m = basis_m.L*exp(j*k_i_r_m)*(rho_m*theta_i_u);
    xi_p = basis_m.L*exp(j*k_i_r_p)*(rho_p*theta_i_u);
    complex_t I_TM=E_TM*(xi_m+xi_p);
    xi_m = basis_m.L*exp(j*k_i_r_m)*(rho_m*phi_i_u);
    xi_p = basis_m.L*exp(j*k_i_r_p)*(rho_p*phi_i_u);
    complex_t I_TE=E_TE*(xi_m+xi_p);
    return I_TM+I_TE;
}

complex_t V_m_plane_wave_2d(const basis_2d_t basis_m, const complex_t E_TM, const complex_t E_TE, 
    const complex_t k, const real_t theta_i, const real_t phi_i, quadl_domain_t quadl){
    int flag=false;
    integrand_2d_args args;
    args.basis_m = basis_m;
    args.k = k;
    args.E_TM = E_TM;
    args.E_TE = E_TE;
    args.theta_i = theta_i;
    args.phi_i = phi_i;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    complex_t I_V=quadl.integral_2d(integrand_V_plane_wave_2d, &args, triangle, flag);
    if(flag){flag=false; print("warning: no convergence!\n");}
    return I_V;
}

// near-field terms

complex_t integrand_E_2d(const complex_t alpha_n, const complex_t beta_n, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    vector_t<real_t> p=args->p;
    vector_t<real_t> direction=args->direction;
    vector_t<real_t> rho_m, rho_p;
    basis_2d_t basis_n;
    engine_2d_t *engine=args->engine;
    complex_t k=engine->k;
    complex_t eta=engine->eta;
    complex_t sum=0.0;
    for (size_t n=0; n<engine->N; n++){
        basis_n = engine->shape.get_basis_2d(n);
        rho_m = +1.0*(real(alpha_n)*basis_n.L_m1+real(beta_n)*basis_n.L_m2);
        rho_p = -1.0*(real(alpha_n)*basis_n.L_p1+real(beta_n)*basis_n.L_p2);
        real_t R_m, R_p;
        R_m = mag(p-(basis_n.r_m+rho_m));
        R_p = mag(p-(basis_n.r_p-rho_p));
        const complex_t j=complex_t(0.0, 1.0);
        complex_t g_m=exp(-j*k*R_m)/(4.0*pi*R_m);
        complex_t g_p=exp(-j*k*R_p)/(4.0*pi*R_p);
        vector_t<real_t> R_m_u=unit(p-(basis_n.r_m+rho_m));
        vector_t<real_t> R_p_u=unit(p-(basis_n.r_p-rho_p));
        complex_t I_A_m, I_A_p, I_B_m, I_B_p;
        I_A_m = (direction*rho_m)*g_m;
        I_A_p = (direction*rho_p)*g_p;
        I_B_m = -2.0*((1.0+j*k*R_m)/R_m)*(direction*R_m_u)*g_m;
        I_B_p = -2.0*((1.0+j*k*R_p)/R_p)*(direction*R_p_u)*g_p;
        complex_t A, B;
        A = (I_A_m+I_A_p);
        B = (I_B_m-I_B_p);
        sum+=-j*k*eta*basis_n.L*(A+B/(k*k))*engine->I_n(n, 0);
    }
    return sum;
}

complex_t integrand_H_2d(const complex_t alpha_n, const complex_t beta_n, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    vector_t<real_t> p=args->p;
    vector_t<real_t> direction=args->direction;
    vector_t<real_t> rho_m, rho_p;
    basis_2d_t basis_n;
    engine_2d_t *engine=args->engine;
    complex_t k=engine->k;
    complex_t sum=0.0;
    for (size_t n=0; n<engine->N; n++){
        basis_n = engine->shape.get_basis_2d(n);
        rho_m = +1.0*(real(alpha_n)*basis_n.L_m1+real(beta_n)*basis_n.L_m2);
        rho_p = -1.0*(real(alpha_n)*basis_n.L_p1+real(beta_n)*basis_n.L_p2);
        real_t R_m, R_p;
        R_m = mag(p-(basis_n.r_m+rho_m));
        R_p = mag(p-(basis_n.r_p-rho_p));
        const complex_t j=complex_t(0.0, 1.0);
        complex_t g_m=exp(-j*k*R_m)/(4.0*pi*R_m);
        complex_t g_p=exp(-j*k*R_p)/(4.0*pi*R_p);
        vector_t<real_t> R_m_u=unit(p-(basis_n.r_m+rho_m));
        vector_t<real_t> R_p_u=unit(p-(basis_n.r_p-rho_p));
        complex_t I_B_m, I_B_p;
        I_B_m = -((1.0+j*k*R_m)/R_m)*(direction*(R_m_u^rho_m))*g_m;
        I_B_p = -((1.0+j*k*R_p)/R_p)*(direction*(R_p_u^rho_p))*g_p;
        complex_t B;
        B = (I_B_m+I_B_p);
        sum+=basis_n.L*B*engine->I_n(n, 0);
    }
    return sum;
}

field_2d_t engine_2d_t::compute_near_field(const vector_t<real_t> p){
    this->shape.check();
    assert_error(is_I_n_available, "no I_n solutions found");
    integrand_2d_args args;
    args.engine = this;
    args.k = k;
    args.p = p;
    field_2d_t field;
    vector_t<real_t>x_u(1.0, 0.0, 0.0);
    vector_t<real_t>y_u(0.0, 1.0, 0.0);
    vector_t<real_t>z_u(0.0, 0.0, 1.0);
    //
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    this->quadl.set_2d(k_max, tol);
    int flag=false;
    args.direction = x_u;
    field.E.x = this->quadl.integral_2d(integrand_E_2d, &args, triangle, flag);
    if(flag){flag=false; print("warning: no convergence!\n");}
    field.H.x = this->quadl.integral_2d(integrand_H_2d, &args, triangle, flag);
    if(flag){flag=false; print("warning: no convergence!\n");}
    args.direction = y_u;
    field.E.y = this->quadl.integral_2d(integrand_E_2d, &args, triangle, flag);
    if(flag){flag=false; print("warning: no convergence!\n");}
    field.H.y = this->quadl.integral_2d(integrand_H_2d, &args, triangle, flag);
    if(flag){flag=false; print("warning: no convergence!\n");}
    args.direction = z_u;
    field.E.z = this->quadl.integral_2d(integrand_E_2d, &args, triangle, flag);
    if(flag){flag=false; print("warning: no convergence!\n");}
    field.H.z = this->quadl.integral_2d(integrand_H_2d, &args, triangle, flag);
    if(flag){flag=false; print("warning: no convergence!\n");}
    this->quadl.unset_2d();
    return field;
}

//

void engine_2d_t::save_Z_mn(const char *filename){
    binary_file_t file;
    assert(filename!=null);
    file.open(filename, 'w');
    file.write(&this->N);
    print("saving Z_mn solutions...");
    for (size_t m=0; m<this->N; m++){
        for (size_t n=0; n<this->N; n++){
            file.write(&this->Z_mn(m, n));
        }
    }
    file.close();
    print(", done\n");
}

void engine_2d_t::load_Z_mn(const char *filename){
    binary_file_t file;
    assert(filename!=null);
    file.open(filename, 'r');
    file.read(&this->N);
    this->Z_mn.unset();
    this->Z_mn.set(this->N, this->N);
    print("loading Z_mn solutions...");
    for (size_t m=0; m<this->N; m++){
        for (size_t n=0; n<this->N; n++){
            file.read(&this->Z_mn(m, n));
        }
    }
    file.close();
    print(", done\n");
    this->is_Z_mn_available = true;
}

// RCS

complex_t integrand_RCS_theta_2d(const complex_t alpha_n, const complex_t beta_n, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_n=args->basis_n;
    real_t theta_s=args->theta_s;
    real_t phi_s=args->phi_s;
    complex_t k=args->k;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> rho_m=+1.0*(real(alpha_n)*basis_n.L_m1+real(beta_n)*basis_n.L_m2);
    vector_t<real_t> rho_p=-1.0*(real(alpha_n)*basis_n.L_p1+real(beta_n)*basis_n.L_p2);
    vector_t<real_t> k_s(sin(theta_s)*cos(phi_s), 
                         sin(theta_s)*sin(phi_s), 
                         cos(theta_s));
    complex_t k_s_r_m=k*(k_s*(basis_n.r_m+rho_m));
    complex_t k_s_r_p=k*(k_s*(basis_n.r_p-rho_p));
    vector_t<real_t> theta_s_u(cos(theta_s)*cos(phi_s),
                               cos(theta_s)*sin(phi_s),
                               -sin(theta_s));
    //
    complex_t xi_m=basis_n.L*exp(j*k_s_r_m)*(rho_m*theta_s_u);
    complex_t xi_p=basis_n.L*exp(j*k_s_r_p)*(rho_p*theta_s_u);
    return xi_m+xi_p;
}

complex_t integrand_RCS_phi_2d(const complex_t alpha_n, const complex_t beta_n, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    basis_2d_t basis_n=args->basis_n;
    real_t theta_s=args->theta_s;
    real_t phi_s=args->phi_s;
    complex_t k=args->k;
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> rho_m=+1.0*(real(alpha_n)*basis_n.L_m1+real(beta_n)*basis_n.L_m2);
    vector_t<real_t> rho_p=-1.0*(real(alpha_n)*basis_n.L_p1+real(beta_n)*basis_n.L_p2);
    vector_t<real_t> k_s(sin(theta_s)*cos(phi_s), 
                         sin(theta_s)*sin(phi_s), 
                         cos(theta_s));
    complex_t k_s_r_m=k*(k_s*(basis_n.r_m+rho_m));
    complex_t k_s_r_p=k*(k_s*(basis_n.r_p-rho_p));
    vector_t<real_t> phi_s_u(-sin(phi_s), cos(phi_s), 0.0);
    //
    complex_t xi_m=basis_n.L*exp(j*k_s_r_m)*(rho_m*phi_s_u);
    complex_t xi_p=basis_n.L*exp(j*k_s_r_p)*(rho_p*phi_s_u);
    return xi_m+xi_p;
}

// engine

void  engine_2d_t::set_medium(const complex_t mu, const complex_t eps, const real_t freq){
    this->shape.set_medium(mu, eps, freq);
}

void engine_2d_t::mesh(const char *filename, const real_t clmax){
    assert_error(this->is_mesh_obtained==false, "mesh already exists");
    this->shape.mesh_2d(filename, clmax);
    this->shape.get_mesh();
    this->shape.get_basis_functions();
    this->shape.load_basis_functions();
    //
    shape_info_t shape_info=this->shape.get_shape_info();
    assert_error(shape_info.is_basis_2d_list_allocated, "no 2d elements were found");
    this->N = shape_info.N_2d_basis;
    this->k = shape_info.k;
    this->eta = shape_info.eta;
    this->freq = shape_info.freq;
    this->lambda = shape_info.lambda;
    this->Z_mn.set(N, N);
    this->V_m.set(N, 1);
    this->I_n.set(N, 1);
    this->is_mesh_obtained = true;
}

void engine_2d_t::solve_currents(){
    this->load_Z_mn("data/Z_mn.bin");
    assert_error(this->is_V_m_available, "no V_m solutions found");
    print("solving for I_n...");
    this->quadl.set_2d(k_max, tol);
    this->Z_mn.lup();
    this->I_n.set(this->N, 1);
    this->Z_mn.solve(this->V_m, this->I_n);
    this->quadl.unset_2d();
    this->is_I_n_available = true;
    print(", done\n");
}

void engine_2d_t::compute_Z_mn(){
    this->shape.check();
    this->quadl.set_2d(k_max, tol);
    basis_2d_t basis_m, basis_n;
    char *msg=(char*)calloc(line_max, sizeof(char));
    size_t count=0;
    timer_lib_t timer;
    timer.set();
    for (size_t m=0; m<1; m++){
        basis_m = this->shape.get_basis_2d(m);
        for (size_t n=m; n<N; n++){
            basis_n = this->shape.get_basis_2d(n);
            this->Z_mn(m, n) = Z_mn_2d(basis_m, basis_n, k, lambda, eta, this->quadl);
            count++;
        }
        for (size_t n=m+1; n<N; n++){
            basis_n = this->shape.get_basis_2d(n);
            this->Z_mn(n, m) = this->Z_mn(m, n);
        }
    }
    timer.unset_silent();
    real_t dt=timer.get_elapsed();
    dt = dt/N;
    dt = dt*N*(N-1.0)/2.0;
    if (dt<3600.0){
        print("estimated time: %3.1f mintues\n", dt/60.0);
    }else{
        print("estimated time: %3.1f hours\n", dt/3600.0);
    }
    for (size_t m=1; m<N; m++){
        basis_m = this->shape.get_basis_2d(m);
        for (size_t n=m; n<N; n++){
            basis_n = this->shape.get_basis_2d(n);
            sprintf(msg, "Z_mn (%zu, %zu)", m, n);
            progress_bar(count, N*(N+1)/2, msg);
            this->Z_mn(m, n) = Z_mn_2d(basis_m, basis_n, k, lambda, eta, this->quadl);
            count++;
        }
        for (size_t n=m+1; n<N; n++){
            basis_n = this->shape.get_basis_2d(n);
            this->Z_mn(n, m) = this->Z_mn(m, n);
        }
    }
    free(msg);
    this->quadl.unset_2d();
    engine_2d_t::save_Z_mn("data/Z_mn.bin");
    this->is_Z_mn_available = true;
}

void engine_2d_t::compute_V_m_plane_wave(const complex_t E_TM, const complex_t E_TE,
    const real_t theta_i, const real_t phi_i){
    this->shape.check();
    this->quadl.set_2d(k_max, tol);
    basis_2d_t basis_m;
    size_t count=0;
    for (size_t m=0; m<N; m++){
        basis_m = this->shape.get_basis_2d(m);
        this->V_m(m, 0) = V_m_plane_wave_2d(basis_m, E_TM, E_TE, k, theta_i, phi_i, this->quadl);
        count++;
    }
    this->quadl.unset_2d();
    this->is_V_m_available = true;
}

RCS_2d_t engine_2d_t::RCS_plane_wave_2d(const real_t theta_s, const real_t phi_s){
    this->shape.check();
    this->quadl.set_2d(k_max, tol);
    assert_error(is_I_n_available, "no I_n solutions found");
    int flag=false;
    integrand_2d_args args;
    args.k = k;
    args.theta_s = theta_s;
    args.phi_s = phi_s;
    RCS_2d_t RCS;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    complex_t sigma_theta=0.0;
    complex_t sigma_phi=0.0;
    for (size_t n=0; n<this->N; n++){
        args.basis_n = this->shape.get_basis_2d(n);
        sigma_theta+=this->quadl.integral_2d(integrand_RCS_theta_2d, &args, triangle, flag)*this->I_n(n, 0);
        if(flag){flag=false; print("warning: no convergence!\n");}
        sigma_phi+=this->quadl.integral_2d(integrand_RCS_phi_2d, &args, triangle, flag)*this->I_n(n, 0);
        if(flag){flag=false; print("warning: no convergence!\n");}
    }
    RCS.sigma_theta = pi*pow(abs(eta*sigma_theta), 2.0);
    RCS.sigma_phi = pi*pow(abs(eta*sigma_phi), 2.0);
    this->quadl.unset_2d();
    return RCS;
}

void engine_2d_t::unset(){
    this->N = 0;
    this->k = 0.0;
    this->eta = 0.0;
    this->Z_mn.unset();
    this->V_m.unset();
    this->I_n.unset();
    this->is_mesh_obtained = false;    
    this->is_Z_mn_available = false;
    this->is_V_m_available = false;
    this->is_I_n_available = false;
}

engine_2d_t::engine_2d_t(){
}

engine_2d_t::~engine_2d_t(){
}

field_2d_t incident_plane_wave_field(const real_t theta_i, const real_t phi_i,
    const complex_t k, const complex_t eta, const complex_t E_TM, const complex_t E_TE,
    const vector_t<real_t> p){
    vector_t<real_t> k_i(sin(theta_i)*cos(phi_i), 
                         sin(theta_i)*sin(phi_i), 
                         cos(theta_i));
    vector_t<real_t> theta_i_u(cos(theta_i)*cos(phi_i),
                               cos(theta_i)*sin(phi_i),
                               -sin(theta_i));
    vector_t<real_t> phi_i_u(-sin(phi_i), cos(phi_i), 0.0);
    field_2d_t field;
    vector_t<real_t>x_u(1.0, 0.0, 0.0);
    vector_t<real_t>y_u(0.0, 1.0, 0.0);
    vector_t<real_t>z_u(0.0, 0.0, 1.0);
    const complex_t j=complex_t(0.0, 1.0);
    complex_t k_i_r=k*(k_i*p);
    field.E.x = (E_TM*(theta_i_u*x_u)+E_TE*(phi_i_u*x_u))*exp(j*k_i_r);
    field.E.y = (E_TM*(theta_i_u*y_u)+E_TE*(phi_i_u*y_u))*exp(j*k_i_r);
    field.E.z = (E_TM*(theta_i_u*z_u)+E_TE*(phi_i_u*z_u))*exp(j*k_i_r);
    field.H.x = (E_TE*(theta_i_u*x_u)-E_TM*(phi_i_u*x_u))*exp(j*k_i_r)/eta;
    field.H.y = (E_TE*(theta_i_u*y_u)-E_TM*(phi_i_u*y_u))*exp(j*k_i_r)/eta;
    field.H.z = (E_TE*(theta_i_u*z_u)-E_TM*(phi_i_u*z_u))*exp(j*k_i_r)/eta;
    return field;
}

field_2d_t engine_2d_t::compute_incident_plane_wave_field(const real_t theta_i, const real_t phi_i,
    const complex_t E_TM, const complex_t E_TE, const vector_t<real_t> p){
    return incident_plane_wave_field(theta_i, phi_i, this->k, this->eta, E_TM, E_TE, p);
}

int is_inside_list(int *list, int index, const size_t N){
    for (size_t i=0; i<N; i++){
        if (index==list[i]){
            return true;
        }
    }
    return false;
}

void engine_2d_t::export_currents(){
    file_t file1, file2, file3;
    basis_2d_t basis_n;
    real_t alpha, beta;
    alpha = beta = 1.0/3.0;
    vector_t<complex_t> J_m, J_p;
    const real_t tol=1.0E-6;
    triangle_t *triangles_list=null, triangle, triangle_;
    //
    size_t N_triangles;
    size_t dummy_size_t;
    file1.open("mesh/mesh/mesh_data.txt", 'r');
    file1.read("%zu", &N_triangles);
    file1.read("%zu", &N_triangles);
    file1.read("%zu", &N_triangles);
    file1.read("%zu", &dummy_size_t);
    file1.close();
    triangles_list = (triangle_t*)calloc(N_triangles, sizeof(triangle_t));
    assert(triangles_list!=null);
    file2.open("mesh/mesh/mesh_2d.txt", 'r');
    int dummy1, dummy2, dummy3;
    for (size_t i=0; i<N_triangles; i++){
        file2.read("%lf %lf %lf %d\n", &triangles_list[i].v[0].x, &triangles_list[i].v[0].y, &triangles_list[i].v[0].z, &dummy1);
        file2.read("%lf %lf %lf %d\n", &triangles_list[i].v[1].x, &triangles_list[i].v[1].y, &triangles_list[i].v[1].z, &dummy2);
        file2.read("%lf %lf %lf %d\n", &triangles_list[i].v[2].x, &triangles_list[i].v[2].y, &triangles_list[i].v[2].z, &dummy3);
        file2.read("\n");
    }
    file2.close();
    //
    print("exporting currents...");
    file3.open("data/current.pos", 'w');
    file3.write("View \"background mesh\" {\n");
    for (size_t i=0; i<N_triangles; i++){
        triangle = triangles_list[i];
        vector_t<complex_t> J;  
        for (size_t n=0; n<N; n++){
            basis_n = engine_2d_t::shape.get_basis_2d(n);
            triangle_.v[0] = basis_n.r_m;
            triangle_.v[1] = basis_n.e_1;
            triangle_.v[2] = basis_n.e_2;
            size_t count=0;
            for (size_t ii=0; ii<3; ii++){
                for (size_t jj=0; jj<3; jj++){
                    if (is_equal(triangle.v[ii], triangle_.v[jj], tol)){
                        count++;
                    }
                    if (count==2){
                        J = J+(alpha*basis_n.L_m1+beta*basis_n.L_m2)*this->I_n(n, 0);
                        count = 0;
                        break;
                    }
                }
            }
            triangle_.v[0] = basis_n.r_p;
            triangle_.v[1] = basis_n.e_1;
            triangle_.v[2] = basis_n.e_2;
            for (size_t ii=0; ii<3; ii++){
                for (size_t jj=0; jj<3; jj++){
                    if (is_equal(triangle.v[ii], triangle_.v[jj], tol)){
                        count++;
                    }
                    if (count==2){
                        J = J-(alpha*basis_n.L_p1+beta*basis_n.L_p2)*this->I_n(n, 0);
                        count = 0;
                        break;
                    }
                }
            }
        }
        file3.write("ST(%21.14E, %21.14E, %21.14E, ", triangle.v[0].x, triangle.v[0].y, triangle.v[0].z);
        file3.write("%21.14E, %21.14E, %21.14E, ", triangle.v[1].x, triangle.v[1].y, triangle.v[1].z);
        file3.write("%21.14E, %21.14E, %21.14E){", triangle.v[2].x, triangle.v[2].y, triangle.v[2].z);
        file3.write("%21.14E, %21.14E, %21.14E};\n", mag(J), mag(J), mag(J));
    }
    file3.write("};\n");
    file3.close();
    print(", done!\n");
    free(triangles_list);
}