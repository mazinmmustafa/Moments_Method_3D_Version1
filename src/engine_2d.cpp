//
#include "engine_2d.hpp"

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

struct integrand_2d_args{
    basis_2d_t basis_m, basis_n;
    complex_t k;
    quadl_domain_t quadl;
    complex_t alpha_m, beta_m;
};

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
    real_t l_m, l_p, R_m, R_p, R0, P0, d;
    vector_t<real_t> P0_u, u;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m=basis_m.r_m+real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2;
    vector_t<real_t> rho_p=basis_m.r_p+real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2;
    projection_2d_t para_mm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_m);
    projection_2d_t para_mp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_m);
    projection_2d_t para_pm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_p);
    projection_2d_t para_pp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_p);
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
    quadl_domain_t quadl, int &flag){
    integrand_2d_args args={basis_m, basis_n, k, quadl, 0, 0};
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
    real_t l_m, l_p, R_m, R_p, R0, P0, d;
    vector_t<real_t> P0_u, u;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m=basis_m.r_m+real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2;
    vector_t<real_t> rho_p=basis_m.r_p+real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2;
    projection_2d_t para_mm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_m);
    projection_2d_t para_mp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_m);
    projection_2d_t para_pm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_p);
    projection_2d_t para_pp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_p);
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
    real_t l_m, l_p, R_m, R_p, R0;
    vector_t<real_t> u;
    real_t I_mm, I_mp, I_pm, I_pp;
    vector_t<real_t> rho_m=basis_m.r_m+real(alpha_m)*basis_m.L_m1+real(beta_m)*basis_m.L_m2;
    vector_t<real_t> rho_p=basis_m.r_p+real(alpha_m)*basis_m.L_p1+real(beta_m)*basis_m.L_p2;
    projection_2d_t para_mm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_m);
    projection_2d_t para_mp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_m);
    projection_2d_t para_pm=get_projection_2d(basis_n.r_m, basis_n.e_1, basis_n.e_2, rho_p);
    projection_2d_t para_pp=get_projection_2d(basis_n.r_p, basis_n.e_2, basis_n.e_1, rho_p);
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
    quadl_domain_t quadl, int &flag){
    integrand_2d_args args={basis_m, basis_n, k, quadl, 0, 0};
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    return args.quadl.integral_2d(integrand_psi_2d_projection_1, &args, triangle, flag)+
            args.quadl.integral_2d(integrand_psi_2d_projection_2, &args, triangle, flag)+
            args.quadl.integral_2d(integrand_psi_2d_outer, &args, triangle, flag);
}

complex_t Z_mn_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    const complex_t eta, quadl_domain_t quadl){
    int flag=false;
    complex_t psi=psi_2d(basis_m, basis_n, k, quadl, flag); if(flag){flag=false; print("warning: no convergence!\n");}
    complex_t phi=+phi_2d(basis_m, basis_n, k, quadl, flag); if(flag){flag=false; print("warning: no convergence!\n");}
    const complex_t j=complex_t(0.0, 1.0);
    return j*k*eta*psi-j*(eta/k)*phi;
}

//

engine_2d_t::engine_2d_t(){

}

engine_2d_t::~engine_2d_t(){

}

void engine_2d_t::compute_Z_mn(){
    shape_info_t shape_info=this->shape.get_shape_info();
    assert_error(shape_info.is_basis_2d_list_allocated, "no 2d elements were found");
    size_t N=shape_info.N_2d_basis;
    complex_t k, eta;
    k = shape_info.k;
    eta = shape_info.eta;
    this->Z_mn.set(N, N);
    basis_2d_t basis_m, basis_n;
    size_t line_max=100;
    char *msg=(char*)calloc(line_max, sizeof(char));
    print("total number of basis functions: %zu\n", N);
    size_t count=0;
    for (size_t m=0; m<N; m++){
        basis_m = this->shape.get_basis_2d(m);
        for (size_t n=0; n<N; n++){
            basis_n = this->shape.get_basis_2d(n);
            sprintf(msg, "element (%zu, %zu)...", m, n);
            progress_bar(count, N*N, msg);
            this->Z_mn(m, n) = Z_mn_2d(basis_m, basis_n, k, eta, this->quadl);
            count++;
        }
    }
    free(msg);
}

// void engine_2d_t::save_Z_mn(const char *filename){

// }
