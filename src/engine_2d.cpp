//
#include "engine_2d.hpp"

void R_mn_2d(const real_t alpha_m, const real_t beta_m, const real_t alpha_n, const real_t beta_n, 
    const basis_2d_t &basis_m, const basis_2d_t &basis_n, 
    real_t &R_mn_mm, real_t &R_mn_mp, real_t &R_mn_pm, real_t &R_mn_pp){
    vector_t<real_t> R_m_m = +1.0*(basis_m.r_m+alpha_m*basis_m.L_m1+beta_m*basis_m.L_m2);
    vector_t<real_t> R_m_p = +1.0*(basis_m.r_p+alpha_m*basis_m.L_p2+beta_m*basis_m.L_p1);
    vector_t<real_t> R_n_m = +1.0*(basis_n.r_m+alpha_n*basis_n.L_m1+beta_n*basis_n.L_m2);
    vector_t<real_t> R_n_p = +1.0*(basis_n.r_p+alpha_n*basis_n.L_p2+beta_n*basis_n.L_p1);
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

complex_t integrand_phi_2d_inner_mm(const complex_t alpha_n, const complex_t beta_n, void *args_){
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
    return (basis_m.L*basis_n.L/pi)*(+I_mm);
}

complex_t integrand_phi_2d_inner_mp(const complex_t alpha_n, const complex_t beta_n, void *args_){
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
    I_mp = -j*k*exp(-j*k*R_mp/2.0)*sinc(k*R_mp/2.0);
    return (basis_m.L*basis_n.L/pi)*(-I_mp);
}

complex_t integrand_phi_2d_inner_pm(const complex_t alpha_n, const complex_t beta_n, void *args_){
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
    I_pm = -j*k*exp(-j*k*R_pm/2.0)*sinc(k*R_pm/2.0);
    return (basis_m.L*basis_n.L/pi)*(-I_pm);
}

complex_t integrand_phi_2d_inner_pp(const complex_t alpha_n, const complex_t beta_n, void *args_){
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
    I_pp = -j*k*exp(-j*k*R_pp/2.0)*sinc(k*R_pp/2.0);
    return (basis_m.L*basis_n.L/pi)*(I_pp);
}

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

complex_t integrand_phi_2d_outer_mm(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_phi_2d_inner_mm, args, triangle, flag);
}

complex_t integrand_phi_2d_outer_mp(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_phi_2d_inner_mp, args, triangle, flag);
}

complex_t integrand_phi_2d_outer_pm(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_phi_2d_inner_pm, args, triangle, flag);
}

complex_t integrand_phi_2d_outer_pp(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_phi_2d_inner_pp, args, triangle, flag);
}

complex_t integrand_phi_2d_outer(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_phi_2d_inner, args, triangle, flag);
}

complex_t phi_2d(const basis_2d_t basis_m, const basis_2d_t basis_n, const complex_t k, 
    quadl_domain_t quadl, int &flag){
    integrand_2d_args args={basis_m, basis_n, k, quadl, 0, 0};
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    // return args.quadl.integral_2d(integrand_phi_2d_outer_mm, &args, triangle, flag)+
    //        args.quadl.integral_2d(integrand_phi_2d_outer_mp, &args, triangle, flag)+
    //        args.quadl.integral_2d(integrand_phi_2d_outer_pm, &args, triangle, flag)+
    //        args.quadl.integral_2d(integrand_phi_2d_outer_pp, &args, triangle, flag);
    return args.quadl.integral_2d(integrand_phi_2d_outer, &args, triangle, flag);
}

complex_t integrand_inner(const complex_t alpha_n, const complex_t beta_n, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    real_t alpha_m=real(args->alpha_m);
    real_t beta_m=real(args->beta_m);
    return cos(alpha_n+beta_m)*cos(beta_n)*alpha_m;
}

complex_t integrand_outer(const complex_t alpha_m, const complex_t beta_m, void *args_){
    integrand_2d_args *args=(integrand_2d_args*)args_;
    args->alpha_m = alpha_m;
    args->beta_m = beta_m;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(2.0, 0.0, 0.0), vector_t<real_t>(0.0, 2.0, 0.0)};
    int flag;
    return args->quadl.integral_2d(integrand_inner, args, triangle, flag);
}

complex_t get_integral(quadl_domain_t quadl, int &flag){
    integrand_2d_args args;
    args.quadl = quadl;
    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    return args.quadl.integral_2d(integrand_outer, &args, triangle, flag);
}