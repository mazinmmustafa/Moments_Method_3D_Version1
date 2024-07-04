//
#include "testbench.hpp"

void test_utilities(){

    print("Hello, world!\n");
    print(complex_t(1.0, -2.0));
    print((real_t)3.2);
    print((int_t)45);
    print(c_0);
    print(eps_0);
    print(h_bar);
    print(pi);
    print(true);
    print(false);
    print(m_e);
    print(complex_t(-1.5, +0.7));
    assert_error(true, "this is an error message");

    const int Ns=1000;
    for (int i=0; i<Ns; i++){
        progress_bar(i, Ns, "doing something");
    }

    file_t file;
    file.open("data/test.txt", 'w');
    file.open("data/test.txt", 'w');
    file.write("This is something else\n");
    file.open("data/test.txt", 'w');
    file.write("This is something %d %d\n", 1, 2);
    file.write("This is something else\n");
    file.close();
    file.close();
    file.close();
    file.open("data/test.txt", 'a');
    file.write("This is something additional!\n");

    file.open("data/test.txt", 'r');

    char string[100];
    file.read("%s", string);
    print("%s\n", string);

    //
    const size_t Ns_=11;
    range_t range;
    
    range.set(-4.0, +6.0, Ns_);
    range.linspace();
    file.open("data/linspace.dat", 'w');
    for (size_t i=0; i<Ns_; i++){
        file.write("%21.14E\n", range(i));
    }
    file.close();

    range.set(1.0E-1, 1.0E+2, Ns_);
    range.logspace();
    file.open("data/lospace.dat", 'w');
    for (size_t i=0; i<Ns_; i++){
        file.write("%21.14E\n", range(i));
    }
    file.close();

}

void test_vector(){

    const complex_t j=complex_t(0.0, 1.0);
    vector_t<complex_t> A, B(1.0*j, -2.4, 1.7);
    A.x = 2.4-j;
    A.y = -0.4+1.2*j;
    A.z = 0.1;
    print(A*B);
    print(B/2);

    vector_t<real_t> a(1.2, -1.6, 0.2), b(0.1, 1.8, -2.1);
    print((a^b)+(0.4*a^unit(b)));

}

void test_read_write_binary_files(){

    const size_t Ns=10;
    vector_t<real_t> data;
    binary_file_t file;
    file.open("data/test.bin", 'w');
    file.write(&Ns);
    print("The written data are:\n");
    for (size_t i=0; i<Ns; i++){
        data.x = 1.0+i;
        data.y = 2.0+i;
        data.z = 3.0+i;
        print(data);
        file.write(&data);
    }
    file.close();

    size_t Ns_=0;
    file.open("data/test.bin", 'r');
    file.read(&Ns_);
    print("I have found %zu samples:\n", Ns_);
    for (size_t i=0; i<Ns_; i++){
        file.read(&data);
        print(data);
    }
    file.close();

}


void test_matrix(){

    matrix_t<real_t> A, B;
    A.set(3, 4);
    A.eye();
    A(0, 0) = -1.0;
    A.disp();

    B.copy(&A);
    B.disp();

    matrix_t<complex_t> C, D;
    C.set(100, 100);
    C.ones();
    C.save("data/matrix.bin");
    D.load("data/matrix.bin");
    D.disp();

    const size_t N=3;
    A.set(N, N);
    A(0, 0) = +0.0; A(0, 1) = +1.0; A(0, 2) = +2.0; 
    A(1, 0) = -4.0; A(1, 1) = +0.0; A(1, 2) = +0.0; 
    A(2, 0) = +6.0; A(2, 1) = +0.0; A(2, 2) = +1.0; 
    matrix_t<real_t> A_LUP;
    A_LUP.copy(&A);
    A_LUP.lup();
    print("\n\n");
    A_LUP.disp();

    matrix_t<real_t> b, x;
    x.set(3, 1);
    b.set(3, 1);
    b(0, 0) = +1.0; b(1, 0) = +1.0; b(2, 0) = +1.0;
    A_LUP.solve(b, x);
    print("\n\n");
    x.disp();
    print(A_LUP.det());

    A_LUP.inv();
    print("\n\n\n");
    A_LUP.disp();

    A.set(3, 3);
    B.set(3, 3);
    matrix_t<real_t> C_;
    C_.set(3, 3);
    A.ones();
    B.eye();

    add_matrix(A, B, C_);
    print("\n\n\n");
    C_.disp();
}

complex_t func_1d(const complex_t x, const real_t beta){
    return exp(-(beta*beta*x*x))+sin(x);
}; 

struct func_1d_args{
    real_t beta=0.0;
};

complex_t func_1d_wrapper(const complex_t x, void *args_){
    func_1d_args *args=(func_1d_args*)args_;
    return func_1d(x, args->beta);
}

complex_t func_2d(const complex_t x, const complex_t y, 
    const real_t alpha, const real_t beta){
    return cos(x+alpha)*log(beta-y);
}; 

struct func_2d_args{
    real_t alpha=0.0;
    real_t beta=0.0;
};

complex_t func_2d_wrapper(const complex_t x, const complex_t y, void *args_){
    func_2d_args *args=(func_2d_args*)args_;
    return func_2d(x, y, args->alpha, args->beta);
}

complex_t func_3d(const complex_t x, const complex_t y, const complex_t z, void *args){
    assert(args==null);
    return x+0.0*(y+z);
}; 

complex_t func_1d_line(const complex_t x, void *args){
    assert(args==null);
    return x-x*x+x*x*x;
}; 

complex_t func_2d_triangle(const complex_t x, const complex_t y, void *args){
    assert(args==null);
    return x-x*x+y*y-y*y*y;
}; 

complex_t func_3d_tetrahedron(const complex_t x, const complex_t y, const complex_t z, void *args){
    assert(args==null);
    return x-x*x+y*y+z*z*z;
}; 

void test_quadl(){

    // quadl_t quadl;
    // const size_t N_quadl=16;
    // const size_t k_max=15;
    // const real_t tol=1.0E-10;

    // quadl.set(N_quadl, k_max, tol);
    // quadl.disp();
    // print("\n\n");

    // int flag;

    // real_t beta=10.0;
    // func_1d_args args_1d={beta};
    // // +4.14742169407021E-01
    // print(quadl.integral_1d(func_1d_wrapper, &args_1d, -2.0, +4.0, flag));
    // if(flag){print("I_1d: no convergence!\n");}
    
    // real_t alpha=0.2; 
    // beta = 2.0;
    // func_2d_args args_2d={alpha, beta};
    // // +2.13734707679366E+00
    // print(quadl.integral_2d(func_2d_wrapper, &args_2d, -1.0, +1.0, -1.0, +1.0, flag));
    // if(flag){print("I_2d: no convergence!\n");}

    // // +3.14159265358979E+00
    // print(quadl.integral_3d(func_3d, null, 0.0, 1.0, 0.0, 2.0*pi, -0.5, +0.5, flag));  
    // if(flag){print("I_3d: no convergence!\n");}

    // //
    // quadl_domain_t quadl_domain;
    // quadl_domain.set(2500, 1.0E-4);

    // line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
    // print(quadl_domain.integral_1d(func_1d_line, null, line, flag));
    // if(flag){print("I_1d: no convergence!\n");}

    // triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), 
    // vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    // print(quadl_domain.integral_2d(func_2d_triangle, null, triangle, flag));
    // if(flag){print("I_2d: no convergence!\n");}

    // tetrahedral_domain_t tetrahedron={
    //     vector_t<real_t>(0.0, 0.0, 0.0),
    //     vector_t<real_t>(1.0, 0.0, 0.0),
    //     vector_t<real_t>(0.0, 1.0, 0.0),
    //     vector_t<real_t>(0.0, 0.0, 1.0)
    // };
    // print(quadl_domain.integral_3d(func_3d_tetrahedron, null, tetrahedron, flag));
    // if(flag){print("I_3d: no convergence!\n");}
    
}

void test_shape(){

    shape_t shape;
    shape.get_mesh();
    
    // enum physical_groups{patch=1, ground=2, port=3, substrate=4};  
    // const complex_t eps_substrate=4.3;
    // shape.assign_volume_properties(eps_substrate, substrate);

    shape.get_basis_functions();

}

complex_t dummy(const complex_t z){
    return 2.0*z;
}

void test_engine_2d(){

    engine_2d_t engine;
    timer_lib_t timer;

    quadl_domain_t quadl;
    const size_t k_max=10;
    const real_t tol=1.0E-4;

    quadl.set_2d(k_max, tol);
    const complex_t k=2.0*pi;
    const complex_t eta=sqrt(mu_0/eps_0);

    vector_t<real_t> v1_m, v2_m, v3_m, v4_m;
    vector_t<real_t> v1_n, v2_n, v3_n, v4_n;

    // Scenario I

    v1_m.x = +0.3/10.0; 
    v1_m.y = -0.4/10.0; 
    v1_m.z = +0.3/10.0;

    v2_m.x = +0.0/10.0; 
    v2_m.y = +0.2/10.0; 
    v2_m.z = +0.0/10.0;

    v3_m.x = +0.1/10.0; 
    v3_m.y = +0.0/10.0; 
    v3_m.z = -0.7/10.0;

    v4_m.x = -0.3/10.0; 
    v4_m.y = +0.0/10.0; 
    v4_m.z = +0.0/10.0;

    
    v1_n.x = +0.0; 
    v1_n.y = -0.4; 
    v1_n.z = +0.0;

    v2_n.x = +0.0; 
    v2_n.y = +0.2; 
    v2_n.z = +0.0;

    v3_n.x = +0.1; 
    v3_n.y = +0.0; 
    v3_n.z = +0.0;

    v4_n.x = -0.3; 
    v4_n.y = +0.0; 
    v4_n.z = +0.0;

    basis_2d_t basis_m(v1_m, v2_m, v3_m, v4_m), basis_n(v1_m, v2_m, v3_m, v4_m);

    // Scenario II Case I

    // v1_m.x = +0.0; 
    // v1_m.y = -0.1; 
    // v1_m.z = +0.0;

    // v2_m.x = +0.0; 
    // v2_m.y = +0.2; 
    // v2_m.z = +0.0;

    // v3_m.x = +0.3; 
    // v3_m.y = +0.0; 
    // v3_m.z = +0.0;

    // v4_m.x = -0.2; 
    // v4_m.y = +0.0; 
    // v4_m.z = +0.0;

    
    // v1_n.x = +0.2; 
    // v1_n.y = +0.3; 
    // v1_n.z = +0.0;

    // v2_n.x = -0.2; 
    // v2_n.y = +0.0; 
    // v2_n.z = +0.0;

    // v3_n.x = +0.0; 
    // v3_n.y = +0.2; 
    // v3_n.z = +0.0;

    // v4_n.x = +0.3; 
    // v4_n.y = +0.0; 
    // v4_n.z = +0.0;

    // Scenario II Case II

    // v1_m.x = +0.0; 
    // v1_m.y = -0.1; 
    // v1_m.z = +0.0;

    // v2_m.x = +0.0; 
    // v2_m.y = +0.2; 
    // v2_m.z = +0.0;

    // v3_m.x = +0.3; 
    // v3_m.y = +0.0; 
    // v3_m.z = +0.0;

    // v4_m.x = -0.2; 
    // v4_m.y = +0.0; 
    // v4_m.z = +0.0;

    
    // v1_n.x = -0.1; 
    // v1_n.y = +0.2; 
    // v1_n.z = +0.0;

    // v2_n.x = +0.3; 
    // v2_n.y = +0.0; 
    // v2_n.z = +0.0;

    // v3_n.x = -0.2; 
    // v3_n.y = +0.0; 
    // v3_n.z = +0.0;

    // v4_n.x = +0.0; 
    // v4_n.y = +0.2; 
    // v4_n.z = +0.0;


    // v1_m.x = +4.58016833999998E-02; 
    // v1_m.y = +8.06880706999997E-02; 
    // v1_m.z = -4.06501183599998E-01;

    // v2_m.x = +2.52743795999999E-02; 
    // v2_m.y = -7.89161683999997E-02; 
    // v2_m.z = -4.08638000399998E-01;

    // v3_m.x = +0.00000000000000E+00; 
    // v3_m.y = +0.00000000000000E+00; 
    // v3_m.z = -4.16955118999998E-01;

    // v4_m.x = -8.29421324999997E-02; 
    // v4_m.y = +6.85692369999997E-03; 
    // v4_m.z = -4.08564751899998E-01;
    
    
    // basis_2d_t basis_m(v1_m, v2_m, v3_m, v4_m), basis_n(v1_m, v2_m, v3_m, v4_m);

    timer.set();
    print(Z_mn_2d(basis_m, basis_n, k, eta, quadl));
    timer.unset();

    quadl.unset_2d();

}


void test_Z_mn_2d(){

    timer_lib_t timer;

    const real_t GHz=1.0E+9;
    real_t freq=0.25*GHz;
    real_t mu=1.0, eps=1.0;
    real_t lambda=c_0/freq;

    engine_2d_t engine;  
    engine.shape.set_medium(mu, eps, freq);

    const real_t clmax=0.1*lambda;
    engine.shape.mesh_2d("FreeCAD/test_sphere.geo", clmax);

    engine.shape.get_mesh();
    engine.shape.get_basis_functions();
    engine.shape.load_basis_functions();
    
    const size_t k_max=25;
    const real_t tol=1.0E-3;
    engine.quadl.set_2d(k_max, tol);

    // find Z_mn
    // timer.set();
    // engine.compute_Z_mn();
    // timer.unset();
    // engine.save_Z_mn("data/Z_mn.bin");

    real_t theta_i, phi_i;
    complex_t E_TM, E_TE;

    RCS_2d RCS;
    file_t file;
    real_t phi_s;
    real_t theta_s_min;
    real_t theta_s_max;
    range_t theta_s;
    const size_t Ns=401;

    // RCS theta-theta
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    E_TM = 1.0;
    E_TE = 0.0;

    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);

    print("solving for I_n...");
    engine.load_Z_mn("data/Z_mn.bin");
    engine.Z_mn.lup();
    engine.I_n.set(engine.N, 1);
    engine.Z_mn.solve(engine.V_m, engine.I_n);
    print(", done!\n");

    phi_s = deg2rad(0.0);
    theta_s_min = deg2rad(-180.0);
    theta_s_max = deg2rad(+180.0);
    theta_s.set(theta_s_min, theta_s_max, Ns);
    theta_s.linspace();
    
    file.open("data/RCS1.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing RCS...");
        RCS = engine.RCS_plane_wave_2d(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", theta_s(i), RCS.sigma_theta, 
            RCS.sigma_phi);
    }
    file.close();

    // RCS phi-phi
    theta_i = deg2rad(0.0);
    phi_i = deg2rad(0.0);
    E_TM = 0.0;
    E_TE = 1.0;

    engine.compute_V_m_plane_wave(E_TM, E_TE, theta_i, phi_i);

    print("solving for I_n...");
    engine.load_Z_mn("data/Z_mn.bin");
    engine.Z_mn.lup();
    engine.I_n.set(engine.N, 1);
    engine.Z_mn.solve(engine.V_m, engine.I_n);
    print(", done!\n");

    phi_s = deg2rad(0.0);
    theta_s_min = deg2rad(-180.0);
    theta_s_max = deg2rad(+180.0);
    theta_s.set(theta_s_min, theta_s_max, Ns);
    theta_s.linspace();
    
    file.open("data/RCS2.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing RCS...");
        RCS = engine.RCS_plane_wave_2d(theta_s(i), phi_s);
        file.write("%21.14E %21.14E %21.14E\n", theta_s(i), RCS.sigma_theta, 
            RCS.sigma_phi);
    }
    file.close();


    engine.quadl.unset_2d();

}






