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

    quadl_t quadl;
    const size_t N_quadl=16;
    const size_t k_max=15;
    const real_t tol=1.0E-10;

    quadl.set(N_quadl, k_max, tol);
    quadl.disp();
    print("\n\n");

    int flag;

    real_t beta=10.0;
    func_1d_args args_1d={beta};
    // +4.14742169407021E-01
    print(quadl.integral_1d(func_1d_wrapper, &args_1d, -2.0, +4.0, flag));
    if(flag){print("I_1d: no convergence!\n");}
    
    real_t alpha=0.2; 
    beta = 2.0;
    func_2d_args args_2d={alpha, beta};
    // +2.13734707679366E+00
    print(quadl.integral_2d(func_2d_wrapper, &args_2d, -1.0, +1.0, -1.0, +1.0, flag));
    if(flag){print("I_2d: no convergence!\n");}

    // +3.14159265358979E+00
    print(quadl.integral_3d(func_3d, null, 0.0, 1.0, 0.0, 2.0*pi, -0.5, +0.5, flag));  
    if(flag){print("I_3d: no convergence!\n");}

    //
    quadl_domain_t quadl_domain;
    quadl_domain.set(2500, 1.0E-4);

    line_domain_t line={vector_t<real_t>(0.0, 0.0, 0.0), vector_t<real_t>(1.0, 0.0, 0.0)};
    print(quadl_domain.integral_1d(func_1d_line, null, line, flag));
    if(flag){print("I_1d: no convergence!\n");}

    triangle_domain_t triangle={vector_t<real_t>(0.0, 0.0, 0.0), 
    vector_t<real_t>(1.0, 0.0, 0.0), vector_t<real_t>(0.0, 1.0, 0.0)};
    print(quadl_domain.integral_2d(func_2d_triangle, null, triangle, flag));
    if(flag){print("I_2d: no convergence!\n");}

    tetrahedral_domain_t tetrahedron={
        vector_t<real_t>(0.0, 0.0, 0.0),
        vector_t<real_t>(1.0, 0.0, 0.0),
        vector_t<real_t>(0.0, 1.0, 0.0),
        vector_t<real_t>(0.0, 0.0, 1.0)
    };
    print(quadl_domain.integral_3d(func_3d_tetrahedron, null, tetrahedron, flag));
    if(flag){print("I_3d: no convergence!\n");}
    
}

void test_shape(){

    shape_t shape;
    shape.get_mesh();
    shape.load_mesh();
    shape.log_mesh();

    // enum physical_groups{patch=1, ground=2, port=3, substrate=4};  
    // const complex_t eps_substrate=4.3;
    // shape.assign_volume_properties(eps_substrate, substrate);

    shape.get_basis_functions();
    shape.load_basis_functions();

}






