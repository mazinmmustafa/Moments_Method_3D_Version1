#ifndef __QUADL_HPP__
#define __QUADL_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"

// Definitions
class quadl_t{
    private:
        size_t N=0;
        size_t k_max=0;
        real_t *x=NULL, *w=NULL;
        int is_allocated=false;
        real_t tol=1.0E-4;
        //
        complex_t quadl_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const real_t a, const real_t b);
        complex_t quadl_1d_(complex_t (*func)(const complex_t, void*), 
            void *args, const real_t a, const real_t b, size_t &k, const complex_t I_p);
        complex_t quadl_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y);
        complex_t quadl_2d_(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
            size_t &k, const complex_t I_p);
        complex_t quadl_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
            const real_t a_z, const real_t b_z);
        complex_t quadl_3d_(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
            const real_t a_z, const real_t b_z, size_t &k, const complex_t I_p);
    public:
        quadl_t();
        ~quadl_t();
        void set(const size_t N, const size_t k_max, const real_t tol);
        void unset();
        void disp();
        complex_t integral_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const real_t a, const real_t b, int &flag);
        complex_t integral_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, int &flag);
        complex_t integral_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
            const real_t a_z, const real_t b_z, int &flag);
};

struct line_domain_t{ // domain x
    vector_t<real_t> v1, v2;
    real_t length(){
        return abs(v2.x-v1.x);
    }
};

struct triangle_domain_t{ // domain x, y
    vector_t<real_t> v1, v2, v3;
    real_t area(){
        return mag((v2-v1)^(v3-v1))/2.0;
    }
};

struct tetrahedral_domain_t{ // domain x, y, z
    vector_t<real_t> v1, v2, v3, v4;
    real_t volume(){
        return ((v2-v1)^(v3-v1))*(v4-v1)/6.0;
    }
};

class quadl_domain_t{
    private:
        size_t k_max=0;
        real_t tol=1.0E-4;
        //
        complex_t quadl_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const line_domain_t line);
        complex_t quadl_1d_(complex_t (*func)(const complex_t, void*), 
            void *args, line_domain_t line, size_t &k, const complex_t I_p);
        complex_t quadl_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, triangle_domain_t triangle);
        complex_t quadl_2d_(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const triangle_domain_t triangle, size_t &k, const complex_t I_p);
        complex_t quadl_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, tetrahedral_domain_t tetrahedron);
        complex_t quadl_3d_(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const tetrahedral_domain_t tetrahedron, size_t &k, const complex_t I_p);
    public:
        void set(const size_t k_max, const real_t tol){
            this->k_max = k_max;
            this->tol = tol;
        }
        complex_t integral_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const line_domain_t line, int &flag);
        complex_t integral_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const triangle_domain_t triangle, int &flag);
        complex_t integral_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const tetrahedral_domain_t tetrahedron, int &flag);
};

#endif
