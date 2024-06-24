#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"

// Definitions

const real_t tol_vertex=1.0E-6;

struct edge_t{
    vector_t<real_t> v1, v2;
    int_t physical_group=0;
    real_t length=0.0;
    size_t N_adjacents=0;
    edge_t(){}
    edge_t(const vector_t<real_t> v1, const vector_t<real_t> v2, const int_t physical_group){
        this->v1 = v1;
        this->v2 = v2;
        get_length();
        this->physical_group = physical_group;
    }   
    void get_length(){
        this->length = mag(v1-v2);
    }
};

struct triangle_t{
    vector_t<real_t> v1, v2, v3;
    int_t physical_group=0;
    real_t area=0.0;
    vector_t<real_t> n;
    size_t N_adjacents=0;
    triangle_t(){}
    triangle_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const int_t physical_group){
        this->v1 = v1;
        this->v2 = v2;
        this->v3 = v3;
        get_area();
        this->physical_group = physical_group;
    }
    void get_area(){
        this->area = mag((v2-v1)^(v3-v1))/2.0;
        this->n = unit((v2-v1)^(v3-v1));
    }
};

struct tetrahedron_t{
    vector_t<real_t> v1, v2, v3, v4;
    int_t physical_group=0;
    real_t volume=0.0;
    complex_t eps=1.0;
    size_t N_adjacents=0;
    tetrahedron_t(){}
    tetrahedron_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const vector_t<real_t> v4, const int_t physical_group){
        this->v1 = v1;
        this->v2 = v2;
        this->v3 = v3;
        this->v4 = v4;
        get_volume();
        this->physical_group = physical_group;
    }
    void get_volume(){
        this->volume = ((v2-v1)^(v3-v1))*(v4-v1)/6.0;
    }
};

struct basis_2d_t{
    vector_t<real_t> r_m, r_p, e1, e2, n_m, n_p;
    real_t L, A_m, A_p;
    vector_t<real_t> L_m1, L_m2, L_p1, L_p2;
    int_t physical_group_m=0, physical_group_p=0;
    basis_2d_t(){}
    basis_2d_t(const vector_t<real_t> r_m, const vector_t<real_t> r_p, 
        const vector_t<real_t> e1, const vector_t<real_t> e2){
        this->r_m = r_m;
        this->r_p = r_p;
        this->e1 = e1;
        this->e2 = e2;
        basis_2d_t::get_values();
    }
    void get_values(){
        this->L_m1 = +1.0*(this->e1-this->r_m);
        this->L_m2 = +1.0*(this->e2-this->r_m);
        this->L_m1 = -1.0*(this->e1-this->r_p);
        this->L_m2 = -1.0*(this->e2-this->r_p);
        this->L = mag(this->e1-this->e2);
        this->A_m = mag((this->e1-this->r_m)^(this->e2-this->r_m))/2.0;
        this->A_p = mag((this->e2-this->r_p)^(this->e1-this->r_p))/2.0;
        this->n_m = unit((this->e1-this->r_m)^(this->e2-this->r_m));
        this->n_p = unit((this->e2-this->r_p)^(this->e1-this->r_p));
    }
};

class shape_t{
    private:
        size_t N_points=0, N_edges=0, N_triangles=0, N_tetrahedrons=0;
        triangle_t *triangle_data=null;
        tetrahedron_t *tetrahedron_data=null;
        int is_triangle_allocated=false;
        int is_tetrahedron_allocated=false;
        int is_basis_2d_list_allocated=false;
        void set();
        void unset();
        size_t N_1d_basis=0, N_2d_basis=0, N_3d_basis=0;
        basis_2d_t *basis_2d_list=null;
    public:
        shape_t();
        ~shape_t();
        void get_mesh();
        void load_mesh();
        void log_mesh();
        triangle_t get_triangle_element(const size_t index);
        tetrahedron_t get_tetrahedron_element(const size_t index);
        void assign_volume_properties(const complex_t eps, const int_t physical_group);
        void get_basis_functions();
};

// Functions


#endif
