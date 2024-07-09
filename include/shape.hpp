#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"

// Definitions

const size_t __max_system__=20000;

struct edge_t{
    vector_t<real_t> v[2];
    int_t physical_group=0;
    real_t length=0.0;
    size_t N_adjacents=0;
    edge_t(){}
    edge_t(const vector_t<real_t> v1, const vector_t<real_t> v2, const int_t physical_group){
        this->v[0] = v1;
        this->v[1] = v2;
        get_length();
        this->physical_group = physical_group;
    }   
    void get_length(){
        this->length = mag(v[0]-v[1]);
    }
};

struct triangle_t{
    vector_t<real_t> v[3];
    int_t physical_group=0;
    real_t area=0.0;
    vector_t<real_t> n;
    size_t N_adjacents=0;
    triangle_t(){}
    triangle_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const int_t physical_group){
        this->v[0] = v1;
        this->v[1] = v2;
        this->v[2] = v3;
        get_area();
        this->physical_group = physical_group;
    }
    void get_area(){
        vector_t<real_t> vector=(v[1]-v[0])^(v[2]-v[0]);
        this->area = mag(vector)/2.0;
        this->n = unit(vector);
    }
};

struct tetrahedron_t{
    vector_t<real_t> v[4];
    int_t physical_group=0;
    real_t volume=0.0;
    complex_t eps=complex_t(+1.0, -0.0);
    complex_t mu=complex_t(+1.0, -0.0);
    size_t N_adjacents=0;
    tetrahedron_t(){}
    tetrahedron_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const vector_t<real_t> v4, const int_t physical_group){
        this->v[0] = v1;
        this->v[1] = v2;
        this->v[2] = v3;
        this->v[3] = v4;
        get_volume();
        this->physical_group = physical_group;
    }
    void get_volume(){
        this->volume = ((v[1]-v[0])^(v[2]-v[0]))*(v[3]-v[0])/6.0;
    }
};

struct basis_1d_t{
    vector_t<real_t> r_m, r_p, e_1;
    vector_t<real_t> L_m, L_p;
    int_t physical_group_m=0, physical_group_p=0;
    basis_1d_t(){}
    basis_1d_t(const vector_t<real_t> r_m, const vector_t<real_t> r_p, const vector_t<real_t> e_1){
        this->r_m = r_m;
        this->r_p = r_p;
        this->e_1 = e_1;
        basis_1d_t::get_values();
    }
    void get_values(){
        this->L_m = +1.0*(e_1-r_m);
        this->L_p = -1.0*(e_1-r_p);
    }
};

struct basis_2d_t{
    vector_t<real_t> r_m, r_p, e_1, e_2, n_m, n_p;
    real_t L, A_m, A_p;
    vector_t<real_t> L_m1, L_m2, L_p1, L_p2;
    int_t physical_group_m=0, physical_group_p=0;
    basis_2d_t(){}
    basis_2d_t(const vector_t<real_t> r_m, const vector_t<real_t> r_p, 
        const vector_t<real_t> e_1, const vector_t<real_t> e_2){
        this->r_m = r_m;
        this->r_p = r_p;
        this->e_1 = e_1;
        this->e_2 = e_2;
        basis_2d_t::get_values();
    }
    void get_values(){
        this->L_m1 = +1.0*(e_1-r_m);
        this->L_m2 = +1.0*(e_2-r_m);
        this->L_p1 = +1.0*(e_1-r_p);
        this->L_p2 = +1.0*(e_2-r_p);
        this->L = mag(e_1-e_2);
        this->A_m = mag(L_m1^L_m2)/2.0;
        this->A_p = mag(L_p2^L_p1)/2.0;
        this->n_m = unit(L_m1^L_m2);
        this->n_p = unit(L_p2^L_p1);
    }
};

struct basis_3d_t{
    vector_t<real_t> r_m, r_p, e_1, e_2, e_3;
    vector_t<real_t> n, n_m1, n_m2, n_m3, n_p1, n_p2, n_p3;
    real_t A, V_m, V_p;
    vector_t<real_t> L_m1, L_m2, L_m3, L_p1, L_p2, L_p3;
    int_t physical_group_m=0, physical_group_p=0;
    basis_3d_t(){}
    basis_3d_t(const vector_t<real_t> r_m, const vector_t<real_t> r_p, 
        const vector_t<real_t> e_1, const vector_t<real_t> e_2, const vector_t<real_t> e_3){
        this->r_m = r_m;
        this->r_p = r_p;
        this->e_1 = e_1;
        this->e_2 = e_2;
        this->e_3 = e_3;
        basis_3d_t::get_values();
    }
    void get_values(){
        this->L_m1 = +1.0*(this->e_1-this->r_m);
        this->L_m2 = +1.0*(this->e_2-this->r_m);
        this->L_m3 = +1.0*(this->e_3-this->r_m);
        this->L_p1 = -1.0*(this->e_1-this->r_p);
        this->L_p2 = -1.0*(this->e_2-this->r_p);
        this->L_p3 = -1.0*(this->e_3-this->r_p);
        vector_t<real_t> vector;
        vector=(e_2-e_1)^(e_3-e_1);
        this->A = mag(vector)/2.0;
        this->V_m = +1.0*(L_m1^L_m2)*L_m3/6.0;
        this->V_p = -1.0*(L_p1^L_p2)*L_p3/6.0;
    }
};

struct shape_info_t{
    size_t N_1d_basis, N_2d_basis, N_3d_basis;
    int is_basis_1d_list_allocated; 
    int is_basis_2d_list_allocated;  
    int is_basis_3d_list_allocated;
    complex_t k, eta;
    real_t freq, lambda;
};

class shape_t{
    private:
        size_t N_points=0, N_edges=0, N_triangles=0, N_tetrahedrons=0;
        edge_t *edge_data=null;
        triangle_t *triangle_data=null;
        tetrahedron_t *tetrahedron_data=null;
        int is_edge_allocated=false;
        int is_triangle_allocated=false;
        int is_tetrahedron_allocated=false;
        int is_basis_1d_list_allocated=false;
        int is_basis_2d_list_allocated=false;
        int is_basis_3d_list_allocated=false;
        void set();
        size_t N_1d_basis=0, N_2d_basis=0, N_3d_basis=0;
        basis_1d_t *basis_1d_list=null;
        basis_2d_t *basis_2d_list=null;
        basis_3d_t *basis_3d_list=null;
        void free_basic_elements();
        complex_t mu, eps;
        complex_t k, eta;
        real_t freq, lambda;
        int is_medium_set=false;
    public:
        void unset();
        shape_t();
        ~shape_t();
        void get_mesh();
        void load_mesh();
        void log_mesh();
        edge_t get_edge_element(const size_t index);
        triangle_t get_triangle_element(const size_t index);
        tetrahedron_t get_tetrahedron_element(const size_t index);
        basis_1d_t get_basis_1d(const size_t index);
        basis_2d_t get_basis_2d(const size_t index);
        basis_3d_t get_basis_3d(const size_t index);
        void assign_volume_properties(const complex_t eps, const int_t physical_group);
        void get_basis_functions();
        void load_basis_functions();
        shape_info_t get_shape_info();
        void set_medium(const complex_t mu, const complex_t eps, const real_t freq);
        void mesh_1d(const char *gmsh_filename, const real_t clmax);
        void mesh_2d(const char *gmsh_filename, const real_t clmax);
        void mesh_3d(const char *gmsh_filename, const real_t clmax);
        void check();
};

// Functions


#endif
