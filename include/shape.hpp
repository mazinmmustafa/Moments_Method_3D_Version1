#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"

// Definitions

struct edge_t{
    vector_t<real_t> v1, v2;
    size_t physical_group=0;
    real_t length=0.0;
    edge_t(){}
    edge_t(const vector_t<real_t> v1, const vector_t<real_t> v2, const size_t physical_group){
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
    vector_t<real_t> v1, v2, v3, v4;
    size_t physical_group=0;
    real_t area=0.0;
    vector_t<real_t> n;
    triangle_t(){}
    triangle_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const size_t physical_group){
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
    tetrahedron_t(){}
    tetrahedron_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const vector_t<real_t> v4, const size_t physical_group){
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

class shape_t{
    private:
        size_t N_points=0, N_edges=0, N_triangles=0, N_tetrahedrons=0;
        triangle_t *triangle_data=null;
        tetrahedron_t *tetrahedron_data=null;
        int is_triangle_allocated=false;
        int is_tetrahedron_allocated=false;
        void set();
        void unset();
    public:
        shape_t();
        ~shape_t();
        void load_mesh();
        void log_mesh();
        triangle_t get_triangle_element(const size_t index);
        tetrahedron_t get_tetrahedron_element(const size_t index);
};

// Functions


#endif
