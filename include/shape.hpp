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
        this->length = mag(v1-v2);
        this->physical_group = physical_group;
    }   
};

struct triangle_t{
    edge_t edge_12, edge_13, edge_23;
    size_t physical_group=0;
    real_t area=0.0;
    vector_t<real_t> n;
    triangle_t(){}
    triangle_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const size_t physical_group){
        this->edge_12 = edge_t(v1, v2, physical_group);
        this->edge_13 = edge_t(v1, v3, physical_group);
        this->edge_23 = edge_t(v2, v3, physical_group);
        this->area = mag((v2-v1)^(v3-v1))/2.0;
        this->n = unit((v2-v1)^(v3-v1));
        this->physical_group = physical_group;
    }
};

struct tetrahedron_t{
    triangle_t triangle_123, triangle_124, triangle_134, triangle_234;
    size_t physical_group=0;
    real_t volume=0.0;
    tetrahedron_t(){}
    tetrahedron_t(const vector_t<real_t> v1, const vector_t<real_t> v2, 
        const vector_t<real_t> v3, const vector_t<real_t> v4, const size_t physical_group){
        this->triangle_123 = triangle_t(v1, v3, v2, physical_group);
        this->triangle_124 = triangle_t(v1, v2, v4, physical_group);
        this->triangle_134 = triangle_t(v1, v4, v3, physical_group);
        this->triangle_234 = triangle_t(v2, v3, v4, physical_group);
        this->volume = ((v2-v1)^(v3-v1))*(v4-v1)/6.0;
        this->physical_group = physical_group;
    }
};

class shape_t{
    private:
        size_t N_edges=0, N_triangles=0, N_tetrahedrons=0;
        edge_t *edge_data=null;
        triangle_t *triangle_data=null;
        tetrahedron_t *tetrahedron_data=null;
        int is_edge_allocated=false;
        int is_triangle_allocated=false;
        int is_tetrahedron_allocated=false;
        void unset();
    public:
        shape_t();
        ~shape_t();
        void set(const size_t N_edges, const size_t N_triangles, const size_t N_tetrahedrons);
        void add_edge(const edge_t edge, const size_t index);
        void add_triangle(const triangle_t triangle, const size_t index);
        void add_tetrahedron(const tetrahedron_t tetrahedron, const size_t index);
        edge_t get_edge(const size_t index);
        triangle_t get_triangle(const size_t index);
        tetrahedron_t get_tetrahedron(const size_t index);
};

// Functions


#endif
