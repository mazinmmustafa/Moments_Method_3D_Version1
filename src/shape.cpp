//
#include "shape.hpp"

shape_t::shape_t(){

}

shape_t::~shape_t(){
    shape_t::unset();
}

void shape_t::unset(){
    if (this->is_edge_allocated){
        free(this->edge_data);
    }
    if (this->is_triangle_allocated){
        free(this->triangle_data);
    }
    if (this->is_tetrahedron_allocated){
        free(this->tetrahedron_data);
    }
    this->N_edges = 0;
    this->N_triangles = 0;
    this->N_tetrahedrons = 0;
}

void shape_t::set(const size_t N_edges, const size_t N_triangles, const size_t N_tetrahedrons){
    this->N_edges = N_edges;
    this->N_triangles = N_triangles;
    this->N_tetrahedrons = N_tetrahedrons;
    shape_t::unset();
    this->edge_data = (edge_t*)calloc(this->N_edges, sizeof(edge_t));
    this->triangle_data = (triangle_t*)calloc(this->N_triangles, sizeof(triangle_t));
    this->tetrahedron_data = (tetrahedron_t*)calloc(this->N_tetrahedrons, sizeof(tetrahedron_t));
    assert(this->edge_data!=null);
    assert(this->triangle_data!=null);
    assert(this->tetrahedron_data!=null);
}

void shape_t::add_edge(const edge_t edge, const size_t index){
    assert_error(index<this->N_edges, "edge index is out of range");
    this->edge_data[index] = edge;
}

void shape_t::add_triangle(const triangle_t triangle, const size_t index){
    assert_error(index<this->N_triangles, "triangle index is out of range");
    this->triangle_data[index] = triangle;
}

void shape_t::add_tetrahedron(const tetrahedron_t tetrahedron, const size_t index){
    assert_error(index<this->N_tetrahedrons, "tetrahedron index is out of range");
    this->tetrahedron_data[index] = tetrahedron;
}

edge_t shape_t::get_edge(const size_t index){
    assert_error(index<this->N_edges, "edge index is out of range");
    return this->edge_data[index];
}

triangle_t shape_t::get_triangle(const size_t index){
    assert_error(index<this->N_triangles, "triangle index is out of range");
    return this->triangle_data[index];
}

tetrahedron_t shape_t::get_tetrahedron(const size_t index){
    assert_error(index<this->N_tetrahedrons, "tetrahedron index is out of range");
    return this->tetrahedron_data[index];
}
