//
#include "shape.hpp"

shape_t::shape_t(){

}

shape_t::~shape_t(){
    shape_t::unset();
}

void shape_t::unset(){
    if (this->is_triangle_allocated){
        free(this->triangle_data);
    }
    if (this->is_tetrahedron_allocated){
        free(this->tetrahedron_data);
    }
    this->N_points = 0;
    this->N_edges = 0;
    this->N_triangles = 0;
    this->N_tetrahedrons = 0;
    this->is_triangle_allocated = false;
    this->is_tetrahedron_allocated = false;
}

void shape_t::set(){
    this->triangle_data = (triangle_t*)calloc(this->N_triangles, sizeof(triangle_t));
    this->tetrahedron_data = (tetrahedron_t*)calloc(this->N_tetrahedrons, sizeof(tetrahedron_t));
    assert(this->triangle_data!=null);
    assert(this->tetrahedron_data!=null);
    this->is_triangle_allocated = true;
    this->is_tetrahedron_allocated = true;
}

void shape_t::log_mesh(){
    file_t file;
    //
    file.open("mesh/mesh/mesh_2d_elements_log.txt", 'w');
    file.write("number of triangles: %zu\n\n", this->N_triangles);
    for (size_t i=0; i<this->N_triangles; i++){
        file.write("triangle: %zu: group: %d\n", i, 
            this->triangle_data[i].physical_group);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->triangle_data[i].v1.x,
            this->triangle_data[i].v1.y,
            this->triangle_data[i].v1.z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->triangle_data[i].v2.x,
            this->triangle_data[i].v2.y,
            this->triangle_data[i].v2.z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->triangle_data[i].v3.x,
            this->triangle_data[i].v3.y,
            this->triangle_data[i].v3.z);
        file.write("\n");
    }
    file.close();
    //
    file.open("mesh/mesh/mesh_3d_elements_log.txt", 'w');
    file.write("number of tetrahedrons: %zu\n\n", this->N_tetrahedrons);
    for (size_t i=0; i<this->N_tetrahedrons; i++){
        file.write("tetrahedron: %zu: group: %d\n", i, 
            this->tetrahedron_data[i].physical_group);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->tetrahedron_data[i].v1.x,
            this->tetrahedron_data[i].v1.y,
            this->tetrahedron_data[i].v1.z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->tetrahedron_data[i].v2.x,
            this->tetrahedron_data[i].v2.y,
            this->tetrahedron_data[i].v2.z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->tetrahedron_data[i].v3.x,
            this->tetrahedron_data[i].v3.y,
            this->tetrahedron_data[i].v3.z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->tetrahedron_data[i].v4.x,
            this->tetrahedron_data[i].v4.y,
            this->tetrahedron_data[i].v4.z);
        file.write("\n");
    }
    file.close();
}

void shape_t::load_mesh(){
    shape_t::unset();
    file_t file;
    //
    file.open("mesh/mesh/mesh_data.txt", 'r');
    file.read("%zu", &this->N_points);
    file.read("%zu", &this->N_edges);
    file.read("%zu", &this->N_triangles);
    file.read("%zu", &this->N_tetrahedrons);
    file.close();
    shape_t::set();
    print("loading mesh information...");
    vector_t<real_t> v1, v2, v3, v4;
    int_t pg1, pg2, pg3, pg4;
    triangle_t triagle;
    tetrahedron_t tetrahedron;
    file.open("mesh/mesh/mesh_2d.txt", 'r');
    for (size_t i=0; i<this->N_triangles; i++){
        file.read("%lf %lf %lf %d\n", &v1.x, &v1.y, &v1.z, &pg1);
        file.read("%lf %lf %lf %d\n", &v2.x, &v2.y, &v2.z, &pg2);
        file.read("%lf %lf %lf %d\n", &v3.x, &v3.y, &v3.z, &pg3);
        file.read("\n");
        triagle.v1 = v1; triagle.v2 = v2; triagle.v3 = v3; triagle.get_area();
        assert_error((pg1==pg2)&&(pg2==pg3), "invalid physical groups");
        triagle.physical_group = pg1;
        this->triangle_data[i] = triagle;
    }
    file.close();
    file.open("mesh/mesh/mesh_3d.txt", 'r');
    for (size_t i=0; i<this->N_tetrahedrons; i++){
        file.read("%lf %lf %lf %d\n", &v1.x, &v1.y, &v1.z, &pg1);
        file.read("%lf %lf %lf %d\n", &v2.x, &v2.y, &v2.z, &pg2);
        file.read("%lf %lf %lf %d\n", &v3.x, &v3.y, &v3.z, &pg3);
        file.read("%lf %lf %lf %d\n", &v4.x, &v4.y, &v4.z, &pg4);
        file.read("\n");
        tetrahedron.v1 = v1; tetrahedron.v2 = v2; tetrahedron.v3 = v3; tetrahedron.v4 = v4;
        tetrahedron.get_volume();
        assert_error((pg1==pg2)&&(pg2==pg3)&&(pg3==pg4), "invalid physical groups");
        triagle.physical_group = pg1;
        this->tetrahedron_data[i] = tetrahedron;
    }
    file.close();
    print(", done!\n");
}

triangle_t shape_t::get_triangle_element(const size_t index){
    assert_error(index<this->N_triangles, "triangle index is out of range");
    return this->triangle_data[index];
}

tetrahedron_t shape_t::get_tetrahedron_element(const size_t index){
    assert_error(index<this->N_tetrahedrons, "tetrahedron index is out of range");
    return this->tetrahedron_data[index];
}
