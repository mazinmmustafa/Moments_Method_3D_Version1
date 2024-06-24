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
    if (this->is_basis_2d_list_allocated){
        free(this->basis_2d_list);
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
    triangle_t triangle;
    tetrahedron_t tetrahedron;
    file.open("mesh/mesh/mesh_2d.txt", 'r');
    for (size_t i=0; i<this->N_triangles; i++){
        file.read("%lf %lf %lf %d\n", &v1.x, &v1.y, &v1.z, &pg1);
        file.read("%lf %lf %lf %d\n", &v2.x, &v2.y, &v2.z, &pg2);
        file.read("%lf %lf %lf %d\n", &v3.x, &v3.y, &v3.z, &pg3);
        file.read("\n");
        triangle.v1 = v1; triangle.v2 = v2; triangle.v3 = v3; triangle.get_area();
        assert_error(triangle.area>0.0, "invalid triangle element");
        assert_error((pg1==pg2)&&(pg2==pg3), "invalid physical groups");
        triangle.physical_group = pg1;
        this->triangle_data[i] = triangle;
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
        assert_error(tetrahedron.volume>0.0, "invalid tetrahedron element");
        assert_error((pg1==pg2)&&(pg2==pg3)&&(pg3==pg4), "invalid physical groups");
        tetrahedron.physical_group = pg1;
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

void shape_t::get_mesh(){
    int_t return_value;
    #ifdef __windows__
        return_value = system("python mesh/generate_mesh.py");
    #endif
    #ifdef __linux__
        return_value = system("python3 mesh/generate_mesh.py");
    #endif
    assert_error(return_value==0, "failed to call python");
}

void shape_t::assign_volume_properties(const complex_t eps, const int_t physical_group){
    for (size_t i=0; i<this->N_tetrahedrons; i++){
        if (this->tetrahedron_data[i].physical_group==physical_group){
            this->tetrahedron_data[i].eps = eps;
        }
    }
}

void print_log(const triangle_t tri_s, const triangle_t tri_d){
    print("I found basis:\n");
    print(tri_s.v1);
    print(tri_s.v2);
    print(tri_s.v3);
    print(tri_d.v1);
    print(tri_d.v2);
    print(tri_d.v3);
    print("\n");
}

void shape_t::get_basis_functions(){
    assert_error(this->is_triangle_allocated, "no 2d shape elements were found");
    assert_error(!this->is_basis_2d_list_allocated, "2d basis functions were already allocated");
    file_t file;
    file.open("mesh/basis/basis_2d.txt", 'w');
    vector_t<real_t> v1_s, v2_s, v3_s, v1_d, v2_d, v3_d;
    triangle_t triangle_s, triangle_d;
    for (size_t i=0; i<this->N_triangles; i++){
        v1_s = this->triangle_data[i].v1;
        v2_s = this->triangle_data[i].v2;
        v3_s = this->triangle_data[i].v3;
        triangle_s = this->triangle_data[i];
        for (size_t j=(i+1); j<this->N_triangles; j++){
            v1_d = this->triangle_data[j].v1;
            v2_d = this->triangle_data[j].v2;
            v3_d = this->triangle_data[j].v3;
            triangle_d = this->triangle_data[j];
            // group 1
            if (is_equal(v2_s, v2_d, tol_vertex)&&is_equal(v3_s, v1_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v3.x, triangle_d.v3.y, triangle_d.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            //
            if (is_equal(v2_s, v1_d, tol_vertex)&&is_equal(v3_s, v3_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v2.x, triangle_d.v2.y, triangle_d.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            // 
            if (is_equal(v2_s, v3_d, tol_vertex)&&is_equal(v3_s, v2_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v1.x, triangle_d.v1.y, triangle_d.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            // group 2
            if (is_equal(v3_s, v3_d, tol_vertex)&&is_equal(v1_s, v2_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v1.x, triangle_d.v1.y, triangle_d.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            //
            if (is_equal(v3_s, v2_d, tol_vertex)&&is_equal(v1_s, v1_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v3.x, triangle_d.v3.y, triangle_d.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            //
            if (is_equal(v3_s, v1_d, tol_vertex)&&is_equal(v1_s, v3_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v2.x, triangle_d.v2.y, triangle_d.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            // group 3
            if (is_equal(v1_s, v1_d, tol_vertex)&&is_equal(v2_s, v3_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v2.x, triangle_d.v2.y, triangle_d.v2.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            //
            if (is_equal(v1_s, v3_d, tol_vertex)&&is_equal(v2_s, v2_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v1.x, triangle_d.v1.y, triangle_d.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
            //
            if (is_equal(v1_s, v2_d, tol_vertex)&&is_equal(v2_s, v1_d, tol_vertex)){
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v3.x, triangle_s.v3.y, triangle_s.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v1.x, triangle_s.v1.y, triangle_s.v1.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_d.v3.x, triangle_d.v3.y, triangle_d.v3.z);
                file.write("%21.14E %21.14E %21.14E ", triangle_s.v2.x, triangle_s.v2.y, triangle_s.v2.z);
                file.write("%d %zu\n", triangle_s.physical_group, triangle_d.physical_group);
                triangle_s.N_adjacents++;
                triangle_d.N_adjacents++;
                this->N_2d_basis++;
            }
        }
        print_log(triangle_s, triangle_d);
        char message[100];
        sprintf(message, "invalid triangle %zu found", i);
        assert_error(triangle_s.N_adjacents>0&&triangle_s.N_adjacents<4, message);
    }

    file.close();
}