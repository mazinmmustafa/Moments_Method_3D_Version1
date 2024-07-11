//
#include "shape.hpp"

const real_t tol_vertex=1.0E-4;

shape_t::shape_t(){
}

shape_t::~shape_t(){
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
    if (this->is_basis_1d_list_allocated){
        free(this->basis_1d_list);
    }
    if (this->is_basis_2d_list_allocated){
        free(this->basis_2d_list);
    }
    if (this->is_basis_3d_list_allocated){
        free(this->basis_3d_list);
    }
    this->N_points = 0;
    this->N_edges = 0;
    this->N_triangles = 0;
    this->N_tetrahedrons = 0;
    this->is_edge_allocated = false;
    this->is_triangle_allocated = false;
    this->is_tetrahedron_allocated = false;
}

void shape_t::set(){
    this->edge_data = (edge_t*)calloc(this->N_edges, sizeof(edge_t));
    this->triangle_data = (triangle_t*)calloc(this->N_triangles, sizeof(triangle_t));
    this->tetrahedron_data = (tetrahedron_t*)calloc(this->N_tetrahedrons, sizeof(tetrahedron_t));
    assert(this->edge_data!=null);
    assert(this->triangle_data!=null);
    assert(this->tetrahedron_data!=null);
    this->is_edge_allocated = true;
    this->is_triangle_allocated = true;
    this->is_tetrahedron_allocated = true;
}

void shape_t::log_mesh(){
    file_t file;
    //
    file.open("mesh/mesh/mesh_1d_elements_log.txt", 'w');
    file.write("number of edge: %zu\n\n", this->N_edges);
    for (size_t i=0; i<this->N_edges; i++){
        file.write("edge: %zu: group: %d\n", i, 
            this->edge_data[i].physical_group);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->edge_data[i].v[0].x,
            this->edge_data[i].v[0].y,
            this->edge_data[i].v[0].z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->edge_data[i].v[1].x,
            this->edge_data[i].v[1].y,
            this->edge_data[i].v[1].z);
        file.write("\n");
    }
    file.close();
    //
    file.open("mesh/mesh/mesh_2d_elements_log.txt", 'w');
    file.write("number of triangles: %zu\n\n", this->N_triangles);
    for (size_t i=0; i<this->N_triangles; i++){
        file.write("triangle: %zu: group: %d\n", i, 
            this->triangle_data[i].physical_group);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->triangle_data[i].v[0].x,
            this->triangle_data[i].v[0].y,
            this->triangle_data[i].v[0].z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->triangle_data[i].v[1].x,
            this->triangle_data[i].v[1].y,
            this->triangle_data[i].v[1].z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->triangle_data[i].v[2].x,
            this->triangle_data[i].v[2].y,
            this->triangle_data[i].v[2].z);
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
            this->tetrahedron_data[i].v[0].x,
            this->tetrahedron_data[i].v[0].y,
            this->tetrahedron_data[i].v[0].z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->tetrahedron_data[i].v[1].x,
            this->tetrahedron_data[i].v[1].y,
            this->tetrahedron_data[i].v[1].z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->tetrahedron_data[i].v[2].x,
            this->tetrahedron_data[i].v[2].y,
            this->tetrahedron_data[i].v[2].z);
        file.write("%21.14E, %21.14E, %21.14E\n", 
            this->tetrahedron_data[i].v[3].x,
            this->tetrahedron_data[i].v[3].y,
            this->tetrahedron_data[i].v[3].z);
        file.write("\n");
    }
    file.close();
}

void shape_t::load_mesh(){
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
    edge_t edge;
    triangle_t triangle;
    tetrahedron_t tetrahedron;
    file.open("mesh/mesh/mesh_1d.txt", 'r');
    for (size_t i=0; i<this->N_edges; i++){
        file.read("%lf %lf %lf %d\n", &v1.x, &v1.y, &v1.z, &pg1);
        file.read("%lf %lf %lf %d\n", &v2.x, &v2.y, &v2.z, &pg2);
        file.read("\n");
        edge.v[0] = v1; edge.v[1] = v2; edge.get_length();
        assert_error(edge.length>0.0, "invalid edge element");
        assert_error((pg1==pg2), "invalid physical groups");
        edge.physical_group = pg1;
        this->edge_data[i] = edge;
    }
    file.close();
    file.open("mesh/mesh/mesh_2d.txt", 'r');
    for (size_t i=0; i<this->N_triangles; i++){
        file.read("%lf %lf %lf %d\n", &v1.x, &v1.y, &v1.z, &pg1);
        file.read("%lf %lf %lf %d\n", &v2.x, &v2.y, &v2.z, &pg2);
        file.read("%lf %lf %lf %d\n", &v3.x, &v3.y, &v3.z, &pg3);
        file.read("\n");
        triangle.v[0] = v1; triangle.v[1] = v2; triangle.v[2] = v3; triangle.get_area();
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
        tetrahedron.v[0] = v1; tetrahedron.v[1] = v2; tetrahedron.v[2] = v3; tetrahedron.v[3] = v4;
        tetrahedron.get_volume();
        assert_error(tetrahedron.volume>0.0, "invalid tetrahedron element");
        assert_error((pg1==pg2)&&(pg2==pg3)&&(pg3==pg4), "invalid physical groups");
        tetrahedron.physical_group = pg1;
        this->tetrahedron_data[i] = tetrahedron;
    }
    file.close();
    print(", done!\n");
}

edge_t shape_t::get_edge_element(const size_t index){
    assert_error(index<this->N_edges, "edge index is out of range");
    return this->edge_data[index];
}

triangle_t shape_t::get_triangle_element(const size_t index){
    assert_error(index<this->N_triangles, "triangle index is out of range");
    return this->triangle_data[index];
}

tetrahedron_t shape_t::get_tetrahedron_element(const size_t index){
    assert_error(index<this->N_tetrahedrons, "tetrahedron index is out of range");
    return this->tetrahedron_data[index];
}

basis_1d_t shape_t::get_basis_1d(const size_t index){
    assert_error(index<this->N_1d_basis, "basis index is out of range");
    return this->basis_1d_list[index];
}

basis_2d_t shape_t::get_basis_2d(const size_t index){
    assert_error(index<this->N_2d_basis, "basis index is out of range");
    return this->basis_2d_list[index];
}

basis_3d_t shape_t::get_basis_3d(const size_t index){
    assert_error(index<this->N_3d_basis, "basis index is out of range");
    return this->basis_3d_list[index];
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
    shape_t::load_mesh();
    shape_t::log_mesh();
}

void shape_t::assign_volume_properties(const complex_t eps, const int_t physical_group){
    for (size_t i=0; i<this->N_tetrahedrons; i++){
        if (this->tetrahedron_data[i].physical_group==physical_group){
            this->tetrahedron_data[i].eps = eps;
        }
    }
}

size_t mod_1d(const size_t a){
    return a==1 ? 0 : 1;
}

size_t mod_2d(const size_t a){
    size_t ans=3-(size_t)(a%3);
    return ans==3 ? 0 : ans;
}

size_t mod_3d(const size_t a){
    size_t ans=6-(size_t)(a%6);
    return ans==6 ? 0 : ans;
}

void shape_t::get_basis_functions(const real_t unit_length){
    assert_error(this->is_triangle_allocated, "no 2d shape elements were found");
    assert_error(!this->is_basis_2d_list_allocated, "2d basis functions were already allocated");
    file_t file;
    // 1d basis
    file.open("mesh/basis/basis_1d.txt", 'w');
    edge_t edge_s, edge_d;
    size_t new_edge[1];
    size_t index_edge_s;
    size_t index_edge_d;
    for (size_t i=0; i<this->N_edges; i++){
        edge_s = this->edge_data[i];
        if (edge_s.physical_group>0){
            for (size_t j=(i+1); j<this->N_edges; j++){
                edge_d = this->edge_data[j];
                // 
                size_t count=0;
                index_edge_s = 0;
                index_edge_d = 0;
                for (size_t ii=0; ii<2; ii++){
                    for (size_t jj=0; jj<2; jj++){
                        if (is_equal(edge_s.v[ii], edge_d.v[jj], this->lambda*tol_vertex/unit_length)){
                            new_edge[count] = ii;
                            index_edge_s+=ii;
                            index_edge_d+=jj;
                            count++;
                            if (count==1){
                                break;
                            }
                        }
                    }
                }
                if (count==1&&edge_d.physical_group>0){
                    index_edge_s = mod_1d(index_edge_s);
                    index_edge_d = mod_1d(index_edge_d);
                    file.write("%21.14E %21.14E %21.14E ", 
                        edge_s.v[index_edge_s].x*unit_length,
                        edge_s.v[index_edge_s].y*unit_length,
                        edge_s.v[index_edge_s].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        edge_s.v[new_edge[0]].x*unit_length,
                        edge_s.v[new_edge[0]].y*unit_length,
                        edge_s.v[new_edge[0]].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        edge_d.v[index_edge_d].x*unit_length,
                        edge_d.v[index_edge_d].y*unit_length,
                        edge_d.v[index_edge_d].z*unit_length);
                    file.write("%d %d\n", edge_s.physical_group, edge_d.physical_group);
                    count = 0;
                    edge_s.N_adjacents++;
                    edge_d.N_adjacents++;
                    this->N_1d_basis++;
                }
            }
        }
    }
    file.close();
    // 2d basis
    file.open("mesh/basis/basis_2d.txt", 'w');
    triangle_t triangle_s, triangle_d;
    size_t new_triangle[2];
    size_t index_triangle_s;
    size_t index_triangle_d;
    for (size_t i=0; i<this->N_triangles; i++){
        triangle_s = this->triangle_data[i];
        if (triangle_s.physical_group>0){
            for (size_t j=(i+1); j<this->N_triangles; j++){
                triangle_d = this->triangle_data[j];
                // 
                size_t count=0;
                index_triangle_s = 0;
                index_triangle_d = 0;
                for (size_t ii=0; ii<3; ii++){
                    for (size_t jj=0; jj<3; jj++){
                        if (is_equal(triangle_s.v[ii], triangle_d.v[jj], this->lambda*tol_vertex/unit_length)){
                            new_triangle[count] = ii;
                            index_triangle_s+=ii;
                            index_triangle_d+=jj;
                            count++;
                            if (count==2){
                                break;
                            }
                        }
                    }
                }
                if (count==2&&triangle_d.physical_group>0){
                    index_triangle_s = mod_2d(index_triangle_s);
                    index_triangle_d = mod_2d(index_triangle_d);
                    file.write("%21.14E %21.14E %21.14E ", 
                        triangle_s.v[index_triangle_s].x*unit_length,
                        triangle_s.v[index_triangle_s].y*unit_length,
                        triangle_s.v[index_triangle_s].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        triangle_s.v[new_triangle[0]].x*unit_length,
                        triangle_s.v[new_triangle[0]].y*unit_length,
                        triangle_s.v[new_triangle[0]].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        triangle_s.v[new_triangle[1]].x*unit_length,
                        triangle_s.v[new_triangle[1]].y*unit_length,
                        triangle_s.v[new_triangle[1]].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        triangle_d.v[index_triangle_d].x*unit_length,
                        triangle_d.v[index_triangle_d].y*unit_length,
                        triangle_d.v[index_triangle_d].z*unit_length);
                    file.write("%d %d\n", triangle_s.physical_group, triangle_d.physical_group);
                    count = 0;
                    triangle_s.N_adjacents++;
                    triangle_d.N_adjacents++;
                    this->N_2d_basis++;
                }
            }
        }
    }
    file.close();
    // 3d basis
    file.open("mesh/basis/basis_3d.txt", 'w');
    tetrahedron_t tetrahedron_s, tetrahedron_d;
    size_t new_tetrahedron[3];
    size_t index_tetrahedron_s;
    size_t index_tetrahedron_d;
    for (size_t i=0; i<this->N_tetrahedrons; i++){
        tetrahedron_s = this->tetrahedron_data[i];
        if (tetrahedron_s.physical_group>0){
            for (size_t j=(i+1); j<this->N_tetrahedrons; j++){
                tetrahedron_d = this->tetrahedron_data[j];
                // 
                size_t count=0;
                index_tetrahedron_s = 0;
                index_tetrahedron_d = 0;
                for (size_t ii=0; ii<4; ii++){
                    for (size_t jj=0; jj<4; jj++){
                        if (is_equal(tetrahedron_s.v[ii], tetrahedron_d.v[jj], this->lambda*tol_vertex/unit_length)){
                            new_tetrahedron[count] = ii;
                            index_tetrahedron_s+=ii;
                            index_tetrahedron_d+=jj;
                            count++;
                            if (count==3){
                                break;
                            }
                        }
                    }
                }
                if (count==3&&tetrahedron_d.physical_group>0){
                    index_tetrahedron_s = mod_3d(index_tetrahedron_s);
                    index_tetrahedron_d = mod_3d(index_tetrahedron_d);
                    file.write("%21.14E %21.14E %21.14E ", 
                        tetrahedron_s.v[index_tetrahedron_s].x*unit_length,
                        tetrahedron_s.v[index_tetrahedron_s].y*unit_length,
                        tetrahedron_s.v[index_tetrahedron_s].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        tetrahedron_s.v[new_tetrahedron[0]].x*unit_length,
                        tetrahedron_s.v[new_tetrahedron[0]].y*unit_length,
                        tetrahedron_s.v[new_tetrahedron[0]].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        tetrahedron_s.v[new_tetrahedron[1]].x*unit_length,
                        tetrahedron_s.v[new_tetrahedron[1]].y*unit_length,
                        tetrahedron_s.v[new_tetrahedron[1]].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        tetrahedron_s.v[new_tetrahedron[2]].x*unit_length,
                        tetrahedron_s.v[new_tetrahedron[2]].y*unit_length,
                        tetrahedron_s.v[new_tetrahedron[2]].z*unit_length);
                    file.write("%21.14E %21.14E %21.14E ", 
                        tetrahedron_d.v[index_tetrahedron_d].x*unit_length,
                        tetrahedron_d.v[index_tetrahedron_d].y*unit_length,
                        tetrahedron_d.v[index_tetrahedron_d].z*unit_length);
                    file.write("%d %d\n", tetrahedron_s.physical_group, tetrahedron_d.physical_group);
                    count = 0;
                    tetrahedron_s.N_adjacents++;
                    tetrahedron_d.N_adjacents++;
                    this->N_3d_basis++;
                }
            }
        }
    }
    file.close();
    //
    print("number of 1d basis functions: ");
    print(this->N_1d_basis);
    print("number of 2d basis functions: ");
    print(this->N_2d_basis);
    print("number of 3d basis functions: ");
    print(this->N_3d_basis);
    print("total number of basis functions: ");
    print(this->N_1d_basis+this->N_2d_basis+this->N_3d_basis);
    assert_error((this->N_1d_basis+this->N_2d_basis+this->N_3d_basis)<__max_system__, 
        "too many basis functions");
    file.open("mesh/basis/basis_info.txt", 'w');
    file.write("%zu\n%zu\n%zu", this->N_1d_basis, this->N_2d_basis, this->N_3d_basis);
    file.close();
    shape_t::free_basic_elements();
    shape_t::load_basis_functions();
}

void shape_t::free_basic_elements(){
    if (this->is_edge_allocated){
        free(this->edge_data); 
        this->is_edge_allocated = false;
    }
    if (this->is_triangle_allocated){
        free(this->triangle_data); 
        this->is_triangle_allocated = false;
    }
    if (this->is_tetrahedron_allocated){
        free(this->tetrahedron_data); 
        this->is_tetrahedron_allocated = false;
    }
}

void shape_t::load_basis_functions(){   
    assert_error(this->is_medium_set, "medium parameters not set yet");
    file_t file;
    file.open("mesh/basis/basis_info.txt", 'r');
    file.read("%zu", &this->N_1d_basis);
    file.read("%zu", &this->N_2d_basis);
    file.read("%zu", &this->N_3d_basis);
    file.close();
    if (this->N_1d_basis>0&&this->is_basis_1d_list_allocated==false){
        this->basis_1d_list = (basis_1d_t*)calloc(this->N_1d_basis, sizeof(basis_1d_t));
        assert(this->basis_1d_list!=null);
        this->is_basis_1d_list_allocated = true;
    }
    if (this->N_2d_basis>0&&this->is_basis_2d_list_allocated==false){
        this->basis_2d_list = (basis_2d_t*)calloc(this->N_2d_basis, sizeof(basis_2d_t));
        assert(this->basis_2d_list!=null);
        this->is_basis_2d_list_allocated = true;
    }
    if (this->N_3d_basis>0&&this->is_basis_3d_list_allocated==false){
        this->basis_3d_list = (basis_3d_t*)calloc(this->N_3d_basis, sizeof(basis_3d_t));
        assert(this->basis_3d_list!=null);
        this->is_basis_3d_list_allocated = true;
    }
    basis_1d_t basis_1d;
    basis_2d_t basis_2d;
    basis_3d_t basis_3d;
    //
    file.open("mesh/basis/basis_1d.txt", 'r');
    for (size_t i=0; i<this->N_1d_basis; i++){
        file.read("%lf", &basis_1d.r_m.x);
        file.read("%lf", &basis_1d.r_m.y);
        file.read("%lf", &basis_1d.r_m.z);
        file.read("%lf", &basis_1d.e_1.x);
        file.read("%lf", &basis_1d.e_1.y);
        file.read("%lf", &basis_1d.e_1.z);
        file.read("%lf", &basis_1d.r_p.x);
        file.read("%lf", &basis_1d.r_p.y);
        file.read("%lf", &basis_1d.r_p.z);
        file.read("%d", &basis_1d.physical_group_m);
        file.read("%d", &basis_1d.physical_group_p);
        basis_1d.get_values();
        this->basis_1d_list[i] = basis_1d;
    }
    file.close();
    //
    file.open("mesh/basis/basis_2d.txt", 'r');
    for (size_t i=0; i<this->N_2d_basis; i++){
        basis_2d_t basis_2d;
        file.read("%lf", &basis_2d.r_m.x);
        file.read("%lf", &basis_2d.r_m.y);
        file.read("%lf", &basis_2d.r_m.z);
        file.read("%lf", &basis_2d.e_1.x);
        file.read("%lf", &basis_2d.e_1.y);
        file.read("%lf", &basis_2d.e_1.z);
        file.read("%lf", &basis_2d.e_2.x);
        file.read("%lf", &basis_2d.e_2.y);
        file.read("%lf", &basis_2d.e_2.z);
        file.read("%lf", &basis_2d.r_p.x);
        file.read("%lf", &basis_2d.r_p.y);
        file.read("%lf", &basis_2d.r_p.z);
        file.read("%d", &basis_2d.physical_group_m);
        file.read("%d", &basis_2d.physical_group_p);
        basis_2d.get_values();
        this->basis_2d_list[i] = basis_2d;
    }
    file.close();
    //
    file.open("mesh/basis/basis_3d.txt", 'r');
    for (size_t i=0; i<this->N_3d_basis; i++){
        basis_3d_t basis_3d;
        file.read("%lf", &basis_3d.r_m.x);
        file.read("%lf", &basis_3d.r_m.y);
        file.read("%lf", &basis_3d.r_m.z);
        file.read("%lf", &basis_3d.e_1.x);
        file.read("%lf", &basis_3d.e_1.y);
        file.read("%lf", &basis_3d.e_1.z);
        file.read("%lf", &basis_3d.e_2.x);
        file.read("%lf", &basis_3d.e_2.y);
        file.read("%lf", &basis_3d.e_2.z);
        file.read("%lf", &basis_3d.e_3.x);
        file.read("%lf", &basis_3d.e_3.y);
        file.read("%lf", &basis_3d.e_3.z);
        file.read("%lf", &basis_3d.r_p.x);
        file.read("%lf", &basis_3d.r_p.y);
        file.read("%lf", &basis_3d.r_p.z);
        file.read("%d", &basis_3d.physical_group_m);
        file.read("%d", &basis_3d.physical_group_p);
        basis_3d.get_values();
        this->basis_3d_list[i] = basis_3d;
    }
    file.close();
}

shape_info_t shape_t::get_shape_info(){
    shape_info_t shape_info={this->N_1d_basis, this->N_2d_basis, this->N_3d_basis, 
        this->is_basis_1d_list_allocated, this->is_basis_2d_list_allocated, this->is_basis_3d_list_allocated,
        this->k, this->eta, this->freq, this->lambda};
    return shape_info;
}

void shape_t::set_medium(const complex_t mu, const complex_t eps, const real_t freq){
    this->mu = mu;
    this->eps = eps;
    this->k = 2.0*pi*freq*sqrt(mu_0*eps_0)*sqrt(mu*eps);
    this->eta = sqrt(mu_0/eps_0)*sqrt(mu/eps);
    assert_error(freq>0.0, "invalid frequency");
    this->freq = freq;
    this->lambda = 2.0*pi/real(this->k);
    this->is_medium_set = true;
}

void shape_t::mesh_1d(const char *gmsh_filename, const real_t clmax){
    assert_error(this->is_medium_set, "medium is not set yet");
    int_t return_value;
    char *cmd=(char*)calloc(200, sizeof(char));
    sprintf(cmd, "gmsh mesh/%s -1 -clmax %0.4f -format vtk -save_all -o mesh/shape.vtk", 
        gmsh_filename, clmax);
    return_value = system(cmd);
    free(cmd);
    assert_error(return_value==0, "failed to call gmsh");
}

void shape_t::mesh_2d(const char *gmsh_filename, const real_t clmax){
    assert_error(this->is_medium_set, "medium is not set yet");
    int_t return_value;
    char *cmd=(char*)calloc(200, sizeof(char));
    sprintf(cmd, "gmsh mesh/%s -2 -clmax %0.4f -format vtk -save_all -o mesh/shape.vtk", 
        gmsh_filename, clmax);
    return_value = system(cmd);
    free(cmd);
    assert_error(return_value==0, "failed to call gmsh");
}

void shape_t::mesh_3d(const char *gmsh_filename, const real_t clmax){
    assert_error(this->is_medium_set, "medium is not set yet");
    int_t return_value;
    char *cmd=(char*)calloc(200, sizeof(char));
    sprintf(cmd, "gmsh mesh/%s -3 -clmax %0.4f -format vtk -save_all -o mesh/shape.vtk", 
        gmsh_filename, clmax);
    return_value = system(cmd);
    free(cmd);
    assert_error(return_value==0, "failed to call gmsh");
}

void shape_t::check(){
    // assert_error(this->is_basis_1d_list_allocated, "no 1d basis functions found");
    assert_error(this->is_basis_2d_list_allocated, "no 2d basis functions found");
    // assert_error(this->is_basis_3d_list_allocated, "no 3d basis functions found");
    assert_error(this->is_medium_set, "medium paramters not set yet");
}