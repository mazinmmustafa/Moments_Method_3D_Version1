//
#include "lib_basic.hpp"
#include "testbench.hpp"

int main(){

    // test_utilities();
    // test_vector();
    // test_read_write_binary_files();
    // test_matrix();
    // test_quadl();
    test_shape();

    return 0;

}

// gmsh mesh/shape.brep -3 -clmax 0.05 -format vtk -save_all -o mesh/shape.vtk