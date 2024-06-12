//
#include "lib_basic.hpp"
#include "testbench.hpp"

int main(){

    test_utilities();
    test_vector();
    test_read_write_binary_files();
    test_matrix();
    test_quadl();

    return 0;

    // gmsh mesh/shape.geo -2 -format stl -save_all -o mesh/shape.stl

}