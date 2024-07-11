#ifndef __TESTBENCH_HPP__
#define __TESTBENCH_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
//
#include "shape.hpp"
#include "projection.hpp"
#include "engine_2d.hpp"

// Definitions

// Functions
void test_utilities();
void test_vector();
void test_read_write_binary_files();
void test_matrix();
void test_quadl();
void test_shape();
void test_engine_2d();
void test_Z_mn_2d();
void test_RCS_sphere_2d();
void test_RCS_shape_2d();
void test_near_field_2d();
void test_near_field_heat_map_2d();
void test_current_2d();
void test_far_field_2d();
void test_far_field_antenna_2d();

#endif