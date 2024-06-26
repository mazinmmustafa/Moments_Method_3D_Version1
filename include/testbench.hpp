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
#include "engine.hpp"

// Definitions

// Functions
void test_utilities();
void test_vector();
void test_read_write_binary_files();
void test_matrix();
void test_quadl();
void test_shape();
void test_engine();

#endif