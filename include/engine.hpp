#ifndef __ENGINE_HPP__
#define __ENGINE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "quadl.hpp"
//
#include "shape.hpp"

// Definitions
class engine_t{
    private:
    public:
        engine_t();
        ~engine_t();
};

// Functions
complex_t sinc(const complex_t x);
real_t sinc(const real_t x);

#endif