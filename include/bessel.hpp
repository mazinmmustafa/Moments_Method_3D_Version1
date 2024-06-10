#ifndef __BESSEL_HPP__
#define __BESSEL_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"

// Definitions

// Functions
complex_t besselj(const real_t n, const complex_t z);
complex_t bessely(const real_t n, const complex_t z);
complex_t besseli(const real_t n, const complex_t z);
complex_t besselk(const real_t n, const complex_t z);
complex_t besselh(const int_t k, const real_t n, const complex_t z);

#endif