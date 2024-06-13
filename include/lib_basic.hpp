#ifndef __LIB_BASIC_HPP__
#define __LIB_BASIC_HPP__

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <stdint.h>

// Definitions
#define true 1
#define false 0
#define null NULL

#define int_t int32_t
#define real_t double 
#ifndef real_t
#define real_t float
#endif

#define complex_t std::complex<real_t>
#define pi 3.141592653589793E0
#include "physics_constants.hpp"

// Functions

#endif