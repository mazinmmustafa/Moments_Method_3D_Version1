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

#ifdef _WIN64
#define __windows__
#endif

#define int_t int32_t
#define real_t double
#define complex_t std::complex<real_t>
#include "physics_constants.hpp"

// Functions

#endif