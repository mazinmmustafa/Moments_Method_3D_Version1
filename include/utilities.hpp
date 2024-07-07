#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

// Libraries
#include "lib_basic.hpp"

// Definitions
#define assert_error(check, msg) __assert_error((check), (msg), __FILE__, __LINE__);

class range_t{
    private:
        real_t *data=null;
        real_t x_min=0.0, x_max=0.0;
        size_t Ns=0;
        int is_allocated=false;
    public:
        range_t();
        ~range_t();
        void linspace();
        void logspace();
        real_t operator() (const size_t index) const;
        void set(const real_t x_min, const real_t x_max, const size_t Ns);
        void unset();
        void get_info(real_t *x_min, real_t *x_max, size_t *Ns);
};

class timer_lib_t{
private:
    int is_set=false;
    #ifdef __windows__
    time_t start, stop;
    #endif
    #ifdef __linux__
    struct timespec start, stop;
    #endif
    double elapsed=0.0;
public:
    timer_lib_t();
    ~timer_lib_t();
    void set();
    void unset();
    void unset_silent();
    double get_elapsed();
};

// Functions
void print(const char *format, ...);
void print(const int_t n);
void print(const size_t n);
void print(const real_t x);
void print(const complex_t z);
void __assert_error(const int_t condition, const char *error_msg, 
    const char* filename, const size_t line);
void progress_bar(const size_t i, const size_t N, const char *msg);
complex_t sinc(const complex_t x);
real_t sinc(const real_t x);
complex_t sinc_dx(const complex_t x);
real_t sinc_dx(const real_t x);
real_t deg2rad(const real_t theta);
real_t rad2deg(const real_t theta);

#endif