#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"

// Definitions
class matrix_t{
    private:
        size_t rows=0, cols=0;
        complex_t **data=NULL;
        int is_allocated=false;
        int is_P_allocated=false;
        size_t *P=NULL;
    public:
        matrix_t();
        ~matrix_t();
        void set(const size_t rows, const size_t cols);
        void unset();
        complex_t operator() (const size_t i, const size_t j)const;
        complex_t& operator() (const size_t i, const size_t j);
        size_t P_data(const size_t i);
        void copy(matrix_t *matrix); // copy from matrix to this instance
        void zeros();
        void eye();
        void ones();
        void get_dimensions(size_t *rows, size_t *cols);
        int is_set();
        int is_P_set();
        void save(const char *filename);
        void load(const char *filename);
        void scale(const complex_t scalar);
        void lup();
        complex_t det();
        void solve(matrix_t &b, matrix_t &x);
        void inv();
};

// Functions
void print(matrix_t &matrix);
void add_matrix(matrix_t &A, matrix_t &B, matrix_t &C);
void sub_matrix(matrix_t &A, matrix_t &B, matrix_t &C);
void mult_matrix(matrix_t &A, matrix_t &B, matrix_t &C);

#endif