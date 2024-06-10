//
#include "matrix.hpp"

int check_dimensions_equal(const size_t rows_A, const size_t cols_A,
    const size_t rows_B, const size_t cols_B){
    return (rows_A==rows_B)&&(cols_A==cols_B);
}
