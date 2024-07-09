#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"
#include "file.hpp"

// Definitions
template <typename type_t>
class matrix_t{
    private:
        size_t rows=0, cols=0;
        type_t **data=null;
        int is_allocated=false;
        int is_P_allocated=false;
        size_t *P=null;
    public:
        matrix_t(){}
        ~matrix_t(){
        }
    void set(const size_t rows, const size_t cols){
        assert_error(this->is_allocated==false, "matrix is already set");
        assert_error(this->is_P_allocated==false, "matrix is already set");
        assert_error(rows>0&&cols>0, "invalid dimensions");
        this->rows = rows;
        this->cols = cols;
        this->data = (type_t**)calloc(this->rows, sizeof(type_t*));
        assert(this->data!=null);
        for (size_t i=0; i<this->rows; i++){
            this->data[i] = (type_t*)calloc(this->cols, sizeof(type_t));
            assert(this->data[i]!=null);
        }
        this->is_allocated = true;
        matrix_t::zeros();
    }
    void unset(){
        if (this->is_allocated){
            for (size_t i=0; i<this->rows; i++){
                free(this->data[i]);
            }
            free(this->data);
            this->is_allocated = false;
        }
        if (this->is_P_allocated){
            free(this->P);
            this->is_P_allocated = false;
        }
    }
    type_t operator() (const size_t i, const size_t j)const{
        assert_error(this->is_allocated, "matrix is not set yet");
        assert_error(i<this->rows&&j<this->cols, "invalid index range");
        return this->data[i][j];
    }
    type_t& operator() (const size_t i, const size_t j){
        assert_error(this->is_allocated, "matrix is not set yet");
        assert_error(i<this->rows&&j<this->cols, "invalid index range");
        return this->data[i][j];
    }
    void copy(matrix_t *matrix){
        assert_error(matrix->is_set(), "matrix is not set yet");
        size_t rows, cols;
        (*matrix).get_dimensions(&rows, &cols);
        this->set(rows, cols);
        for (size_t i=0; i<this->rows; i++){
            for (size_t j=0; j<this->cols; j++){
                this->data[i][j] = (*matrix)(i, j);
            }
        }
    }
    void zeros(){
        assert_error(this->is_allocated, "matrix is not set yet");
        for (size_t i=0; i<this->rows; i++){
            for (size_t j=0; j<this->cols; j++){
                this->data[i][j] = (type_t) 0.0;
            }
        }
    }
    void eye(){
        assert_error(this->is_allocated, "matrix is not set yet");
        for (size_t i=0; i<this->rows; i++){
            for (size_t j=0; j<this->cols; j++){
                if (i==j){
                    this->data[i][j] = (type_t) 1.0;
                }else{
                    this->data[i][j] = (type_t) 0.0;
                }
            }
        }
    }
    void ones(){
        assert_error(this->is_allocated, "matrix is not set yet");
        for (size_t i=0; i<this->rows; i++){
            for (size_t j=0; j<this->cols; j++){
                this->data[i][j] = (type_t) 1.0;
            }
        }
    }
    void get_dimensions(size_t *rows, size_t *cols){
        assert_error(this->is_allocated, "matrix is not set yet");
        (*rows) = this->rows;
        (*cols) = this->cols;
    }
    int is_set(){
        return this->is_allocated;
    }
    int is_P_set(){
        return this->is_P_allocated;
    }
    void save(const char *filename){
        assert_error(this->is_allocated, "matrix is not set yet");
        binary_file_t file;
        file.open(filename, 'w');
        file.write(&(this->rows));
        file.write(&(this->cols));
        for (size_t i=0; i<this->rows; i++){
            for (size_t j=0; j<this->cols; j++){
                file.write(&(this->data[i][j]));
            }
        }
        file.close();
    }
    void load(const char *filename){
        matrix_t::unset();
        binary_file_t file;
        file.open(filename, 'r');
        file.read(&(this->rows));
        file.read(&(this->cols));
        matrix_t::set(rows, cols);
        for (size_t i=0; i<this->rows; i++){
            for (size_t j=0; j<this->cols; j++){
                file.read(&(this->data[i][j]));
            }
        }
        file.close();
    }
    void scale(const complex_t scalar){
        for (size_t i=0; i<this->rows; i++){
            for (size_t j=0; j<this->cols; j++){
                this->data[i][j]*=scalar;
            }
        }
    }
    void lup(){
        assert_error(this->is_allocated, "matrix is not allocated");
        assert_error(this->rows==this->cols, "matrix is not square");
        size_t N=this->cols;
        type_t *ptr;
        size_t i_max, i_temp; 
        real_t max_element, abs_element;
        this->P=(size_t*)calloc((N+1), sizeof(size_t));
        assert(this->P!=null);
        for (size_t i=0; i<=N; i++){
            this->P[i]=i;
        }
        for (size_t i=0; i<N; i++){
            max_element = 0.0;
            i_max = i;
            for (size_t k=i; k<N; k++){
                if ((abs_element= abs(this->data[k][i]))>max_element){ 
                    max_element = abs_element;
                    i_max = k;
                }
            }
            if (max_element==0.0){
                assert_error(false, "matrix is singular");
            }
            if (i_max!=i) {
                i_temp = this->P[i];
                this->P[i] = this->P[i_max];
                this->P[i_max] = i_temp;
                ptr = this->data[i];
                this->data[i] = this->data[i_max];
                this->data[i_max] = ptr;
                P[N]++;
            }
            for (size_t j=i+1; j<N; j++){
                this->data[j][i]/=this->data[i][i];
                for (size_t k=i+1; k<N; k++){
                    this->data[j][k]-=this->data[j][i]*this->data[i][k];
                }    
            }
        }
        this->is_P_allocated = true;
    }
    size_t P_data(const size_t i){
        return this->P[i];
    }
    complex_t det(){
        assert_error(this->is_P_set(), "matrix A is not lup decomposed");
        const size_t N=this->rows;
        complex_t sum=this->data[0][0];
        for (size_t i=1; i<N; i++){
            sum*=this->data[i][i];
        }
        return (this->P[N]-N)%2 == 0 ? sum : -sum;
    }
    void solve(matrix_t &b, matrix_t &x){
        assert_error(b.is_set(), "matrix b is not set yet");
        assert_error(this->is_P_set(), "matrix A is not lup decomposed");
        size_t rows_b, cols_b;
        b.get_dimensions(&rows_b, &cols_b);
        assert_error(this->rows==rows_b, "invalid matrix dimensions");
        assert_error(cols_b==1, "invalid matrix dimensions");
        const size_t N=this->rows;
        x.set(N, 1);
        for (size_t i=0; i<N; i++){
            x(i, 0) = b(this->P[i], 0);
            for (size_t k=0; k<i; k++){
                x(i, 0)-=this->data[i][k]*x(k, 0);
            }
        }
        for (int_t i=N-1; i>=0; i--){
            for (size_t k=(size_t)i+1; k<N; k++){
                x((size_t)i, 0)-=this->data[(size_t)i][k]*x(k, 0);
            }
            x((size_t)i, 0)/=data[(size_t)i][i];
        }
    }
    void inv(){
        assert_error(this->is_P_set(), "matrix A is not lup decomposed");
        size_t N=this->rows;
        matrix_t LUP;
        LUP.copy(this);
        for (size_t j=0; j<N; j++) {
            for (size_t i=0; i<N; i++) {
                this->data[i][j] = this->P[i] == j ? 1.0 : 0.0;
                for (size_t k=0; k<i; k++){
                    this->data[i][j]-=LUP(i, k)*this->data[k][j];
                }
            }
            for (int_t i=N-1; i>=0; i--) {
                for (size_t k=(size_t)i+1; k<N; k++){
                    this->data[(size_t)i][j]-=LUP((size_t)i, k)*this->data[k][j];
                }
                this->data[(size_t)i][j]/=LUP((size_t)i, (size_t)i);
            }
        }
        LUP.unset();
    }
    void disp(){
        const size_t dim_max=4;
        assert_error(this->is_set(), "matrix is not set yet");
        for (size_t i=0; i<(this->rows<dim_max?this->rows:dim_max); i++){
            for (size_t j=0; j<(this->cols<dim_max?this->cols:dim_max); j++){
                print("element %zu, %zu:\n", i, j);
                print(this->data[i][j]);
            }
        }
    }
};

// Functions

int check_dimensions_equal(const size_t rows_A, const size_t cols_A,
    const size_t rows_B, const size_t cols_B);

template <typename type_t>
void add_matrix(matrix_t<type_t> &A, matrix_t<type_t> &B, matrix_t<type_t> &C){
    size_t rows_A, cols_A;
    size_t rows_B, cols_B;
    A.get_dimensions(&rows_A, &cols_A);
    B.get_dimensions(&rows_B, &cols_B);
    assert_error(check_dimensions_equal(rows_A, cols_A, rows_B, cols_B), 
        "invalid matrices dimensions");
    C.set(rows_A, cols_A);
    for (size_t i=0; i<rows_A; i++){
        for (size_t j=0; j<cols_A; j++){
            C(i, j) = A(i, j)+B(i, j);
        }
    }
}

template <typename type_t>
void sub_matrix(matrix_t<type_t> &A, matrix_t<type_t> &B, matrix_t<type_t> &C){
    size_t rows_A, cols_A;
    size_t rows_B, cols_B;
    A.get_dimensions(&rows_A, &cols_A);
    B.get_dimensions(&rows_B, &cols_B);
    assert_error(check_dimensions_equal(rows_A, cols_A, rows_B, cols_B), 
        "invalid matrices dimensions");
    C.set(rows_A, cols_A);
    for (size_t i=0; i<rows_A; i++){
        for (size_t j=0; j<cols_A; j++){
            C(i, j) = A(i, j)-B(i, j);
        }
    }
}

template <typename type_t>
void mult_matrix(matrix_t<type_t> &A, matrix_t<type_t> &B, matrix_t<type_t> &C){
    size_t rows_A, cols_A;
    size_t rows_B, cols_B;
    A.get_dimensions(&rows_A, &cols_A);
    B.get_dimensions(&rows_B, &cols_B);
    assert_error(cols_A==rows_B, "invalid matrices dimensions");
    C.set(rows_A, cols_B);
    complex_t temp;
    for (size_t i=0; i<rows_A; i++){
        for (size_t j=0; j<cols_B; j++){
            temp = 0.0;
            for (size_t k=0; k<cols_A; k++){
                temp+=A(i, k)*B(k, j);
            }
            C(i, j) = temp;
        }
    }
}

#endif