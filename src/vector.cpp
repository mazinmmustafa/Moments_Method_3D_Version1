//
#include "vector.hpp"

void print(vector_t<complex_t> A){
    printf("(%21.14E, %21.14E),\n(%21.14E, %21.14E),\n(%21.14E, %21.14E)\n",
        real(A.x), imag(A.x), 
        real(A.y), imag(A.y), 
        real(A.z), imag(A.z));
}

void print(vector_t<real_t> A){
    printf("(%21.14E, %21.14E, %21.14E)\n", A.x, A.y, A.z);
}

vector_t<complex_t> operator + (const vector_t<complex_t> A, const vector_t<complex_t> B){
    return vector_t<complex_t>(A.x+B.x, A.y+B.y, A.z+B.z);
}

vector_t<real_t> operator + (const vector_t<real_t> A, const vector_t<real_t> B){
    return vector_t<real_t>(A.x+B.x, A.y+B.y, A.z+B.z);
}

vector_t<complex_t> operator - (const vector_t<complex_t> A, const vector_t<complex_t> B){
    return vector_t<complex_t>(A.x-B.x, A.y-B.y, A.z-B.z);
}

vector_t<real_t> operator - (const vector_t<real_t> A, const vector_t<real_t> B){
    return vector_t<real_t>(A.x-B.x, A.y-B.y, A.z-B.z);
}

complex_t operator * (const vector_t<complex_t> A, const vector_t<complex_t> B){
    return A.x*B.x+A.y*B.y+A.z*B.z;
}

real_t operator * (const vector_t<real_t> A, const vector_t<real_t> B){
    return A.x*B.x+A.y*B.y+A.z*B.z;
}

vector_t<complex_t> operator ^ (const vector_t<complex_t> A, const vector_t<complex_t> B){
    return vector_t<complex_t>(A.y*B.z-A.z*B.y,
                               A.z*B.x-A.x*B.z,
                               A.x*B.y-A.y*B.x);
}

vector_t<real_t> operator ^ (const vector_t<real_t> A, const vector_t<real_t> B){
    return vector_t<real_t>(A.y*B.z-A.z*B.y,
                            A.z*B.x-A.x*B.z,
                            A.x*B.y-A.y*B.x);
}

vector_t<complex_t> operator * (const complex_t a, const vector_t<complex_t> A){
    return vector_t<complex_t>(a*A.x, a*A.y, a*A.z);
}

vector_t<complex_t> operator * (const real_t a, const vector_t<complex_t> A){
    return vector_t<complex_t>(a*A.x, a*A.y, a*A.z);
}

vector_t<real_t> operator * (const real_t a, const vector_t<real_t> A){
    return vector_t<real_t>(a*A.x, a*A.y, a*A.z);
}

vector_t<complex_t> operator * (const vector_t<complex_t> A, const complex_t a){
    return vector_t<complex_t>(a*A.x, a*A.y, a*A.z);
}

vector_t<complex_t> operator * (const vector_t<complex_t> A, const real_t a){
    return vector_t<complex_t>(a*A.x, a*A.y, a*A.z);
}

vector_t<real_t> operator * (const vector_t<real_t> A, const real_t a){
    return vector_t<real_t>(a*A.x, a*A.y, a*A.z);
}

vector_t<complex_t> operator / (const vector_t<complex_t> A, const complex_t a){
    return vector_t<complex_t>(A.x/a, A.y/a, A.z/a);
}

vector_t<complex_t> operator / (const vector_t<complex_t> A, const real_t a){
    return vector_t<complex_t>(A.x/a, A.y/a, A.z/a);
}

vector_t<real_t> operator / (const vector_t<real_t> A, const real_t a){
    return vector_t<real_t>(A.x/a, A.y/a, A.z/a);
}

real_t mag(const vector_t<real_t> A){
    return sqrt(A*A);
}

vector_t<real_t> unit(const vector_t<real_t> A){
    real_t A_mag=mag(A);
    return vector_t<real_t>(A.x/A_mag, A.y/A_mag, A.z/A_mag);
}