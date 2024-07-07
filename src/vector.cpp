//
#include "vector.hpp"

void print(const vector_t<complex_t> A){
    printf("(%21.14E, %21.14E),\n(%21.14E, %21.14E),\n(%21.14E, %21.14E)\n",
        (real_t)real(A.x), (real_t)imag(A.x), 
        (real_t)real(A.y), (real_t)imag(A.y), 
        (real_t)real(A.z), (real_t)imag(A.z));
}

void print(const vector_t<real_t> A){
    printf("(%21.14E, %21.14E, %21.14E)\n", (real_t)A.x, (real_t)A.y, (real_t)A.z);
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

real_t mag(const vector_t<complex_t> A){
    return sqrt(abs(A*A));
}

vector_t<real_t> unit(const vector_t<real_t> A){
    real_t A_mag=mag(A);
    return vector_t<real_t>(A.x/A_mag, A.y/A_mag, A.z/A_mag);
}

int is_equal(const vector_t<real_t> A, const vector_t<real_t> B, const real_t tol){
    if (abs(A.x-B.x)<=abs(A.x*tol) && abs(A.y-B.y)<=abs(A.y*tol) && abs(A.z-B.z)<=abs(A.z*tol)){
        return true;
    }else{
        return false;
    }
}

vector_t<real_t> real_v(const vector_t<complex_t> A){
    return vector_t<real_t>(real(A.x), real(A.y), real(A.z));
}

vector_t<real_t> imag_v(const vector_t<complex_t> A){
    return vector_t<real_t>(imag(A.x), imag(A.y), imag(A.z));
}