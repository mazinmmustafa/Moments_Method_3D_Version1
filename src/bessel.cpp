//
#include "bessel.hpp"


// Fortran functions
extern "C" void zbesj_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesy_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesi_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesk_f77_(double *fnu, double *zr, double *zi, double *cyr, double *cyi);
extern "C" void zbesh_f77_(int *k, double *fnu, double *zr, double *zi, double *cyr, double *cyi);

complex_t besselj(const real_t n, const complex_t z){
	real_t zr=real(z);
	real_t zi=imag(z);
	real_t cyr, cyi;
	complex_t j=complex_t(0.0, 1.0);
	real_t n_=n;
	zbesj_f77_(&n_, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t bessely(const real_t n, const complex_t z){
	real_t zr=real(z);
	real_t zi=imag(z);
	real_t cyr, cyi;
	complex_t j=complex_t(0.0, 1.0);
	real_t n_=n;
	zbesy_f77_(&n_, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t besseli(const real_t n, const complex_t z){
	real_t zr=real(z);
	real_t zi=imag(z);
	real_t cyr, cyi;
	complex_t j=complex_t(0.0, 1.0);
	real_t n_=n;
	zbesi_f77_(&n_, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t besselk(const real_t n, const complex_t z){
	real_t zr=real(z);
	real_t zi=imag(z);
	real_t cyr, cyi;
	complex_t j=complex_t(0.0, 1.0);
	real_t n_=n;
	zbesk_f77_(&n_, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t besselh(const int k, const real_t n, const complex_t z){
	assert(k==1||k==2);
	real_t zr=real(z);
	real_t zi=imag(z);
	real_t cyr, cyi;
	complex_t j=complex_t(0.0, 1.0);
	int k_=k;
	real_t n_=n;
	zbesh_f77_(&k_, &n_, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}
