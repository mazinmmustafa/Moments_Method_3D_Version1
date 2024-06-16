//
#include "utilities.hpp"

void print(const char *format, ...){
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}

void print(const int_t n){
    print("%d\n", n);
}

void print(const size_t n){
    print("%zu\n", n);
}

void print(const real_t x){
    print("%21.14E\n", (real_t)x);
}

void print(const complex_t z){
    print("(%21.14E, %21.14E)\n", real(z), imag(z));
}

void __assert_error(const int_t condition, const char *error_msg, const char* filename, const size_t line){
    if (condition!=true){
        print("%s:%zu: error: %s!\n", filename, line, error_msg);
        exit(1);
    }
}

void progress_bar(const size_t i, const size_t N, const char *msg){
    const size_t N_bar=15;
    size_t n=i;
    n++;
    print("\rProgress: [");
    for (size_t i=0; i<=n*N_bar/N; i++){
        print("#");
    }
    for (size_t i=n*N_bar/N; i<N_bar; i++){
        print("-");
    }
    print("] %3zu%%,", 100*n/N);
    print(" %s", msg);
    if (n==N){
        print(", done!\n");
    }else
    if (n>N){
        print("\n");
    }
    fflush(stdout);
}

//

range_t::range_t(){

}

range_t::~range_t(){
    range_t::unset();
}

void range_t::linspace(){
    assert_error(this->is_allocated, "range info are not set");
    const real_t dx=(this->x_max-this->x_min)/(this->Ns-1.0);
    for (size_t i=0; i<this->Ns; i++){
        this->data[i] = this->x_min+i*dx;
    }
}

void range_t::logspace(){
    assert_error(this->is_allocated, "range info are not set");
    assert_error(this->x_min>0.0, "invalid range");
    const real_t dx=(log10(this->x_max)-log10(this->x_min))/(this->Ns-1.0);
    print(this->x_min);
    print(this->x_max);
    for (size_t i=0; i<this->Ns; i++){
        this->data[i] = pow(10.0, log10(this->x_min)+i*dx);
    }
}

real_t range_t::operator() (const size_t index) const{
    assert_error(index<this->Ns, "invalid range index");
    return this->data[index];
}

void range_t::set(const real_t x_min, const real_t x_max, const size_t Ns){
    if (this->is_allocated){
        range_t::unset();
    }
    assert_error(Ns>0, "invalid number of samples");
    this->Ns = Ns;
    assert_error(x_max>x_min, "invalid range index");
    this->x_min = x_min;
    this->x_max = x_max;
    this->data = (real_t*)calloc(Ns, sizeof(real_t));
    assert(this->data!=null);
    this->is_allocated = true;
}

void range_t::unset(){
    if (this->is_allocated){
        free(this->data);
        this->Ns = 0;
    }
}

void range_t::get_info(real_t *x_min, real_t *x_max, size_t *Ns){
    (*x_min) = this->x_min;
    (*x_max) = this->x_max;
    (*Ns) = this->Ns;
}

