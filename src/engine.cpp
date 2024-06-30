//
#include "engine.hpp"



engine_t::engine_t(){

}

engine_t::~engine_t(){

}



//

const real_t eps_sinc=1.0E-15;

complex_t sinc(const complex_t x){
    return abs(x)<eps_sinc ? 1.0 : sin(x)/x;
}

real_t sinc(const real_t x){
    return abs(x)<eps_sinc ? 1.0 : sin(x)/x;
}

