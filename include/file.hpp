#ifndef __FILE_HPP__
#define __FILE_HPP__

// Libraries
#include "lib_basic.hpp"
#include "utilities.hpp"

// Definitions
class file_t{
    protected:
        FILE *file_ptr=NULL;
        char mode='w';
        const size_t max_length=150;
        char *filename=NULL;
        int is_open=false;
    public:
        file_t();
        ~file_t();
        void open(const char *filename, const char mode);
        void close();
        void write(const char *format, ...);
        int read(const char *format, ...);
};

class binary_file_t : public file_t{
    private:
    public:
        binary_file_t(){}
        ~binary_file_t(){
            binary_file_t::close();
        }
        void open(const char *filename, const char mode);
        template <typename type_t>
        int write(const type_t *data){
            return fwrite(data, sizeof(*data), 1, this->file_ptr);
        }
        template <typename type_t>
        int read(type_t *data){
            return fread(data, sizeof(*data), 1, this->file_ptr);
        }
};

// Functions


#endif

