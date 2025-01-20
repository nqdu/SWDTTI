#include <iostream>

inline void __myfwrite(const void *__ptr, size_t __size, size_t __nitems, FILE *__stream)
{
    size_t size = fwrite(__ptr,__size,__nitems,__stream);
    if(size != __nitems) {
        fprintf(stderr,"cannot write to binary!\n");
        exit(1);
    }
}


template<typename T>
void 
write_binary_f(FILE *fp, const T *data, size_t n)
{
    // write integers of the size 
    int size = (int)(n * sizeof(T));

    // integer front
    __myfwrite(&size,sizeof(int),1,fp);

    // data
    __myfwrite(data,sizeof(T),n,fp);

    // integer back
    __myfwrite(&size,sizeof(int),1,fp); 
}