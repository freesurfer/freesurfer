#ifndef HMEM_H
#define HMEM_H
                            
#ifndef _MSDOS                            
#define hmemset   memset
#define hmemcpy   memcpy
#define hcalloc   calloc
#else 
                                               
#define hmemset(buf, val, count)   huge_memset((huge char *)buf, (char)val, (long)count)
#define hcalloc(n,size)            huge_alloc((hsize_t)n*size, 1)
#define hmemcpy(Dst, Src, Count)   huge_memcpy((huge char *)Dst, (huge char *)Src, (long)Count)

#endif

#endif
