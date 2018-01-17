// See https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
//

static float fnv_hash(int i, int k, int* random_counter, double LO, double HI) {

    unsigned int combined[3]; combined[0] = i; combined[1] = k; combined[2] = (*random_counter)++; 
    unsigned char* p = (unsigned char*)&combined[0];

    unsigned long hash = 0x811c9dc5;
    int c;
    for (c = 0; c < sizeof(combined); c++) {
        hash ^= p[c];
        hash *= 16777619L;
    }

    double f = (double)(hash & (unsigned long)0xffffffff) 
             / (double)        (unsigned long)0xffffffff;

    return (float)((HI-LO)*f + LO);
}
