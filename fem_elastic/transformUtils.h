
#ifndef h_gmp_utils_transform_h
#define h_gmp_utils_transform_h

float* read_transform(const char* fname);

void write_transform(float* transform, const char* fname);

void inv_transform(float* t, float** pinv_t);

#endif
