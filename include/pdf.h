
#ifndef PDF_INC
#define PDF_INC

unsigned long PDFtodSeed(void);
double PDFgaussian(void);
double PDFerlang(int order);
int PDFloadCDF(char *fname, double **xcdf, double **cdf, int *ncdf);
double PDFsampleCDF(double *xcdf, double *cdf, int ncdf);


#endif //#ifndef PDF_INC
