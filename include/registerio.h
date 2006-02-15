/**********************************************************
Routines for handling register.dat files
***********************************************************/


#ifndef REGISTERIO_H_INC
#define REGISTERIO_H_INC

#include <stdio.h>
#include "matrix.h"

int regio_read_register(char *regfile, char **subject, float *inplaneres, 
      float *betplaneres, float *intensity,  MATRIX **R,
      int *float2int);

int regio_print_register(FILE *fp, char *subject, float inplaneres, 
       float betplaneres, float intensity, MATRIX *R,
       int float2int);

int regio_write_register(char *regfile, char *subject, float inplaneres, 
       float betplaneres, float intensity, MATRIX *R, 
       int float2int);

int regio_read_mincxfm(char *xfmfile, MATRIX **R, char **fileinfo);
int regio_write_mincxfm(char *xfmfile, MATRIX *R, char *fileinfo);


int regio_read_xfm4(char *xfmfile, MATRIX **R);
int regio_read_xfm(char *xfmfile, MATRIX **R);


#endif /*BF_H_INC*/
