/*
  @(#)runfuncs.h  1.1
  4/4/94
*/
/*------------------------------------------------------------------------
      File Name:  runfuncs.h

         Author:  Bruce Fischl

        Created:  Jan. 1993

    Description:  

------------------------------------------------------------------------*/
#ifndef RUNFUNCS_H
#define RUNFUNCS_H

typedef int (*fwd_func)(unsigned char **outPtr) ;
typedef void (*inv_func)(unsigned char **outPtr, unsigned char val) ;

void runFuncInit(fwd_func *fwd_array, inv_func *inv_array) ;


#endif
