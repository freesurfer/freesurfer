#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DLEN 1000

int
main(int argc, char *argv[])
{
  double result ;
  int    i ;

  srand48(0L) ;

  for (result = 0.0, i = 0 ; i < DLEN ; i++)
    result += log(drand48()) ;

  result = exp(-(result/(DLEN/8))) ;
  printf("result = %f\n", result) ;
  return(0) ;
}
  
