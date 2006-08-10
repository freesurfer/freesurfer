#include <stdlib.h>
#include <time.h>

#include "nr_wrapper_open_source.h"

int main(int argc, char *argv[]) {
  // just trying to make sure that we can call the functions from c
  long seed = -1L * (long)( abs( (int)time(NULL) ) );  
  
  printf("random number: %f\n", OpenRan1( &seed ));
  
  return 0;  
}
