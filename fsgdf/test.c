#include "fsgdf.h"

char *Progname="fsgdf test app";

int main (int argc, char** argv)
{
  FSGD *gd;
  char *env;
  char fnTestPath[1000];
  char fnTest[1000];

  /* Build a proper path for our test data. */
  env = getenv("FSDEV_TEST_DATA");
  if(NULL != env)
    {
      strcpy(fnTestPath, env);
    }
  else 
    {
      strcpy(fnTestPath, "/space/lyon/1/fsdev/test_data");
    }
  sprintf(fnTest, "%s/fsgdf/y-lh.fsgd", fnTestPath);

  /* Just test the read function, first without reading the data and
     then with. */
  gd = gdfRead(fnTest,0);
  if(NULL == gd)
    {
      printf("ERROR: gdfRead(%s,0) returned NULL",fnTest);
      exit(1);
    }
  gdfFree(&gd);
  gd = gdfRead(fnTest,1);
  if(NULL == gd)
    {
      printf("ERROR: gdfRead(%s,1) returned NULL",fnTest);
      exit(1);
    }
  gdfFree(&gd);


  return(0);
}
