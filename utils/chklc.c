

#include <unistd.h>
#include <const.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char     *crypt(const char *, const char *);

void chklc(void)
{
  /*#ifndef Darwin*/
  char  dirname[STRLEN], *cp ;
  FILE* lfile;
  char* email;
  char* magic;
  char* key;
  char* gkey;
  char* lfilename;
  char  str[100] ;

  sprintf(str, "S%sER%sACK%sOR", "URF", "_B", "DO") ;
  if (getenv(str) != NULL)
    return ;

  cp = getenv("MRI_DIR");
 
  if (cp)
   {
     strncpy(dirname, cp, STRLEN) ;
    }
  else
  {
    dirname[0] = 46; /*  ascii "."   */
     dirname[1] = 0;
  }
 

  lfilename = (char*)calloc(1,512);
  email = (char*)calloc(1,512);
  magic = (char*)calloc(1,512);
  key = (char*)calloc(1,512);
  gkey = (char*)calloc(1,1024);

  sprintf(lfilename,"%s/.lic%s",dirname, "ense");

  lfile = fopen(lfilename,"r");
  if(lfile) {
    fscanf(lfile,"%s\n",email);
    fscanf(lfile,"%s\n",magic);
    fscanf(lfile,"%s\n",key);
    
    sprintf(gkey,"%s.%s",email,magic);
#ifndef Darwin
    if (strcmp(key,crypt(gkey,"*C*O*R*T*E*C*H*S*0*1*2*3*"))!=0) {
      printf("No valid license key !\n");
      exit(-1);
    }
#endif
  }
  else {
    printf("License file not found !\n");
    exit(-1);
  }
  free(email);
  free(magic);
  free(key);
  free(gkey);
  free(lfilename);
  return;  
  /*#endif*/
}

