/////////////////////////////////////////////////////////////////////////
/* Bruker.c                               */
/* created by : y.tosa                    */
/* date       :8/27/2003                  */
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: tosa $
// Revision Date  : $Date: 2003/08/29 16:44:42 $
// Revision       : $Revision: 1.1 $

// there are many files present in Bruker directory
//
// fid             ... raw data
// method          ... store similar info as acqp  <--- used
// acqp            ... used for reconstruction     <--- used
// pulseprogram    ...
// spnam()         ...
// gdprog.ax, ay, az . 
// pdata ... reconstructed image strorage
//   |
//   |-> 1   
//       |-> 2dseq  ... bshort image                                   <-- used
//           reco   ... created after reconstruction
//           d3proc ... reconstruted image info (width, height, depth) <-- used
//           procs  ... ???
//           roi    ... ???
//
char *BRUCKER_C_VERSION= "$Revision: 1.1 $";

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h> 
#include <math.h>
#include "mri.h"
#include "Bruker.h"

MRI *brukerRead(char *fname, int read_volume)
{
  char methodFile[1024];
  char acqpFile[1024];
  char dataFile[1024];
  char d3procFile[1024];

  int succeed = 0;
  MRI *mri = 0;

  // first make sure that method, acqp, fid exist
  // use stat to do this.
  if ((succeed = checkBrukerFiles(fname, methodFile, acqpFile, dataFile, d3procFile,1))==0)
    return (MRI*) 0;
 
  // next read method to get TE, TR, volume size, pixel size
  if ((mri = readBrukerMethod(methodFile, d3procFile))==0)
  {
    if (mri)
      MRIfree(&mri);
    return (MRI *) 0;
  }
  // next read acqp to get missing info like direction cosine
  if ((succeed = readBrukerAcqp(mri, acqpFile))==0)
  {
    if (mri)
      MRIfree(&mri);
    return (MRI *) 0;
  }
  // now ready to read volume
  if (read_volume)
  {
    if ((succeed = readBrukerVolume(mri, dataFile))==0)
    {
      if (mri)
	MRIfree(&mri);
      return (MRI *) 0;
    }
  }
  return mri;
}

int checkBrukerFiles(char *fname, char *methodFile, char *acqpFile, char *dataFile, char *d3procFile, int flag)
{
  struct stat stat_buf;
  
  if (stat(fname, &stat_buf) < 0)
  {
    if (flag)
      fprintf(stderr, "ERROR: could not stat %s.\n", fname);
    return 0;
  }
  // first check fname is a directory  
  if (S_ISDIR(stat_buf.st_mode) ==0)
  {
    if (flag)
      fprintf(stderr, "ERROR: %s is not a directory.\n", fname);
    return 0;
  }
  if (flag)
    printf("INFO: directory %s \n", fname);
  // build up methodFile name
  strcpy(methodFile, fname);
  strcat(methodFile, "/method");
  if (stat(methodFile, &stat_buf) < 0)
  {
    if (flag)
      fprintf(stderr, "ERROR: could not stat %s.\n", methodFile);
    return 0;
  }
  if (flag)
    printf("INFO:    method %s\n", methodFile);
  // build up acqpFile name
  strcpy(acqpFile, fname); 
  strcat(acqpFile, "/acqp");
  if (stat(acqpFile, &stat_buf) < 0)
  {
    if (flag)
      fprintf(stderr, "ERROR: could not stat %s.\n", acqpFile);
    return 0;
  }
  if (flag)
    printf("INFO:      acqp %s\n", acqpFile);
  // build up dataFile name

  strcpy(dataFile, fname);
  strcat(dataFile, "/pdata/1/2dseq");
  if (stat(dataFile, &stat_buf) < 0)
  {
    if (flag)
      fprintf(stderr, "ERROR: could not stat %s.\n", dataFile);
    return 0;
  }
  if (flag)
    printf("INFO:     2dseq %s \n", dataFile);
  // build up d3proc file
  strcpy(d3procFile, fname);
  strcat(d3procFile, "/pdata/1/d3proc");
  if (stat(dataFile, &stat_buf) < 0)
  {
    if (flag)
      fprintf(stderr, "ERROR: could not stat %s.\n", d3procFile);
    return 0;
  }
  if (flag)
    printf("INFO:    d3proc %s \n", d3procFile);
  return 1;
}

int separate_parameter_and_value(char *sWholeLine, char *sParameter, char *sValue)
{
  char *P0, *P1;

  P0 = strstr(sWholeLine,"##");
  if ( !P0 ) return(1); /* ignore line */
  P0+=2;                        /* advance past ## */
  strcpy(sParameter,P0);        /* initialize sParameter */
  P0 = sParameter;              /* reset P0 */
  while ( *P0 != '=' ) P0++;    /* search for '=' */
  *P0 = '\0';                   /* terminate string */
  P0++;                         /* step past '=' */
  strcpy(sValue,P0);            /* initialize sValue */
  /* Use P1 (sValue) to copy rest of line (sParameter string) */
  for (P1=sValue; *P0; P0++, P1++)
    *P1 = *P0;
  /* Eliminate parentheses and CR in the value string. */
  for (P0=sValue; *P0; P0++)
    if ( *P0=='(' || *P0==')' || *P0=='\n' ) *P0=' ';

  return(0);
}

int readBrukerD3proc(char *d3procFile, int *px, int *py, int *pz, int *ptype, int *pnframes)
{
  FILE *fp=0;
  char line[512];
  char Value[128];
  char Parameter[256];
  int lRead;
  int ignore;

  fp = fopen(d3procFile,"r");
  if (fp ==0)
  {
    fprintf(stderr, "ERROR: could not open d3proc %s", d3procFile);
    return 0;
  }
  while (fgets(line, sizeof(line), fp))
  {
    if ((ignore = separate_parameter_and_value(line, Parameter, Value)))
      continue;

    // now gets the values
    if ( !strcmp(Parameter,"END") )
      break;
    // get volume size
    else if ( !strcmp( Parameter, "$IM_SIX") )
      lRead = sscanf(Value, "%d", px);
    else if ( !strcmp( Parameter, "$IM_SIY") )
      lRead = sscanf(Value, "%d", py);
    else if ( !strcmp( Parameter, "$IM_SIZ") )
      lRead = sscanf(Value, "%d", pz);
    else if ( !strcmp( Parameter, "$IM_SIT") )
    {
      lRead = sscanf(Value, "%d", pnframes);
      if (*pnframes > 1)
      {
	fprintf(stderr, "ERROR: nframes %d but one is supported.\n", *pnframes);
	return 0;
      }
    }
    else if ( !strcmp( Parameter, "DATTYPE") )
    {
      if (strcmp(Value, "ip_short")==0)
	*ptype = MRI_SHORT;
      else
      {
	fprintf(stderr, "ERROR: currently only short type is supported.\n");
	return 0;
      }
    }
  }
  fclose(fp);
  return 1;
}

MRI* readBrukerMethod(char *methodFile, char *d3procFile)
{
  MRI *mri = 0;
  FILE *fp = 0;
  char line[512];
  char Parameter[256];
  char Value[128];
  int ignore;
  int width=0;
  int height=0;
  int depth=0;
  double xsize, ysize, zsize;
  double TE=0;
  double TR=0;
  int type = MRI_SHORT;
  int nframes=0;
  int lDim=0;
  int lRead=0;

  if (!readBrukerD3proc(d3procFile, &width, &height, &depth, &type, &nframes))
    return (MRI *) 0;

  fp = fopen(methodFile, "r");
  if (fp ==0)
  {
    fprintf(stderr, "ERROR: could not open methodFile %s", methodFile);
    return 0;
  }
  while (fgets(line, sizeof(line), fp))
  {
    if ((ignore = separate_parameter_and_value(line, Parameter, Value)))
      continue;

    // now gets the values
    if ( !strcmp(Parameter,"END") )
      break;

    // get volume size
    else if ( !strcmp( Parameter,"$PVM_Matrix") )
    {
      lRead = sscanf( Value,"%d",&lDim);
      if ( lDim != 3 )
      {
	printf("INFO: Matrix is not 3 but %d.  The image is not of a volume image\n", lDim);
      }
      if (!fgets(line, sizeof(line), fp))
      {
	printf("Error: dimensions should always follow $PVM_Matrix.\n");
	fclose(fp);
	return 0;
      }
      if (lDim == 3)
	lRead = sscanf(line,"%d %d %d", &height, &width, &depth);
      else if (lDim == 2)
      {
	lRead = sscanf(line,"%d %d", &height, &width);
	depth = 3;
      }
    }
    else if ( !strcmp(Parameter,"$PVM_Fov") )
    {
      lRead = sscanf(Value,"%d",&lDim);
      if ( lDim != 3 )
      {
	printf("INFO: Fov is not 3 but %d.  The image is not of a volume image\n", lDim);
      }
      if (!fgets(line, sizeof(line), fp))
      {
	fprintf(stderr, "Error: FOV should always follow $PVM_Fov.\n");
	fclose(fp);
	return 0;
      }
      // note y, x, z order! 
      if (lDim == 3)
	lRead = sscanf(line,"%lf %lf %lf", &ysize, &xsize, &zsize);
      else if (lDim == 2)
      {
	lRead = sscanf(line,"%lf %lf", &ysize, &xsize);
	zsize = xsize; // fake
      }
    }
    else if ( !strcmp(Parameter,"$PVM_EchoTime") )
    {
      lRead = sscanf(Value,"%lf",&TE);
    }
    else if ( !strcmp(Parameter,"$PVM_RepetitionTime") )
    {
      lRead = sscanf(Value,"%lf",&TR);
    }
  }
  fclose(fp);

  if (width > 0 && height > 0 && depth > 0)
  {
    mri = MRIalloc(width, height, depth, type);

    // get real sizes
    xsize /= (double) width;
    ysize /= (double) height;
    zsize /= (double) depth; // in mm

    // TR, TE, flip angle overridden by acqp
    mri->tr = TR;
    mri->te = TE;
    mri->nframes = nframes;
    mri->xsize = xsize;
    mri->ysize = ysize;
    mri->zsize = zsize;
    return mri;
  }
  else
    return 0;
}

int readBrukerAcqp(MRI *mri, char *acqpFile)
{
  FILE *fp = 0;
  char line[512];
  char Parameter[256];
  char Value[128];
  int ignore=0;
  int lRead=0;
  int dim=0;
  double flip_angle=0;
  double TI=0;
  double TR=0;
  double TE=0;
  double x_r, x_a, x_s;
  double y_r, y_a, y_s;
  double z_r, z_a, z_s;
  double c_r, c_a, c_s;

  x_r = x_a = x_s = 0;
  y_r = y_a = y_s = 0;
  z_r = z_a = z_s = 0;


  if (!mri)
  {
    fprintf(stderr, "ERROR: readBrukerMethod() must be called before readBrukerAcqp");
    return 0;
  }

  fp = fopen(acqpFile, "r");
  if (fp ==0)
  {
    fprintf(stderr, "ERROR: could not open acqpFile %s", acqpFile);
    return 0;
  }
  while (fgets(line, sizeof(line), fp))
  {
    if ((ignore = separate_parameter_and_value(line, Parameter, Value)))
      continue;

    // now gets the values
    if ( !strcmp(Parameter,"END") )
      break;
    else if ( !strcmp( Parameter, "OWNER") )
    {
      fgets(line, sizeof(line), fp);
      printf("INFO: %s", line);
      fgets(line, sizeof(line), fp);
      printf("INFO: %s", line);
    }
    // another check for dimension
    else if ( !strcmp( Parameter,"$ACQ_dim") )
    {
      lRead = sscanf( Value, "%d", &dim);
      if (dim != 3)
      {
	printf("INFO: acqp tells the dimension is not 3 but %d\n", dim);
      }
    }
    // grad_matrix
    else if ( !strcmp( Parameter,"$ACQ_grad_matrix") )
    {
      if (!fgets(line, sizeof(line), fp))
      {
	fprintf(stderr, "ERROR: float value must follow ACQ_grad_matrix");
	fclose(fp);
	return 0;
      }
      // at this time I don't know the relation between grad-matrix and direction cosines
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	     &x_r, &x_a, &x_s,
	     &y_r, &y_a, &y_s,
	     &z_r, &z_a, &z_s);
    }
    // flip angle
    else if ( !strcmp( Parameter, "$ACQ_flip_angle") )
    {
      lRead=sscanf(Value,"%lf", &flip_angle);
      // convert into radians
      flip_angle = flip_angle*3.141592653589793/180.;
    }
    // TR
    else if ( !strcmp( Parameter, "$ACQ_repetition_time") )
    {
      if (!fgets(line, sizeof(line), fp))
      {
	fprintf(stderr, "ERROR: float value must follow ACQ_repetition_time");
	fclose(fp);
	return 0;
      }
      sscanf(line, "%lf", &TR);
    }
    // TE
    else if ( !strcmp( Parameter, "$ACQ_echo_time") )
    {
      if (!fgets(line, sizeof(line), fp))
      {
	fprintf(stderr, "ERROR: float value must follow ACQ_echo_time");
	fclose(fp);
	return 0;
      }
      sscanf(line, "%lf", &TE);
    }
    // TI
    else if ( !strcmp( Parameter, "$ACQ_inversion_time") )
    {
      if (!fgets(line, sizeof(line), fp))
      {
	fprintf(stderr, "ERROR: float value must follow ACQ_inversion_time");
	fclose(fp);
	return 0;
      }
      sscanf(line, "%lf", &TI);
    }
  }
  fclose(fp);
  // override method values
  mri->flip_angle = flip_angle;
  mri->ti = TI;
  mri->te = TE;
  mri->tr = TR;
  // at this time I don't know the relation between grad_matrix and direction cosines
  // use information from George Dai
  // that is, 
  //   grad_matrix         corresponds to 
  //    0 0 1 1 0 0 0 1 0      his real world  x -> -R, y->S, z -> A           
  //                           image produced  x -> -R, y -> -S, z -> A
  c_r = c_a = c_s = 0;
  mri->x_r = -1; mri->x_a = 0, mri->x_s = 0;
  mri->y_r =  0; mri->y_a = 0, mri->y_s = -1;
  mri->z_r =  0; mri->z_a = 1, mri->z_s = 0;
  mri->c_r = c_r; mri->c_a = c_a, mri->c_s = c_s;
  mri->ras_good_flag = 1;
  return 1;
}

#ifdef Linux
extern void swab(const void *from, void *to, size_t n);
#endif

int readBrukerVolume(MRI *mri, char *dataFile)
{
  FILE *fp = 0;
  int k,j;
  int nread;
  int swap_bytes_flag = 0;

  if (!mri)
  {
    fprintf(stderr, "ERROR: readBrukerMethod() must be called before readBrukerVolume");
    return 0;
  }
  // save the data filename
  strcpy(mri->fname, dataFile);
  fp = fopen(dataFile, "r");
  if (fp ==0)
  {
    fprintf(stderr, "ERROR: could not open dataFile %s", dataFile);
    return 0;
  }
  for (k = 0; k < mri->depth; ++k)
    for (j = 0; j < mri->height; ++j)
    {
      nread = fread(mri->slices[k][j], sizeof(short), mri->width, fp);
      if (nread != mri->width)
      {
	fclose(fp);
	MRIfree(&mri);
	mri = 0;
	return 0;
      }
      // this was not needed
      if(swap_bytes_flag)
      {
	swab(mri->slices[k][j], mri->slices[k][j], mri->width *sizeof(short));
      }
    }

  fclose(fp);
  return 1;
}



