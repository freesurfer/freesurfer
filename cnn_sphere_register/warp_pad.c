/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */



/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/



/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <string.h>
//#include <fstream>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"

#include "mrisurf.h"
#include "mri.h"
#include "matfile.h"
#include "matrix.h"
//#include <direct.h>

//using namespace std;

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

bool flag = false;

char *TempVolFile=NULL;
char *subject, *hemi, *SUBJECTS_DIR;

//#define U_DIM 256
//#define V_DIM 512

float resampleatpolarcoord(float **mtx, float theta, float phi)
{
  float uf, vf, du, dv, val;
  int u0, u1, v0, v1;

  uf = 256.0f * phi / PHI_MAX;
  vf = 512.0f * theta / THETA_MAX;
  u0 = floor(uf);
  u1 = ceil(uf);
  v0 = floor(vf);
  v1 = ceil(vf);  // name = */
  du = uf - (float)u0;
  dv = vf - (float)v0;




  if (u0 < 0)  u0 = -u0;
  if (u0 >= 256) u0 = 256 - (u0 - 256 + 1);
  if (u1 < 0) /* enforce spherical topology  */
    u1 = -u1;
  if (u1 >= 256) u1 = 256 - (u1 - 256 + 1);
  if (v0 < 0) v0 += 512;
  if (v0 >= 512) v0 -= 512;
  if (v1 < 0) v1 += 512;
  if (v1 >= 512) v1 -= 512;


  //float *vals = mf->data;

  //val = du * dv * *MATRIX_RELT(mtx, v1, u1) + (1.0f - du) * dv * *MATRIX_RELT(mtx, v1, u0) +
  //       (1.0f - du) * (1.0f - dv) * *MATRIX_RELT(mtx, v0, u0) + du * (1.0f - dv) * *MATRIX_RELT(mtx, v0, u1);

         val = du * dv * mtx[v1][u1] + (1.0f - du) * dv * mtx[v1][u0] +
                (1.0f - du) * (1.0f - dv) * mtx[v0][u0] + du * (1.0f - dv) * mtx[v0][u1];

                if (flag) {printf("u0=%d, u1=%d, v0=%d, v1=%d\n", u0, u1, v0, v1); printf("du=%f, dv=%f\n, val=%f", du, dv, val);flag = false;}

  return val;
}

/*float **getdata(MATFILE *mf)
{

  float **mtx = calloc(mf->mrows * mf->ncols, sizeof(float));
  float *fptr = (float *)mf->data;
  //mtx = matAlloc((int)mf->mrows, (int)mf->ncols);

  //int i;
  //for (i = 0; i < mf->ncols; i++) mtx[i] = new float[mf->ncols];

  int row, col;
  for (row = 0; row < mf->mrows; row++) {
    for (col = 0; col < mf->ncols; col++) {
      mtx[row][col] = *fptr++ ;
    }
  }

  return mtx;
}*/

float **readtxt(const char *fname, int rows, int cols)
{
  /*float **mtx = new float*[rows];
  int i;
  for (i= 0; i<cols; i++) mtx[i] = new float[cols];*/
  //float mtx[512][256];
  //MATRIX *mtx = matAlloc(rows, cols);
  float **mtx;
  mtx = (float **)calloc(rows, sizeof(float *));
  int i;
  for (i= 0; i<rows; i++) mtx[i] = (float *)calloc(cols, sizeof(float));

  FILE *fp = fopen(fname, "rb");
  //ifstream fs(fname);
  printf("alloc room for mtx\n" );

  int u,v;
  for (u = 0; u < rows; u++)
    for (v = 0; v< cols; v++)
      fscanf(fp, "%f", &mtx[u][v]);


  fclose(fp);

  return mtx;
}


/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {



  //char *fname[STRLEN];
  MRI_SURFACE  *mris ;
  //MRI_SP       *mrisp1, *mrisp2;
  //MATFILE *ux = NULL, *uy = NULL;


  //fprintf(stderr, "processing subject %s...\n", ar, hemi, surf_name) ;

  mris = MRISread(argv[1]) ;

  float **ux=NULL, **uy=NULL;

  ux = readtxt(argv[2],512,256);
  uy = readtxt(argv[3],512,256);
  


  //deformation field by rotate 90 around y


    int vno;
    float x, y, z, radius, a, b, c, d, delta_theta, delta_phi, theta, phi , newtheta, newphi;
    VERTEX *vertex/*, *ver_distheta, *ver_disphi*/;

    radius = a = b = c = MRISaverageRadius(mris);

    printf("number of vertex is %d\n", mris->nvertices);

    for (vno = 0; vno < mris->nvertices; vno++) {
      vertex = &mris->vertices[vno];


      if (vno==100) flag = true;

      x = vertex->x;
      y = vertex->y;
      z = vertex->z;
      theta = atan2(vertex->y / b, vertex->x / a);
      if (theta < 0.0f) theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
      d = c * c - z * z;
      if (d < 0.0) d = 0.0;
      phi = atan2(sqrt(d), z);

      //if (fabsf(phi-M_PI/2.0f)<1.4f)
      if (1)
      {
        delta_theta = resampleatpolarcoord(ux, theta, phi);
        delta_phi = resampleatpolarcoord(uy, theta, phi);

        newtheta = theta+delta_theta*2.0f*M_PI/512.0f;
        newphi = phi + delta_phi*M_PI/256.0f;



        // enforce spherical topology
        if (newphi < 0.0f) newphi = -newphi;
        if (newphi > M_PI) newphi = M_PI-(newphi-M_PI);

  //      if (phi > M_PI) { phi = 2*M_PI-phi; theta = theta + M_PI;}

        vertex->phi = newphi;

        if (newtheta < 0.0f)  newtheta = 2.0f * M_PI + newtheta;
        else if (newtheta >= 2.0f*M_PI) newtheta -= 2.0f*M_PI;
        //if (phi < 0.0f) phi = -phi;


        //vertex->phi = phi;
        vertex->theta = newtheta;

        vertex->x = radius * sin(newphi) * cos(newtheta);
        vertex->y = radius * sin(newphi) * sin(newtheta);
        vertex->z = radius * cos(newphi);
      }
     else //high-latitude region
      {
        float x1, y1, z1, theta1, phi1, delta_theta1, delta_phi1, newtheta1, newphi1;
        x1 = -z; y1 = y; z1 = x; // rotate 90 degree around y axis
        theta1 = atan2(y1 / b, x1 / a);
        if (theta1 < 0.0f) theta1 = 2 * M_PI + theta1;
        d = c * c - z1 * z1;
        if (d < 0.0) d = 0.0;
        phi1 = atan2(sqrt(d), z1);

        delta_theta1 = resampleatpolarcoord(ux, theta1, phi1);
        delta_phi1 = resampleatpolarcoord(uy, theta1, phi1);

        newtheta1 = theta1+delta_theta1*2.0f*M_PI/512.0f;
        newphi1 = phi1 + delta_phi1*M_PI/256.0f;

        // enforce spherical topology
        if (newphi1 < 0.0f) newphi1 = -newphi1;
        if (newphi1 > M_PI) newphi1 = M_PI-(newphi1-M_PI);

  //      if (phi > M_PI) { phi = 2*M_PI-phi; theta = theta + M_PI;}

        //vertex->phi = newphi;

        if (newtheta1 < 0.0f)  newtheta1 = 2.0f * M_PI + newtheta1;
        else if (newtheta1 >= 2.0f*M_PI) newtheta1 -= 2.0f*M_PI;

        x1 = radius * sin(newphi1) * cos(newtheta1);
        y1 = radius * sin(newphi1) * sin(newtheta1);
        z1 = radius * cos(newphi1);

        vertex->x = z1;
        vertex->y = y1;
        vertex->z = -x1;

      }


    }



    char *outname;
    outname = argv[4];

    MRISwrite(mris, outname) ;



  if (mris) MRISfree(&mris) ;

  //MRISPfree(&mrisp1);
  //MRISPfree(&mrisp2);

  return 0;
}
/* ------ Doxygen markup starts on the line below ---- */
/*!
\fn int parse_commandline(int argc, char **argv)
\brief Parses the command-line arguments
\param argc - number of command line arguments
\param argv - pointer to a character pointer
*/
/* ------ Doxygen markup ends on the line above ---- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--temp-vol")) {
      if (nargc < 1) CMDargNErr(option,1);
      TempVolFile = pargv[0];
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void usage_exit(void)
\brief Prints usage and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_usage(void)
\brief Prints usage and returns (does not exit)
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --temp-vol volfile : template volume \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_help(void)
\brief Prints help and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_version(void)
\brief Prints version and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void check_options(void)
\brief Checks command-line options
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void check_options(void) {
  return;
}

/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void dump_options(FILE *fp)
\brief Prints command-line options to the given file pointer
\param FILE *fp - file pointer
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  return;
}
