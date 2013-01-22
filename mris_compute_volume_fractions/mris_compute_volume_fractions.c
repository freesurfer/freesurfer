/**
 * @file  mris_compute_volume_fractions.c
 * @brief Computes accurate volume fractions remaining within a surface
 *
 * This program computes an accurate estimate of the fraction of the volume 
 * remaining wihtin a surface. 
 */
/*
 * Original Author: Ender Konukoglu
 * CVS Revision Info:
 *    $Author: enderk $
 *    $Date: 2013/01/22 16:14:14 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mris_compute_volume_fractions.h"

MRI*
MRIcomputeVolumeFractionFromSurface(MRI_SURFACE *mris, double acc, MRI *mri_src, MRI *mri_fractions)
{
  const int width = mri_src->width;
  const int height = mri_src->height;
  const int depth = mri_src->depth;
  int x,y,z, vno;
  double xs, ys, zs, dist; 
  MRIS_HASH_TABLE *mht;
  VERTEX *v; 
  
  /* preparing the output */
  printf("preparing the output\n");
  if (mri_fractions == NULL)
    {
      mri_fractions = MRIalloc(width,height,depth,MRI_FLOAT);
      MRIcopyHeader(mri_src, mri_fractions);
    }
  MRI *mri_shell, *mri_interior; 
  /* creating a shell from the surface */
  printf("computing the shell\n");
  mri_shell = MRIclone(mri_src, NULL);
  mri_shell = MRISshell(mri_src, mris, mri_shell, 1);
  /* creating an interior image from the surface */
  printf("computing an interior image\n");
  mri_interior = MRIclone(mri_src, NULL);
  MRIclear(mri_interior);
  mri_interior = MRISfillInterior(mris, 0.0, mri_interior);
  /* creating the hash table related to the surface vertices */
  printf("computing the hash table\n");
  mht = MHTfillVertexTableRes(mris, NULL, CURRENT_VERTICES, 10);
  /* looping over the nonzero elements of the shell */
  printf("computing the fractions\n");
  volFraction frac;
  octTreeVoxel V; 
  double vox[3], vsize[3];
  vsize[0] = mri_src->xsize; vsize[1] = mri_src->ysize; vsize[2] = mri_src->zsize; 
  for (x = 0; x < width; x++)
    {
      for (y = 0; y < height; y++)
	{
	  for (z = 0; z < depth; z++)
	    {	  
	      if (MRIgetVoxVal (mri_shell, x, y, z, 0) > 125.0)
		{	      
		  /* change of coordinates from image to surface domain */
		  MRIvoxelToSurfaceRAS(mri_shell, x, y, z, &xs, &ys, &zs);
		  /* find the closest vertex to the point */
		  MHTfindClosestVertexGeneric(mht, mris, xs, ys, zs, 10, 2, &v, &vno, &dist);	      
		  /* creating the oct tree voxel structure */
		  vox[0] = xs - vsize[0] / 2.0; 
		  vox[1] = ys - vsize[1] / 2.0; 
		  vox[2] = zs - vsize[2] / 2.0; 
		  V = octTreeVoxelCreate(vox,vsize);
		  /* compute the volume fraction of this voxel */	      
		  frac = MRIcomputeVoxelFractions( V, v, acc, 1, mris);
		  MRIsetVoxVal(mri_fractions,x,y,z,0,frac.frac);
		}
	      else if(MRIgetVoxVal(mri_interior,x,y,z,0) > 0.0)
		MRIsetVoxVal(mri_fractions,x,y,z,0,1.0);

	    }
	}

    }
  return mri_fractions;
}
volFraction MRIcomputeVoxelFractions(octTreeVoxel V, VERTEX *v, double acc, int current_depth, MRI_SURFACE *mris)
{
  /* inputs:
     V: voxel element, which includes center, corners and the vsize
     v: closest vertex to this voxel 
     acc: the accuracy requirement - max error tolerable per voxel 
     current_depth: current depth in the oct-tree structure
     mris: mri surface required for face normal information. 
   */
  int fnum, fiter, j, k; 
  FACE *f;
  double meanNorm[3], meanVert[3];
  /* get the faces the closest vertex to the point (xs,ys,zs) is a part of */ 
  fnum = v->num;   
  /* compute a mean normal based on these faces */ 
  /* underlying assumption is that this mean normal is a good vector 
     to decide if a point is inside or outside a surface
     CAUTION: there might very odd surfaces that can violate this assumption */
  meanNorm[0] = 0.0; meanNorm[1] = 0.0; meanNorm[2] = 0.0;
  for (fiter = 0; fiter < fnum; fiter++)
    {
      f = &mris->faces[v->f[fiter]]; 
      meanNorm[0] = meanNorm[0] + f->nx / (double)fnum;
      meanNorm[1] = meanNorm[1] + f->ny / (double)fnum; 
      meanNorm[2] = meanNorm[2] + f->nz / (double)fnum; 
    }	      
  meanVert[0] = v->x; meanVert[1] = v->y; meanVert[2] = v->z;
  /* detemining if the voxel is completely in or out */
  /* completely in means all the corners and the center of the voxel is inside the surface */  
  int allin = 1;
  int allout = 1;
  double dotProduct;
  for (j = 0; j < 8 ; j++)
    {      
      dotProduct = meanNorm[0] * (V.corn[j][0] - meanVert[0]) + 
	meanNorm[1] * (V.corn[j][1] - meanVert[1]) + meanNorm[2] * (V.corn[j][2] - meanVert[2]);
      if (dotProduct < 0.0) {allout = 0;}
      else {allin = 0;}      
    }

  volFraction frac; frac.frac = 0; frac.err = 0;
  /* relative volume of the voxel to the whole */
  double relativeVolume = 1.0/(double)(pow(2.0,current_depth-1) * pow(2.0,current_depth-1) * pow(2.0,current_depth-1)); 
  if (allin == 1) /* all corners are in return the relative volume of the voxel */
    {            
      frac.frac = relativeVolume;
      frac.err = 0;
    }
  else if (allout == 1) /* all corners are out return 0 */
    {
      frac.frac = 0;
      frac.err = 0;
    }
  else /* some in some out */
    {
      if (relativeVolume < acc) /* the error we are making is small enough */
	{
	  frac.frac = relativeVolume/2;
	  frac.err = relativeVolume/2;
	}
      else /* the error we are making is too big we will redivide and recurse */
	{	  
	  for (k = 0; k < 8 ; k++)
	    {
	      volFraction frac_new = MRIcomputeVoxelFractions( octTreeVoxelDivide(k+1,V), v, acc, current_depth + 1, mris );
	      frac.frac += frac_new.frac;
	      frac.err += frac_new.err;
	    }
	}
    }
  return frac;
}
octTreeVoxel octTreeVoxelCreate (double *vox, double* vsize)
{
  octTreeVoxel v; 
  int k; 
  /*
  for (k = 0; k < 3; k++)
    { v.vox[k] = vox[k]; v.vsize[k] = vsize[k]; }
  */
  for (k = 0; k < 3; k++)
    { v.vox[k] = vox[k]; v.vsize[k] = vsize[k]; }
  /*center*/
  v.cent[0] = v.vox[0] + v.vsize[0] / 2.0; v.cent[1] = v.vox[1] + v.vsize[1] / 2.0; v.cent[2] = v.vox[2] + v.vsize[2] / 2.0;
  /* 000 */
  v.corn[0][0] = v.vox[0];            v.corn[0][1] = v.vox[1];            v.corn[0][2] = v.vox[2];
  /* 100 */
  v.corn[1][0] = v.vox[0] + v.vsize[0]; v.corn[1][1] = v.vox[1];            v.corn[1][2] = v.vox[2];
  /* 010 */
  v.corn[2][0] = v.vox[0];            v.corn[2][1] = v.vox[1] + v.vsize[1]; v.corn[2][2] = v.vox[2];
  /* 001 */
  v.corn[3][0] = v.vox[0];            v.corn[3][1] = v.vox[1];            v.corn[3][2] = v.vox[2] + v.vsize[2];
  /* 110 */
  v.corn[4][0] = v.vox[0] + v.vsize[0]; v.corn[4][1] = v.vox[1] + v.vsize[1]; v.corn[4][2] = v.vox[2];
  /* 011 */
  v.corn[5][0] = v.vox[0];            v.corn[5][1] = v.vox[1] + v.vsize[1]; v.corn[5][2] = v.vox[2] + v.vsize[2];
  /* 101 */
  v.corn[6][0] = v.vox[0] + v.vsize[0]; v.corn[6][1] = v.vox[1];            v.corn[6][2] = v.vox[2] + v.vsize[2];
  /* 111 */
  v.corn[7][0] = v.vox[0] + v.vsize[0]; v.corn[7][1] = v.vox[1] + v.vsize[1]; v.corn[7][2] = v.vox[2] + v.vsize[2];   

  return v; 
}
octTreeVoxel octTreeVoxelDivide (int type, octTreeVoxel v)
{  
  double vsize_new[3];
  double vox_new[3];
  vsize_new[0] = v.vsize[0]/2.0; vsize_new[1] = v.vsize[1]/2.0; vsize_new[2] = v.vsize[2]/2.0;
  switch (type)
    {
    case 1: /* 000 */
      vox_new[0] = v.vox[0]; vox_new[1] = v.vox[1]; vox_new[2] = v.vox[2];
      break;
    case 2: /* 100 */
      vox_new[0] = v.vox[0] + vsize_new[0]; vox_new[1] = v.vox[1]; vox_new[2] = v.vox[2];
      break;
    case 3: /* 010 */
      vox_new[0] = v.vox[0]; vox_new[1] = v.vox[1] + vsize_new[1]; vox_new[2] = v.vox[2];
      break;
    case 4: /* 001 */
      vox_new[0] = v.vox[0]; vox_new[1] = v.vox[1]; vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    case 5: /* 110 */
      vox_new[0] = v.vox[0] + vsize_new[0]; vox_new[1] = v.vox[1] + vsize_new[1]; vox_new[2] = v.vox[2];
      break;
    case 6: /* 101 */
      vox_new[0] = v.vox[0] + vsize_new[0]; vox_new[1] = v.vox[1]; vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    case 7: /* 011 */
      vox_new[0] = v.vox[0]; vox_new[1] = v.vox[1] + vsize_new[1]; vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    case 8: /* 111 */
      vox_new[0] = v.vox[0] + vsize_new[0]; vox_new[1] = v.vox[1] + vsize_new[1]; vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    default: break;
    }
  return octTreeVoxelCreate(vox_new, vsize_new);
}


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

static char vcid[] = "$Id: mris_compute_volume_fractions.c,v 1.2 2013/01/22 16:14:14 enderk Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *VolFile=NULL;
char *SurfFile = NULL;
char *OutFile = NULL;
double Accuracy = -1000.0;

int main(int argc, char *argv[]) {
  
  printf("working!\n");
  int nargs;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");

  if (nargs && argc - nargs == 1) exit (0);

  argc -= nargs;

  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  
  parse_commandline(argc, argv);
  
  check_options();
  if (checkoptsonly) return(0);
  
  dump_options(stdout);

  MRI* mri = MRIread(VolFile);
  MRI_SURFACE* mris = MRISread(SurfFile);    
  printf("running the computation...\n");
  MRI* mri_fractions = MRIcomputeVolumeFractionFromSurface(mris, Accuracy, mri, NULL);
  printf("computation is finished...\n");
  MRIwrite(mri_fractions, OutFile);

  return 0;
}

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

    else if (!strcasecmp(option, "--vol")) {
      if (nargc < 1) CMDargNErr(option,1);
      VolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--surf")){
      if (nargc < 1) CMDargNErr(option,1);
      SurfFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--acc")){
      if (nargc < 1) CMDargNErr(option,1);
      Accuracy = atof(pargv[0]); 
      nargsused = 1; 
    }
    else if (!strcasecmp(option, "--out")){
      if (nargc < 1) CMDargNErr(option,1);
      OutFile = pargv[0];
      nargsused = 1;
    }
    else {
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
  printf("   --vol volume_file : volume \n");
  printf("   --surf surface_file: surface\n");
  printf("   --acc accuracy: required accuracy\n");
  printf("   --out out_file: output volume file for the fractions\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
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
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void check_options(void)
\brief Checks command-line options
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void check_options(void) {
  if(VolFile == NULL || SurfFile == NULL || Accuracy < 0 || OutFile == NULL)
    {
      print_usage(); 
      exit(1);
    }
    
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
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"Working Directory: %s\n",cwd);
  fprintf(fp,"cmdline: %s\n",cmdline);
  /*
  fprintf(fp,"sysname:  %s\n",uts.sysname);
  fprintf(fp,"hostname: %s\n",uts.nodename);
  fprintf(fp,"machine:  %s\n",uts.machine);
  fprintf(fp,"user:     %s\n",VERuser());
  */
  return;
}
