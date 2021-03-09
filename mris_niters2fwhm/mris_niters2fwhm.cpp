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


/*
BEGINHELP

ENDHELP
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "fmriutils.h"
#include "mrisurf.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "annotation.h"
#include "cmdargs.h"
#include "timer.h"
#include "matfile.h"
#include "randomfields.h"
#include "icosahedron.h"
double MRISmeanInterVertexDist(MRIS *surf);

MRI *MRISgaussianSmooth2(MRIS *Surf, MRI *Src, double GStd, MRI *Targ,
                         double TruncFactor);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *subject=NULL, *hemi=NULL, *SUBJECTS_DIR=NULL;
const char *surfname="white";
char *surfpath=NULL;
char tmpstr[2000];

MRIS *surf;
int dof = 100;
int nitersmax = 100;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, nthiter;
  MRI *mri, *var, *mri0, *delta, *deltasm, *xyz;
  double gmax, vrfmn, vrfstd, gstd, fwhm;

  nthiter = 0;
  mri = var = mri0 = delta = deltasm = xyz = NULL;

  nargs = handleVersionOption(argc, argv, "mris_niters2fwhm");
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

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
  surfpath = strcpyalloc(tmpstr);

  if (debug) dump_options(stdout);

  surf = ReadIcoByOrder(7,50);
  //surf = MRISread(surfpath);
  if (surf == NULL) {
    printf("ERROR: could not read %s\n",surfpath);
    exit(1);
  }
  MRIScomputeMetricProperties(surf) ;

  printf("dof %d\n",dof);
  printf("Number of vertices %d\n",surf->nvertices);
  printf("Number of faces    %d\n",surf->nfaces);
  printf("Avg IterVertex     %lf\n",MRISmeanInterVertexDist(surf));

  //----------------------------------------------------------
  // Smooth a delta function and get results
  xyz = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,4);
  MRIcopyMRIS(xyz,surf,0,"x");
  MRIcopyMRIS(xyz,surf,1,"y");
  MRIcopyMRIS(xyz,surf,2,"z");
  MRIcopyMRIS(xyz,surf,3,"area");
  MRIwrite(xyz,"xyz.mgh");

  delta = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);
  MRIsetVoxVal(delta,(int)(surf->nvertices/2),0,0,0,1);
  MRIwrite(delta,"delta.mgh");

  deltasm = MRISgaussianSmooth2(surf, delta, 2, NULL, 5.0);
  //deltasm = MRISsmoothMRI(surf, delta, 2, NULL, deltasm);
  MRIwrite(deltasm,"deltasm.mgh");
  //----------------------------------------------------------
  printf("\n\n");

  mri0 = MRIrandn(surf->nvertices,1,1,dof,0, 1, NULL);
  mri = MRIcopy(mri0,NULL);

  for (nthiter = 2; nthiter <= nitersmax; nthiter++) {
    //MRISsmoothMRI(surf, mri, 1, NULL, mri);
    MRISgaussianSmooth2(surf, mri0, nthiter, mri, 5.0);

    //var = fMRIvariance(mri, dof, 0, var);
    var = fMRIcovariance(mri, 0, mri->nframes-dof, 0, var);
    RFglobalStats(var, NULL, &vrfmn, &vrfstd, &gmax);
    gstd = 1/(2*sqrt(vrfmn*PI));
    fwhm = gstd*sqrt(log(256.0));
    printf("%3d %lf  %lf  %lf %lf\n",nthiter,vrfmn,vrfstd,gstd,fwhm);
    exit(1);
  }

  return 0;
}
/* --------------------------------------------- */
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

    else if (!strcasecmp(option, "--s")) {
      if (nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--h")) {
      if (nargc < 1) CMDargNErr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      surfname = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--dof")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&dof);
      nargsused = 1;
    } else if (!strcasecmp(option, "--niters")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nitersmax);
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
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s --s --h --surf --dof --niters \n",Progname) ;
  printf("\n");
  printf("   --s subject \n");
  printf("   --h hemi \n");
  printf("   --surf surf\n");
  printf("   --dof  dof\n");
  printf("   --niters nitersmax\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"%s\n",Progname);
  fprintf(fp,"FREESURFER_HOME %s\n",getenv("FREESURFER_HOME"));
  fprintf(fp,"SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  fprintf(fp,"cwd       %s\n",cwd);
  fprintf(fp,"cmdline   %s\n",cmdline);
  fprintf(fp,"timestamp %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname   %s\n",uts.sysname);
  fprintf(fp,"hostname  %s\n",uts.nodename);
  fprintf(fp,"machine   %s\n",uts.machine);
  fprintf(fp,"user      %s\n",VERuser());
  if (subject) fprintf(fp,"subject %s\n",subject);
  if (hemi)     fprintf(fp,"hemi     %s\n",hemi);
  if (surfname) fprintf(fp,"surfname %s\n",surfname);
  fprintf(fp,"dof %d\n",dof);
  fprintf(fp,"nitersmax %d\n",nitersmax);

  return;
}
/*---------------------------------------------------------------------*/
double MRISmeanInterVertexDist(MRIS *surf) {
  int vtx, nbrvtx, nnbrs, nthnbr;
  double dx, dy, dz, x0, y0, z0, xn, yn, zn, d;
  double dnbrsum, dnbrmn, dsum;

  dsum = 0.0;
  for (vtx = 0; vtx < surf->nvertices; vtx++) {
    nnbrs = surf->vertices_topology[vtx].vnum;
    x0 = surf->vertices[vtx].x;
    y0 = surf->vertices[vtx].y;
    z0 = surf->vertices[vtx].z;
    dnbrsum=0.0;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      nbrvtx = surf->vertices_topology[vtx].v[nthnbr];
      xn = surf->vertices[nbrvtx].x;
      yn = surf->vertices[nbrvtx].y;
      zn = surf->vertices[nbrvtx].z;
      dx = x0-xn;
      dy = y0-yn;
      dz = z0-zn;
      d = sqrt(dx*dx + dy*dy + dz*dz);
      dnbrsum += d;
    }/* end loop over neighbor */
    dnbrmn = dnbrsum/nnbrs;
    dsum += dnbrmn;
  } /* end loop over vertex */

  return(dsum/surf->nvertices);
}
/*-------------------------------------------------------------------
  MRISgaussianSmooth() - perform gaussian smoothing on a spherical
  surface. The gaussian is defined by stddev GStd and is truncated
  at TruncFactor stddevs. Note: this will change the val2bak of all
  the vertices. See also MRISspatialFilter() and MRISgaussianWeights().
  -------------------------------------------------------------------*/
MRI *MRISgaussianSmooth2(MRIS *Surf, MRI *Src, double GStd, MRI *Targ,
                         double TruncFactor) {
  int vtxno1, vtxno2;
  float val;
  MRI *SrcTmp, *GSum, *GSum2, *nXNbrsMRI, *AreaSum;
  VERTEX *vtx1;
  double Radius, Radius2, dmax, GVar2, f, d, costheta, theta, g, dotprod, ga;
  int n, err, nXNbrs, *XNbrVtxNo, frame;
  double *XNbrDotProd, DotProdThresh;
  double InterVertexDistAvg,InterVertexDistStdDev;
  double VertexRadiusAvg,VertexRadiusStdDev;

  if (Surf->nvertices != Src->width) {
    printf("ERROR: MRISgaussianSmooth: Surf/Src dimension mismatch\n");
    return(NULL);
  }

  if (Targ == NULL) {
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth,
                            MRI_FLOAT, Src->nframes);
    if (Targ==NULL) {
      printf("ERROR: MRISgaussianSmooth: could not alloc\n");
      return(NULL);
    }
  } else {
    if (Src->width   != Targ->width  ||
        Src->height  != Targ->height ||
        Src->depth   != Targ->depth  ||
        Src->nframes != Targ->nframes) {
      printf("ERROR: MRISgaussianSmooth: output dimension mismatch\n");
      return(NULL);
    }
    if (Targ->type != MRI_FLOAT) {
      printf("ERROR: MRISgaussianSmooth: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  /* Make a copy in case it's done in place */
  SrcTmp = MRIcopy(Src,NULL);

  AreaSum = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);//dng

  /* This is for normalizing */
  GSum = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if (GSum==NULL) {
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum\n");
    return(NULL);
  }

  GSum2 = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if (GSum2==NULL) {
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum2\n");
    return(NULL);
  }
  MRIScomputeMetricProperties(Surf);

  nXNbrsMRI = MRIallocSequence(Src->width,Src->height,Src->depth,MRI_FLOAT, 1);

  vtx1 = &Surf->vertices[0] ;
  Radius2 = (vtx1->x * vtx1->x) + (vtx1->y * vtx1->y) + (vtx1->z * vtx1->z);
  Radius  = sqrt(Radius2);
  dmax = TruncFactor*GStd; // truncate after TruncFactor stddevs
  GVar2 = 2*(GStd*GStd);
  f = pow(1/(sqrt(2*M_PI)*GStd),2.0); // squared for 2D
  DotProdThresh = Radius2*cos(dmax/Radius)*(1.0001);

  printf("Radius = %g, gstd = %g, dmax = %g, GVar2 = %g, f = %g, dpt = %g\n",
         Radius,GStd,dmax,GVar2,f,DotProdThresh);

  InterVertexDistAvg    = Surf->avg_vertex_dist;
  InterVertexDistStdDev = Surf->std_vertex_dist;
  VertexRadiusAvg = MRISavgVetexRadius(Surf, &VertexRadiusStdDev);

  printf("Total Area = %g \n",Surf->total_area);
  printf("Dist   = %g +/- %g\n",InterVertexDistAvg,InterVertexDistStdDev);
  printf("Radius = %g +/- %g\n",VertexRadiusAvg,VertexRadiusStdDev);
  printf("nvertices = %d\n",Surf->nvertices);

  /* Initialize */
  for (vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++) {
    MRIFseq_vox(AreaSum,vtxno1,0,0,0)  = 0; //dng
    MRIFseq_vox(GSum,vtxno1,0,0,0)  = 0;
    MRIFseq_vox(GSum2,vtxno1,0,0,0) = 0;
    for (frame = 0; frame < Targ->nframes; frame ++)
      MRIFseq_vox(Targ,vtxno1,0,0,frame) = 0;
    Surf->vertices[vtxno1].val2bak = -1;
  }

  /* These are needed by MRISextendedNeighbors()*/
  XNbrVtxNo   = (int *) calloc(Surf->nvertices,sizeof(int));
  XNbrDotProd = (double *) calloc(Surf->nvertices,sizeof(double));

  if (0) {
    // This will mess up future searches because it sets
    // val2bak to 0
    printf("Starting Search\n");
    err = MRISextendedNeighbors(Surf,0,0,DotProdThresh, XNbrVtxNo,
                                XNbrDotProd, &nXNbrs, Surf->nvertices,1);
    printf("Found %d (err=%d)\n",nXNbrs,err);
    for (n = 0; n < nXNbrs; n++) {
      printf("%d %d %g\n",n,XNbrVtxNo[n],XNbrDotProd[n]);
    }
  }

  // --------------- Loop over target voxel -------------------
  for (vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++) {
    nXNbrs = 0;
    err = MRISextendedNeighbors(Surf,vtxno1,vtxno1,DotProdThresh, XNbrVtxNo,
                                XNbrDotProd, &nXNbrs, Surf->nvertices,1);
    MRIFseq_vox(nXNbrsMRI,vtxno1,0,0,0) = nXNbrs;
    if (vtxno1%10000==0 && Gdiag_no > 0) {
      printf("vtxno1 = %d, nXNbrs = %d\n",vtxno1,nXNbrs);
      fflush(stdout);
    }

    // ------- Loop over neighbors of target voxel --------------
    for (n = 0; n < nXNbrs; n++) {
      vtxno2  = XNbrVtxNo[n];
      dotprod =  XNbrDotProd[n];
      costheta = dotprod/Radius2;

      // cos theta might be slightly > 1 due to precision
      if (costheta > +1.0) costheta = +1.0;
      if (costheta < -1.0) costheta = -1.0;

      // Compute the angle between the vertices
      theta = acos(costheta);

      /* Compute the distance bet vertices along the surface of the sphere */
      d = Radius * theta;

      /* Compute weighting factor for this distance */
      g = f*exp( -(d*d)/(GVar2) );
      ga = g * Surf->vertices[vtxno2].area;

      if (vtxno2 == 81921 && 0) {
        printf("@ %d %d %g %g %g %g %g\n",
               vtxno1,vtxno2,dotprod,costheta,theta,d,g);
      }

      MRIFseq_vox(AreaSum,vtxno1,0,0,0) += Surf->vertices[vtxno2].area;  //dng
      MRIFseq_vox(GSum,vtxno1,0,0,0)  += ga;
      MRIFseq_vox(GSum2,vtxno1,0,0,0) += (ga*ga);

      for (frame = 0; frame < Targ->nframes; frame ++) {
        val = ga*MRIFseq_vox(SrcTmp,vtxno2,0,0,frame);
        MRIFseq_vox(Targ,vtxno1,0,0,frame) += val;
      }

    } /* end loop over vertex2 */

  } /* end loop over vertex1 */

  //MRIwrite(Targ,"ynoscale.mgh");

  /* Normalize */
  if (0) {
    for (vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++) {
      vtx1 = &Surf->vertices[vtxno1] ;
      g = MRIFseq_vox(GSum,vtxno1,0,0,0);
      MRIFseq_vox(GSum2,vtxno1,0,0,0) /= (g*g);
      for (frame = 0; frame < Targ->nframes; frame ++) {
        val = MRIFseq_vox(Targ,vtxno1,0,0,frame);
        MRIFseq_vox(Targ,vtxno1,0,0,frame) = val/g;
        if (vtxno1 == 81921 && 1) {
          printf("%d gsum = %g  src=%g tpre=%g  tpost=%g\n",vtxno1,g,val,
                 MRIFseq_vox(Src,vtxno1,0,0,frame),
                 MRIFseq_vox(Targ,vtxno1,0,0,frame));
        }
      }
    }
  }

  MRIwrite(AreaSum,"areasum.mgh"); //dng
  MRIwrite(GSum,"gsum.mgh");
  MRIwrite(GSum2,"gsum2.mgh");
  MRIwrite(nXNbrsMRI,"nxnbrs.mgh");

  MRIfree(&SrcTmp);
  MRIfree(&GSum);
  MRIfree(&GSum2);
  MRIfree(&nXNbrsMRI);

  free(XNbrVtxNo);
  free(XNbrDotProd);

  return(Targ);
}
