/**
 * @brief apply a specified displacement field to warp a surface
 *
 * applies a specified displacement field specified as a volume (M3Z)
 * or as an "overlay" (MGZ or NIfTI), and checks quality of output
 * surface
 */
/*
 * Original Author: jonathan polimeni
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


#include <stdio.h>
#include <stdlib.h>

// nint
#include "version.h"
#include "macros.h"

#include "error.h"
#include "diag.h"
#include "timer.h"

#include "mri.h"
#include "mrisurf.h"
#include "gcamorph.h"

#include "registerio.h"

#include "resample.h"

// string_to_type
#include "mri_identify.h"


const char *Progname = NULL;
static int  parse_commandline(int argc, char **argv);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;

static void argnerr(char *, int) ;

static int mris_warp__check_deformation(MRI **, unsigned short int);
MATRIX* mris_warp__TkReg2vox(MRI_SURFACE *);


static int debug = 0;



static char *hemi    = NULL;
static char *srcsubject = NULL;

static const char *surf_filename = "white";
static const char *deformvol_filename = "deform_vol_abs.nii.gz";
static const char *deformsurf_filename = "deform_surf_abs.mgz";
static const char *m3z_filename = "deform_vol_rel.m3z";
static const char *reg_filename = "register.dat";
static const char *warpsurf_filename = "white.warp";

int FLAG__abs = 0;


//static char *outfile  = NULL;
//static char *outtypestring = NULL;
static int  outtype = MRI_VOLUME_TYPE_UNKNOWN;

char *basename (char* path)
{
  char *ptr = strrchr (path, '/');
  return ptr ? ptr + 1 : (char*)path;
}


static MATRIX *Qsrc, *Fsrc, *Wsrc, *Dsrc, *vox2ras;
static float ProjFrac = 0.0;
//
//static char  *interpmethod_string = "nearest";
static int  interpmethod = -1;
//
static int  float2int_src;
//static char *float2int_string = "round";
static int  float2int = -1;
static MRI *SrcHitVol;
static int   ProjDistFlag = 0;
//
//static MRI *SrcVol, *SurfVals, *SurfVals2, *overlay;
static MRI *overlay;



int main(int argc, char *argv[])
{

  int          nargs = 0;
  int          err;

  MRI_SURFACE  *mris  = NULL;
  MRI          *mri   = NULL;

  int          msec, minutes, seconds ;
  Timer start ;

  float ipr, bpr, intensity;

  MATRIX       *p_vox, *p_ras, *surfaceRAS;

  GCA_MORPH    *gcam1, *gcam2 ;


  // nargs = handleVersionOption(argc, argv, "mris_warp");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = basename(argv[0]) ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  argc--;
  argv++;

  if (argc == 0)
  {
    usage_exit();
  }

  parse_commandline(argc, argv);


  //==--------------------------------------------------------------
  // 0) read in a surface and optionally a registration matrix


  printf("reading registration file %s\n", reg_filename);

  //  float2int = FLT2INT_ROUND;
  float2int = 0;

  // read in M0_vox2ras matrix as Dsrc
  err = regio_read_register(reg_filename, &srcsubject, &ipr, &bpr,
                            &intensity, &Dsrc, &float2int_src);
  if ( err )
  {
    exit(1);
  }

  printf("subject is '%s'\n", srcsubject);


  mris = MRISread(surf_filename) ;
  if (!mris)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", 
              Progname, surf_filename) ;
  }

  // check input surface for intersections
  {
    MRIS_HASH_TABLE  *mht ;
    int fno ;
    int count_input = 0;
    float grid = 1.0;

    mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, grid);

    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      if (MHTdoesFaceIntersect(mht, mris, fno))
      {
        count_input++;
        //printf("%d\n", fno);
      }
    }
    printf("%d intersections found in INPUT surface at a "
           "grid resolution of %2.2f mm\n", count_input, grid);
    MHTfree(&mht) ;
  }

  //==-----------------------------------------------------------
  // 1) read in a volume warp, i.e., an M3Z file, 
  // OR a surface warp, i.e., an N x 1 x 1 x 3 MGZ


  printf("reading %s -- assuming ABSOLUTE position convention\n", 
         deformvol_filename);
  mri = MRIread(deformvol_filename) ;
  if (mri == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not read warp volume %s\n",
              Progname, deformvol_filename);
  }


// NOTE: mri_convert expects M3Z files that are specified in the
//  RELATIVE position convention
//  printf("reading %s\n", m3z_filename);
//  gcam2 = GCAMread(m3z_filename) ;
//  if (gcam2 == NULL) return(1);


  //==---------------------------------------------------------------
  // 2) if it is a volume warp, find values on surface based on
  //    registration similar to mri_vol2surf


  // find voxels intersecting input surface

  SrcHitVol = MRIallocSequence(mri->width,mri->height,
                               mri->depth,MRI_FLOAT,1);

  MRIcopyHeader(mri,SrcHitVol);


  /* Wsrc: Get the source warping Transform */
  Wsrc = NULL;
  /* Fsrc: Get the source FOV registration matrix */
  Fsrc = NULL;
  // Compute vox2ras for source
  vox2ras = MRIxfmCRS2XYZtkreg(mri);
  // Compute ras2vox (Qsrc: the quantization matrix)
  Qsrc = MatrixInverse(vox2ras,NULL);

  vox2ras = MRIxfmCRS2XYZ(mri, 0);

  // TODO: provide user command line option to specify this
  interpmethod = interpolation_code("trilinear");

  // TODO: provide user command line option to specify this
  // (although not sure that they would want to change this)
  ProjFrac = 0.0;
  ProjDistFlag = 0;

  unsigned short int bIsAbsoluteDeform = 1;

  //* flip volume and scale
  mris_warp__check_deformation(&mri, bIsAbsoluteDeform);

  // TODO: provide user command line option to specify this
  outtype = MRI_GZIPPED;
  //  outtype = string_to_type(outtypestring);
  outtype = string_to_type("mgz");

  overlay =
    vol2surf_linear(mri, Qsrc, Fsrc, Wsrc, Dsrc,
                    mris, ProjFrac, interpmethod, float2int, SrcHitVol,
                    ProjDistFlag, 1);

  int nsrchits, c, r, s;

  /* count the number of source voxels hit */
  nsrchits = 0;
  for (c=0; c < SrcHitVol->width; c++)
  {
    for (r=0; r < SrcHitVol->height; r++)
    {
      for (s=0; s < SrcHitVol->depth; s++)
      {
        if (MRIFseq_vox(SrcHitVol,c,r,s,0) > 0.5)
        {
          nsrchits++;
        }
      }
    }
  }
  printf("number of source voxels hit = %d\n",nsrchits);


  // TODO: provide user command line option to specify this
//  err = MRIwriteType(overlay,"overlay.mgz",outtype);
//  if(err){
//    printf("ERROR: saving \n");
//    exit(1);
//  }

  gcam1 = GCAMalloc(mri->width, mri->height, mri->depth) ;
  GCAMinitVolGeom(gcam1, mri, mri) ;

#if 1
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam1, mri) ;
#else
  GCAMreadWarpFromMRI(gcam1, mri) ;
#endif

  GCAMfree(&gcam1);

  if ( 0 )
  {
    GCAMfree(&gcam2);
  }

  // to trace where the gcam is scaled by the voxel size...
#if 0
  GCA_MORPH_NODE  *gcamn ;

  gcamn = &gcam->nodes[0][0][0] ;
  printf("g ox oy oz: %2.2f %2.2f %2.2f\ng x y z: %2.2f %2.2f %2.2f\n"
         "g nx ny nz: %2.2f %2.2f %2.2f\n",
         gcamn->origx, gcamn->origy, gcamn->origz, 
         gcamn->x, gcamn->y, gcamn->z, 
         gcamn->xn, gcamn->yn, gcamn->zn);
#endif

  //==---------------------------------------------------------
  // 3) reposition surface by replacing RAS XYZ coordinates with 
  // those based on warp

  //    question: how to extrapolate if volume too small?

  MATRIX* iM = NULL;

  iM = mris_warp__TkReg2vox(mris);

  p_vox = MatrixAlloc(4, 1, MATRIX_REAL);
  *MATRIX_RELT(p_vox, 4, 1) = 1;

  // TODO: identify vertices that are outside of deformation volume
  // (whose positions will take the value '0') and leave these vertices
  // unchanged.

  int vno;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX* v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    *MATRIX_RELT(p_vox, 1, 1) = MRIgetVoxVal(overlay,vno,0,0,0);
    *MATRIX_RELT(p_vox, 2, 1) = MRIgetVoxVal(overlay,vno,0,0,1);
    *MATRIX_RELT(p_vox, 3, 1) = MRIgetVoxVal(overlay,vno,0,0,2);

    p_ras = MatrixMultiply(vox2ras, p_vox, NULL);

    // do not need register.dat matrix here
    surfaceRAS = MatrixMultiply(iM, p_ras, NULL);

    int f;
    for(f=0; f < 3; f++)
    {
      MRIsetVoxVal(overlay,vno,0,0,f, surfaceRAS->rptr[f+1][1]);
    }

    // for debugging
    if ( 0 && vno == 0 )
    {
      printf("vx: %2.2f, vy: %2.2f, vz: %2.2f; \n"
             "px: %2.2f, py: %2.2f, pz: %2.2f; \n"
             "px: %2.2f, py: %2.2f, pz: %2.2f; \n"
             "qx: %2.2f, qy: %2.2f, qz: %2.2f\n",
             v->x, v->y, v->z,
             p_vox->rptr[1][1], p_vox->rptr[2][1], p_vox->rptr[3][1],
             p_ras->rptr[1][1], p_ras->rptr[2][1], p_ras->rptr[3][1],
             surfaceRAS->rptr[1][1], 
             surfaceRAS->rptr[2][1], 
             surfaceRAS->rptr[3][1]);
    }

    MRISsetXYZ(mris, vno,
      surfaceRAS->rptr[1][1],
      surfaceRAS->rptr[2][1],
      surfaceRAS->rptr[3][1]);
  }

  MatrixFree(&p_vox);
  MatrixFree(&p_ras);
  MatrixFree(&surfaceRAS);


  err = MRIwriteType(overlay,"overlay_warp.mgz",outtype);
  if(err)
  {
    printf("ERROR: saving \n");
    exit(1);
  }

  // TODO: add command line string to surface header (a la HIPS)
  printf("writing to %s\n", warpsurf_filename);
  MRISwrite(mris, warpsurf_filename);


  //==-------------------------------------------------------------------
  // 4) check quality of output surface ( mrisurf.c:IsMRISselfIntersecting )

  MRI *TrgVol = NULL;

  /* allocate a "volume" to hold the output */
  TrgVol = MRIallocSequence(mris->nvertices,1,1,MRI_FLOAT,1);
  if (TrgVol == NULL)
  {
    printf("error\n");
  }
  MRIcopyHeader(mri,TrgVol);
  // Dims here are meaningless, but setting to 1 means "volume" will be
  // number of vertices.
  TrgVol->xsize = 1;
  TrgVol->ysize = 1;
  TrgVol->zsize = 1;

  int count = 0;

  {
    MRIS_HASH_TABLE  *mht ;
    int fno ;
    float grid = 1.0;

    //
    mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, grid);

    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      // wrapper around "mhtDoesFaceVoxelListIntersect", which seems
      // to test whether multiple faces intersect a 1mm "hash bin"

      // TODO: explore alternate methods for detecting
      // intersections, e.g., tri_tri_inter() of tetgen, or

      if (MHTdoesFaceIntersect(mht, mris, fno))
      {
        count++;
        //            printf("%d\n", fno);

        //                                        r s f v
        MRIsetVoxVal(TrgVol,mris->faces[fno].v[0],0,0,0,1);
        MRIsetVoxVal(TrgVol,mris->faces[fno].v[1],0,0,0,1);
        MRIsetVoxVal(TrgVol,mris->faces[fno].v[2],0,0,0,1);
      }
    }

    printf("%d intersections found in WARPED surface at a grid "
           "resolution of %2.2f mm\n", count, grid);
    MHTfree(&mht) ;
  }

  // TODO: output label file instead
//  err = MRIwriteType(TrgVol,"overlay_intersections.mgz",outtype);
//  if(err){
//    printf("ERROR: saving \n");
//    exit(1);
//  }

  if ( count )
  {
    printf("surface intersections found: %d\n", count);
  }

  MRIfree(&mri);
  MRISfree(&mris);

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "%s took %d minutes and %d seconds.\n", 
          Progname, minutes, seconds) ;

  exit(0);
  return(0);
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1)
  {
    usage_exit();
  }

  nargc = argc;
  pargv = argv;
  while (nargc > 0)
  {
    option = pargv[0];
    if (debug)
    {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))
    {
      print_help() ;
    }
    else if (!strcasecmp(option, "--version"))
    {
      print_version() ;
    }
    else if (!strcasecmp(option, "--debug"))
    {
      debug = 1;
    }

    /*--------------*/

    else if (!strcasecmp(option, "--deformvol"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      deformvol_filename = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--m3z"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      m3z_filename = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--deformsurf"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      deformsurf_filename = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--out"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      warpsurf_filename = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--rel"))
    {
      FLAG__abs = 0;
      nargsused = 0;
    }
    else if (!strcasecmp(option, "--abs"))
    {
      FLAG__abs = 1;
      nargsused = 0;
    }
    else if (!strcasecmp(option, "--reg"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      reg_filename = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--regheader"))
    {
      printf("currently unsupported\n");
      exit(1);
      return(1);
    }


    else if (!strcmp(option, "--surf"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      surf_filename = pargv[0];
      if( !strcmp(surf_filename, "inflated") )
      {
        printf("\nWARNING: do you really want to warp the "
               "*inflated* surface?\n\n");
        exit(1);
      }
      nargsused = 1;
    }
    else if ( !strcmp(option, "--hemi") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      hemi = pargv[0];
      nargsused = 1;
    }

    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}


/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1)
  {
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  }
  else
  {
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  }
  exit(-1);
}


/* --------------------------------------------- */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}


/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --deformvol <filename>   volume containing deformation\n");
  printf("   --deformsurf <filename>  (NOT YET IMPLEMENTED)\n");
  printf("   --m3z <filename>         M3Z file containing deformation\n");
  printf("   --reg <filename>         register.dat file between surface and volume\n");
  printf("   --regheader              (NOT YET IMPLEMENTED)\n");
  printf("   --abs                    absolute coordinate displacement convention (default)\n");
  printf("   --rel                    (NOT YET IMPLEMENTED)\n");
  printf("   --surf <filename>        surface file to warp\n");
  printf("   --out <filename>         name for output surface (if does not contain '/'\n");
  printf("                            outputs to same directory as input surface)\n");
  printf("\n");
  printf("   --help        print out information on how to use this program\n");
  printf("   --version     print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;

  printf(
    "This program will warp a surface using a specified deformation field.\n"
  ) ;

  exit(1) ;
}

/* --------------------------------------------- */
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

/* --------------------------------------------- */
static int mris_warp__check_deformation(MRI **pmri, 
                                        unsigned short int bIsAbsoluteDeform)
{
  /*
   * FSL uses a hybrid convention for the deformations produced by
   * their warp tools. the deformation can be in an absolute or
   * relative position convention, but the deformation direction is
   * always along the axes of the volume AND the units of the
   * deformation are mm---so it is a mix of index-based and
   * spatial-based coordinates.
   *
   * this function returns a deformation as absolute positions in
   * index-based coordinates and that
   */

  MATRIX    *m;

  MRI       *mri2 = NULL;
  int       c=0, r=0, s=0;
  float     val, det=0.0;

  MRI       *mri;

  printf("%s\n", __func__);

  mri = *pmri;

  m = MRIgetVoxelToRasXform(mri) ;

  // NOTE: this assumes a standard siemens image orientation in which
  // case a neurological orientation means that the first frame is
  // flipped

  det = MatrixDeterminant(m);

  if ( det > 0 )
  {
    fprintf(stdout, "non-negative Jacobian determinant -- "
            "converting to radiological ordering\n");
  }

  mri2 = MRIcopy(mri,NULL);
  for(c=0; c < mri->width; c++)
  {
    for(r=0; r < mri->height; r++)
    {
      for(s=0; s < mri->depth; s++)
      {

        val = MRIgetVoxVal(mri,c,r,s,0) / mri->xsize;
        if ( !bIsAbsoluteDeform )
        {
          // make relative -- not yet implemented
        }

        // only flip first frame (by subtracting shift from width)
        if ( det > 0 )
        {
          MRIsetVoxVal(mri2,c,r,s,0,mri->width-val-1);
        }
        else
        {
          MRIsetVoxVal(mri2,c,r,s,0,val);
        }

        val = MRIgetVoxVal(mri,c,r,s,1) / mri->ysize;
        if ( !bIsAbsoluteDeform )
        {
          // make relative -- not yet implemented
        }
        MRIsetVoxVal(mri2,c,r,s,1,val);

        val = MRIgetVoxVal(mri,c,r,s,2) / mri->zsize;
        if ( !bIsAbsoluteDeform )
        {
          // make relative -- not yet implemented
        }

        MRIsetVoxVal(mri2,c,r,s,2,val);
      }
    }
  }

  MRIfree(&mri);
  *pmri = mri2;

  MatrixFree(&m) ;

  return(0);
}

/* --------------------------------------------- */
MATRIX* mris_warp__TkReg2vox(MRI_SURFACE *mris)
{

  MATRIX *iM = NULL;

  if ( mris->vg.valid )
  {
    MRI* tmp = MRIallocHeader(mris->vg.width,
                              mris->vg.height, 
                              mris->vg.depth,
                              MRI_UCHAR,
                              1);

    useVolGeomToMRI(&mris->vg, tmp);
    MATRIX* vox2rasScanner = MRIxfmCRS2XYZ(tmp, 0);
    MATRIX* vo2rasTkReg = MRIxfmCRS2XYZtkreg(tmp);
    MATRIX* vox2rasTkReg_inv = MatrixInverse( vo2rasTkReg, NULL );
    MATRIX* M = MatrixMultiply( vox2rasScanner, vox2rasTkReg_inv, NULL );
    iM = MatrixInverse( M, NULL );

    MRIfree( &tmp );
    MatrixFree( &vox2rasScanner );
    MatrixFree( &vo2rasTkReg );
    MatrixFree( &vox2rasTkReg_inv );
    MatrixFree( &M );
  }
  else
  {
    return NULL;
  }

  return iM;
}
