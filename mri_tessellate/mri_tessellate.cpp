/**
 * @brief tessellate a volume to create a surface
 *
 * "Cortical Surface-Based Analysis I: Segmentation and Surface
 * Reconstruction", Dale, A.M., Fischl, B., Sereno, M.I.
 * (1999) NeuroImage 9(2):179-194
 */
/*
 * Original Author: Bruce Fischl
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


//
// mri_tessellate.c
//
// How it works.
//
// 1. pick the boundary voxel position using the grey value given in input
//
//              O O O O O O O O            O O O O O O O O
//              O O O O O O O O            O O O O O O O O
//              O O x x x x O O    picks   O O B B B B B O
//              0 O x x x x O O      ->    O O B O O O B O
//              O O x x x x O O            O 0 B O O O B O
//              O O x x x x O O            O O B O O O B O
//              O O O O O O O O            O O B B B B B O
//              O O O O O O O O  imnr = 3  O O O O O O O O
//
//    You see that there is a bias toward larger i, j, imnr.
//    If you take the real surface in floating value voxel-index
//    coordinate system in the following way,
//
//              -0.5  0.5   1.5
//                |  0 |  1  |  ....
//
//    then the surface position defined as the voxel boundary is given by
//
//              v_surf = v_boundary - 1/2.
//
// 2. Even though the RAS coordinates of surface is given by
//
//    MRIvoxelToWorld(mri,
//                    vertex[k].j-0.5,
//                    vertex[k].i - 0.5,
//                    vertex[k].imnr - 0.5,
//                    &x, &y, &z);
//
//    currently we use the surface RAS being given by
//    conformed volume with c_(r,a,s) = 0
//
//    MRIvoxelToSurfaceRAS()
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mri.h"
#include "fio.h"
#include "const.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "MRIio_old.h"
#include "version.h"
#include "tags.h"
#include "matrix.h"
#include "transform.h"
#include "cma.h"
#include "diag.h"
#include "mrisurf.h"


/////////////////////////////////////////////
#define MAXV 10000000
static long MAXVERTICES = MAXV;
static long  MAXFACES = (2*MAXV);

////////////////////////////////////////////////
// gather globals
static int remove_non_hippo_voxels(MRI *mri) ;
static int all_flag = 0 ;
static int hippo_flag = 0 ;
static int type_changed = 0;
// orig->surface RAS is not MRIvoxelToWorld(), but more involved one
int compatibility= 1;
////////////////////////////////////////////////

tface_type *face;
int *face_index_table0;
int *face_index_table1;

tvertex_type *vertex;
int *vertex_index_table;

static int value;

int main(int argc, char *argv[]) ;
static MRI *read_images(char *fpref) ;
static void add_face(MRI *mri, int imnr, int i, int j, int f, int prev_flag) ;
static int add_vertex(MRI *mri, int imnr, int i, int j) ;
static int facep(MRI *mri, int im0, int i0, int j0, int im1, int i1, int j1) ;
static void check_face(MRI *mri, int im0, int i0, int j0,
                       int im1, int i1,int j1,
                       int f, int n, int v_ind, int prev_flag) ;
static void make_surface(MRI *mri) ;
static void write_binary_surface(char *fname, MRI *mri, std::string& cmdline) ;
static int get_option(int argc, char *argv[]) ;
static void usage_exit(int code);

const char *Progname ;
int UseMRIStessellate=0;

int main(int argc, char *argv[])
{
  char ofpref[STRLEN] /*,*data_dir*/;
  int  nargs ;
  MRI *mri = 0;
  int xnum, ynum, numimg;

  std::string cmdline = getAllInfo(argc, argv, "mri_tessellate");

  nargs = handleVersionOption(argc, argv, "mri_tessellate");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;
  Progname = argv[0] ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  // we removed the option
  if (argc<3) usage_exit(1);
  if (argc == 3)   // value not specified on cmdline - figure it out from hemi in output surf
  {
    char *cp, fname[STRLEN]  ;

    sprintf(ofpref,"%s",argv[2]);
    FileNameOnly(ofpref, fname) ;
    cp = strstr(fname, "h.") ;
    if (cp == NULL)
    {
      ErrorExit(ERROR_UNSUPPORTED,
                "%s: if no fillval is specified, then output hemi must be",
                Progname) ;
    }
    if (*(cp-1) == 'r')
    {
      value = MRI_RIGHT_HEMISPHERE ;
    }
    else if (*(cp-1) == 'l')
    {
      value = MRI_LEFT_HEMISPHERE ;
    }
    else
    {
      ErrorExit(ERROR_UNSUPPORTED,
                "%s: if no fillval is specified, then output hemi must be",
                Progname) ;
    }
    printf("input fill value assumed to be %d based on output name\n", value) ;
  }
  else  // fill value specified explicitly
  {
    sscanf(argv[2],"%d",&value);  // this assumes that argv[2]
    // can be changed to int
    sprintf(ofpref,"%s",argv[3]); // this assumes argv[3] is the file
  }

  // passing dir/COR-
  mri = read_images(argv[1]);
  if (hippo_flag)
  {
    remove_non_hippo_voxels(mri) ;
    all_flag = 1 ;
  }

  printf("%s\n",getVersion().c_str());
  printf("  %s\n",getVersion().c_str());
  fflush(stdout);

  if(UseMRIStessellate){
    MRIS *mris;
    int err;
    printf("Using MRIStessellate())\n");
    mris = MRIStessellate(mri,value,all_flag);
    if(mris == NULL){
      printf("ERROR: mri_tessellate\n");
      exit(1);
    }
    err = MRISwrite(mris,ofpref);
    if(err) exit(1);
    printf("mri_tessellate done\n");
    exit(0);
  }

  // 4 connected (6 in 3D) neighbors
  xnum = mri->width;
  ynum = mri->height;
  numimg = mri->depth;

  face = (tface_type *)lcalloc(MAXFACES,sizeof(tface_type));

  face_index_table0 = (int *)lcalloc(600*ynum*xnum,sizeof(int));
  face_index_table1 = (int *)lcalloc(600*ynum*xnum,sizeof(int));

  vertex = (tvertex_type *)lcalloc(MAXVERTICES,sizeof(tvertex_type));
  vertex_index_table = (int *)lcalloc(800*ynum*xnum,sizeof(int));

  make_surface(mri);

  write_binary_surface(ofpref, mri, cmdline);

  exit(0) ;
}


static MRI *read_images(char *fpref)
{
  MRI *mri = 0;

  mri = MRIread(fpref);
  if (!mri)
  {
    ErrorExit(ERROR_NOFILE, "could not open %s\n", fpref) ;
  }
  // change to UCHAR if not
  if (mri->type!=MRI_UCHAR)
  {
    MRI *mri_tmp ;

    type_changed = 1;
    printf("changing type of input volume to 8 bits/voxel...\n") ;
    mri_tmp = MRIchangeType(mri, MRI_UCHAR, 0.0, 0.999, TRUE) ;
    MRIfree(&mri) ;
    mri = mri_tmp ;
  }
  else
  {
    type_changed = 0 ;
  }

  return mri;
}

static int face_index, vertex_index;


static void
add_face(MRI *mri, int imnr, int i, int j, int f, int prev_flag)
{
  int xnum, ynum, pack;

  xnum = mri->width;
  ynum = mri->height;
  pack = f*ynum*xnum+i*xnum+j;

  if (face_index >= MAXFACES-1)
    ErrorExit(ERROR_NOMEMORY, "%s: max faces %d exceeded",
              Progname,MAXFACES) ;
  if (prev_flag)
  {
    face_index_table0[pack] = face_index;
  }
  else
  {
    face_index_table1[pack] = face_index;
  }
  face[face_index].imnr = imnr; // z
  face[face_index].i = i;       // y
  face[face_index].j = j;       // x
  face[face_index].f = f;
  face[face_index].num = 0;
  face_index++;
}


static int
add_vertex(MRI *mri, int imnr, int i, int j)
{
  int xnum = mri->width;

  int pack = i*(xnum+1)+j;

  if (vertex_index >= MAXVERTICES-1)
    ErrorExit(ERROR_NOMEMORY, "%s: max vertices %d exceeded",
              Progname,MAXVERTICES) ;
  vertex_index_table[pack] = vertex_index;
  vertex[vertex_index].imnr = imnr; // z
  vertex[vertex_index].i = i;       // y
  vertex[vertex_index].j = j;       // x
  vertex[vertex_index].num = 0;
  return vertex_index++;
}


static int
facep(MRI *mri, int im0, int i0, int j0, int im1, int i1, int j1)
{
  int numimg, imax, imin, jmax, jmin;
  numimg = mri->depth;
  // it is so confusing this guy uses j for width and i for height
  jmax = mri->width;
  jmin = 0;
  imax = mri->height;
  imin = 0;
  return (im0>=0&&im0<numimg&&i0>=imin&&i0<imax&&j0>=jmin&&j0<jmax&&
          im1>=0&&im1<numimg&&i1>=imin&&i1<imax&&j1>=jmin&&j1<jmax&&
          MRIvox(mri, j0, i0, im0) != MRIvox(mri, j1, i1, im1) &&
          ((MRIvox(mri, j0, i0, im0)==value ||
            MRIvox(mri, j1, i1, im1)==value) || all_flag));
}

static void
check_face(MRI *mri, int im0, int i0, int j0, int im1, int i1,int j1,
           int f, int n, int v_ind, int prev_flag)
{
  int xnum, ynum, numimg, f_pack, f_ind;
  int imax, imin, jmax, jmin;
  xnum = mri->width;
  ynum = mri->height;
  numimg = mri->depth;
  f_pack = f*ynum*xnum+i0*xnum+j0; // f= 0, 1, 2, 3, 4, 5

  jmax = mri->width;
  jmin = 0;
  imax = mri->height;
  imin = 0;

  if ((im0>=0&&im0<numimg&&i0>=imin&&i0<imax&&j0>=jmin&&j0<jmax&&
       im1>=0&&im1<numimg&&i1>=imin&&i1<imax&&j1>=jmin&&j1<jmax))
  {
    if ((all_flag &&
         ((MRIvox(mri, j0, i0, im0) != 
           MRIvox(mri, j1, i1, im1))/* && (MRIvox(mri, j1, i1, im1) == 0)*/))
        ||
        (((MRIvox(mri, j0, i0, im0)==value) && 
          (MRIvox(mri, j1, i1, im1)!=value))))
    {
      if (n==0)
      {
        add_face(mri, im0,i0,j0,f,prev_flag);
      }
      if (prev_flag)
      {
        f_ind = face_index_table0[f_pack];
      }
      else
      {
        f_ind = face_index_table1[f_pack];
      }
      if (f_ind == Gdiag_no || f_ind == 311366)
      {
        DiagBreak() ;
      }
      face[f_ind].v[n] = v_ind;
      if (vertex[v_ind].num<9)
      {
        vertex[v_ind].f[vertex[v_ind].num++] = f_ind;
      }
    }
  }
}

static void make_surface(MRI *mri)
{
  int imnr,i,j,f_pack,v_ind,f;
  int xnum, ynum, numimg;

  face_index = 0;
  vertex_index = 0;

  xnum = mri->width;
  ynum = mri->height;
  numimg = mri->depth;

  for (imnr=0; imnr<=numimg; imnr++)
  {
    if ((vertex_index || face_index) && !(imnr % 10))
      printf("slice %d: %d vertices, %d faces\n",
             imnr,vertex_index,face_index);
    if (imnr == Gdiag_no)
    {
      DiagBreak() ;
    }
    // i is for width
    for (i=0; i<=ynum; i++)
      for (j=0; j<=xnum; j++)
      {
        if (j == Gx && i == Gy && imnr == Gz)
        {
          DiagBreak() ;
        }
        //              z, y,  x,     z,   y,   x
        if (facep(mri,imnr,i-1,j-1,imnr-1,i-1,j-1) ||
            facep(mri,imnr,i-1,j,imnr-1,i-1,j) ||
            facep(mri,imnr,i,j,imnr-1,i,j) ||
            facep(mri,imnr,i,j-1,imnr-1,i,j-1) ||
            facep(mri,imnr-1,i,j-1,imnr-1,i-1,j-1) ||
            facep(mri,imnr-1,i,j,imnr-1,i-1,j) ||
            facep(mri,imnr,i,j,imnr,i-1,j) ||
            facep(mri,imnr,i,j-1,imnr,i-1,j-1) ||
            facep(mri,imnr-1,i-1,j,imnr-1,i-1,j-1) ||
            facep(mri,imnr-1,i,j,imnr-1,i,j-1) ||
            facep(mri,imnr,i,j,imnr,i,j-1) ||
            facep(mri,imnr,i-1,j,imnr,i-1,j-1))
        {
          v_ind = add_vertex(mri, imnr,i,j);
          check_face(mri, imnr  ,i-1,j-1,imnr-1,i-1,j-1,0,2,v_ind,0);
          check_face(mri, imnr  ,i-1,j  ,imnr-1,i-1,j  ,0,3,v_ind,0);
          check_face(mri, imnr  ,i  ,j  ,imnr-1,i  ,j  ,0,0,v_ind,0);
          check_face(mri, imnr  ,i  ,j-1,imnr-1,i  ,j-1,0,1,v_ind,0);
          check_face(mri, imnr-1,i  ,j-1,imnr-1,i-1,j-1,2,2,v_ind,1);
          check_face(mri, imnr-1,i  ,j  ,imnr-1,i-1,j  ,2,1,v_ind,1);
          check_face(mri, imnr  ,i  ,j  ,imnr  ,i-1,j  ,2,0,v_ind,0);
          check_face(mri, imnr  ,i  ,j-1,imnr  ,i-1,j-1,2,3,v_ind,0);
          check_face(mri, imnr-1,i-1,j  ,imnr-1,i-1,j-1,4,2,v_ind,1);
          check_face(mri, imnr-1,i  ,j  ,imnr-1,i  ,j-1,4,3,v_ind,1);
          check_face(mri, imnr  ,i  ,j  ,imnr  ,i  ,j-1,4,0,v_ind,0);
          check_face(mri, imnr  ,i-1,j  ,imnr  ,i-1,j-1,4,1,v_ind,0);

          check_face(mri, imnr-1,i-1,j-1,imnr  ,i-1,j-1,1,2,v_ind,1);
          check_face(mri, imnr-1,i-1,j  ,imnr  ,i-1,j  ,1,1,v_ind,1);
          check_face(mri, imnr-1,i  ,j  ,imnr  ,i  ,j  ,1,0,v_ind,1);
          check_face(mri, imnr-1,i  ,j-1,imnr  ,i  ,j-1,1,3,v_ind,1);
          check_face(mri, imnr-1,i-1,j-1,imnr-1,i  ,j-1,3,2,v_ind,1);
          check_face(mri, imnr-1,i-1,j  ,imnr-1,i  ,j  ,3,3,v_ind,1);
          check_face(mri, imnr  ,i-1,j  ,imnr  ,i  ,j  ,3,0,v_ind,0);
          check_face(mri, imnr  ,i-1,j-1,imnr  ,i  ,j-1,3,1,v_ind,0);
          check_face(mri, imnr-1,i-1,j-1,imnr-1,i-1,j  ,5,2,v_ind,1);
          check_face(mri, imnr-1,i  ,j-1,imnr-1,i  ,j  ,5,1,v_ind,1);
          check_face(mri, imnr  ,i  ,j-1,imnr  ,i  ,j  ,5,0,v_ind,0);
          check_face(mri, imnr  ,i-1,j-1,imnr  ,i-1,j  ,5,3,v_ind,0);
        }
      }
    for (i=0; i<ynum; i++)
      for (j=0; j<xnum; j++)
        for (f=0; f<6; f++)
        {
          f_pack = f*ynum*xnum+i*xnum+j;
          face_index_table0[f_pack] = face_index_table1[f_pack];
        }
  }
}

#define V4_LOAD(v, x, y, z, r)  (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z, VECTOR_ELT(v,4)=r) ;

static void write_binary_surface(char *fname, MRI *mri, std::string& cmdline)
{
  int k,n;
  double x,y,z;
  FILE *fp;
  MATRIX *m;
  VECTOR *vw, *vv;
  VOL_GEOM vg;

  int useRealRAS = compatibility ? 0 : 1;
  if (useRealRAS == 1)
  {
    printf("using real RAS to save vertex points...\n");
  }
  else
  {
    printf("using the conformed surface RAS to save vertex points...\n");
  }

  fp = fopen(fname,"w");
  if (fp == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not write a file %s",
              Progname,fname) ;
  else
  {
    fprintf(stdout, "writing %s\n", fname);
  }

  fwrite3(-3,fp); //fwrite3(-2,fp); -2 for MRIS_TRIANGULAR_SURFACE, but breaks
  fwrite3(vertex_index,fp);

  fwrite3(face_index,fp);

  // matrix is the same all the time so cache it
  if (useRealRAS==1)
  {
    m = extract_i_to_r(mri);
  }
  else
  {
    m = surfaceRASFromVoxel_(mri);
  }
  printf("using vox2ras matrix:\n") ;
  MatrixPrint(stdout, m) ;

  vv = VectorAlloc(4, MATRIX_REAL);
  vw = VectorAlloc(4, MATRIX_REAL);
  for (k=0; k<vertex_index; k++)
  {

    V4_LOAD(vv, vertex[k].j-0.5, vertex[k].i-0.5, vertex[k].imnr-0.5, 1);
    MatrixMultiply(m, vv, vw);
    // we are doing the same thing as the following, but we save time in
    // calculating the matrix at every point
    // if (useRealRAS == 1)  // use the physical RAS as the vertex point
    //   MRIvoxelToWorld(mri,
    //                   vertex[k].j-0.5,
    //                   vertex[k].i-0.5,
    //                   vertex[k].imnr-0.5,
    //                   &x, &y, &z);
    // else
    //   MRIvoxelToSurfaceRAS(mri,
    //                        vertex[k].j-0.5,
    //                        vertex[k].i-0.5,
    //                        vertex[k].imnr-0.5,
    //                        &x, &y, &z);
    x = V3_X(vw);
    y = V3_Y(vw);
    z = V3_Z(vw);
    fwriteFloat(x,fp);
    fwriteFloat(y,fp);
    fwriteFloat(z,fp);
  }
  MatrixFree(&m);
  VectorFree(&vv);
  VectorFree(&vw);

  for (k=0; k<face_index; k++)
  {
    for (n=0; n<4; n++)
    {
      fwrite3(face[k].v[n],fp);
    }
  }
  // record whether use the physical RAS or not
#if 0
  fwriteInt(TAG_USEREALRAS, fp); // first tag
  fwriteInt(useRealRAS, fp);     // its value
  // save geometry
  fwriteInt(TAG_SURF_GEOM, fp);
  getVolGeom(mri, &vg);
  writeVolGeom(fp, &vg);
#else
  TAGwrite(fp, TAG_USEREALRAS, &useRealRAS, sizeof(useRealRAS)) ;
  if (useRealRAS == 0) // messes up scuba if realras is true - don't know why (BRF)
  {
    fwriteInt(TAG_OLD_SURF_GEOM, fp);
    getVolGeom(mri, &vg);
    writeVolGeom(fp, &vg);
  }
  TAGwrite(fp, TAG_CMDLINE, &cmdline[0], cmdline.size() + 1) ;
#endif

  fclose(fp);
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "an option")) {}
  else if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    usage_exit(0);
  }
  else if (!strcasecmp(option, "-version"))
  {
    fprintf(stderr, "Version: ##version##\n");
    exit(0);
  }
  else if (!stricmp(option, "seed"))
  {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number generator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "maxv") || !stricmp(option, "max_vertices"))
  {
    MAXVERTICES = atol(argv[2]) ;
    MAXFACES = 2*MAXVERTICES ;
    fprintf(stderr,"setting max vertices = %ld, and max faces = %ld\n", MAXVERTICES, MAXFACES);
    nargs = 2 ;
  }
  else if (!stricmp(option, "new")) UseMRIStessellate=1;
  else switch (toupper(*option))
    {
    case 'H':
      hippo_flag = 1 ;
      compatibility = 0;
      printf
      ("tessellating the surface of all hippocampal voxels with different labels\n");
      break ;
    case 'A':
      all_flag = 1 ;
      printf
      ("tessellating the surface of all voxels with different labels\n");
      break ;
    case 'N':
      compatibility = 0;
      printf
      ("surface saved uses coordinates "
       "in the real RAS where c_(r,a,s) != 0\n");
      break;
    case '?':
    case 'U':
      usage_exit(1);
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
int
is_hippo(int label)
{
  switch (label)
  {
  case left_fimbria:
  case left_subiculum:
  case left_CA2_3:
  case left_CA1:
  case left_CA4_DG:
  case left_presubiculum:
  case left_hippocampal_fissure:
  case left_alveus:
  case right_alveus:
  case right_hippocampal_fissure:
  case right_presubiculum:
  case Left_Hippocampus:
  case Right_Hippocampus:
  case right_fimbria:
  case right_subiculum:
  case right_CA2_3:
  case right_CA1:
  case right_CA4_DG:
    return(1) ;
  default:
    break ;
  }
  return(0) ;
}
static int
remove_non_hippo_voxels(MRI *mri)
{
  int    labels[MAX_LABEL+1], x, y, z, l ;

  memset(labels, 0, sizeof(labels)) ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        l = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
        if (is_hippo(l) == 0)
        {
          labels[l] = 1 ;
        }
      }

  for (l = 0 ; l  <= MAX_LABEL ; l++)
    if (labels[l])
    {
      MRIreplaceValues(mri, mri, l, 0) ;
    }
  return(NO_ERROR) ;
}

#include "mri_tessellate.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_tessellate_help_xml,
                mri_tessellate_help_xml_len);
  exit(code) ;
}
