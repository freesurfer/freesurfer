//
// mri_tessellate.c
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: tosa $
// Revision Date  : $Date: 2003/11/17 20:14:31 $
// Revision       : $Revision: 1.19 $
//
//
// How it works.
//
// 1. pick the boundary voxel position using the grey value given in input 
//
//              O O O O O O O O            O O O O O O O O
//	        O O O O O O O O            O O O O O O O O
//              O O x x x x O O    picks   O O B B B B B O
//	        0 O x x x x O O      ->    O O B O O O B O
//              O O x x x x O O            O 0 B O O O B O  
//              O O x x x x O O            O O B O O O B O
//              O O O O O O O O            O O B B B B B O
//              O O O O O O O O  imnr = 3  O O O O O O O O
//      
//    You see that there is a bias toward larger i, j, imnr.  If you take the real
//    surface in floating value voxel-index coordinate system in the following way,
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
//          MRIvoxelToWorld(mri, vertex[k].j-0.5, vertex[k].i - 0.5, vertex[k].imnr - 0.5,
//                               &x, &y, &z);
//
//    currently we use the surface RAS being given by conformed volume with c_(r,a,s) = 0
//
//          MRIvoxelToSurfaceRAS()
//
char *MRI_TESSELLATE_VERSION = "$Revision: 1.19 $";
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mri.h"
#include "fio.h"
#include "const.h"
#include "diag.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "MRIio_old.h"
#include "version.h"
#include "tags.h"
#include "matrix.h"

#define SQR(x) ((x)*(x))
#define IMGSIZE     256
#define NUMVALS     256
#define MAXIM       256
#define MAXFACES    2000000
#define MAXVERTICES 1000000

static int get_option(int argc, char *argv[]) ;

static int all_flag = 0 ;

// mrisurf.h defines bigger structures (face_type_ and vertex_type_). 
// we don't need big structure here
typedef struct tface_type_
{
  int imnr,i,j,f;
  int num;
  int v[4];
} tface_type;

typedef struct tvertex_type_
{
  int imnr,i,j;
  int num;
  int f[9];
} tvertex_type;

tface_type *face;
int *face_index_table0;
int *face_index_table1;

tvertex_type *vertex;
int *vertex_index_table;

int xnum=256,ynum=256;
unsigned long bufsize;
unsigned char **im[MAXIM];  /* image matrix  */
unsigned char *buf;  /* scratch memory  */
int imnr0,imnr1,numimg;

int imin=0;
int imax=IMGSIZE;
int jmin=0;
int jmax=255;

MRI *mri;
int type_changed = 0;

// orig->surface RAS is not MRIvoxelToWorld(), but more involved one
int compatibility= 1;

static int value;

int main(int argc, char *argv[]) ;
static void read_images(char *fpref) ;
static void add_face(int imnr, int i, int j, int f, int prev_flag) ;
static int add_vertex(int imnr, int i, int j) ;
static int facep(int im0, int i0, int j0, int im1, int i1, int j1) ;
static void check_face(int im0, int i0, int j0, int im1, int i1,int j1, 
                       int f, int n, int v_ind, int prev_flag) ;
static void make_surface(void) ;
static void write_binary_surface2(char *fname) ;
#if 0
static void write_surface(void) ;
static void write_binary_surface(char *fname) ;
#endif

char *Progname ;

int
main(int argc, char *argv[])
{
  char ofpref[STRLEN] /*,*data_dir*/;
  int  nargs ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_tessellate.c,v 1.19 2003/11/17 20:14:31 tosa Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
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
  if (argc<4)
  {
    printf("Usage: %s <option> <input volume> <label-value> <output surface>\n",Progname);
    exit(1);
  }

  sscanf(argv[2],"%d",&value);        // this assumes that argv[2] can be changed to int

  sprintf(ofpref,"%s",argv[3]);      // this assumes argv[3] is the file
  
  face = (tface_type *)lcalloc(MAXFACES,sizeof(tface_type));
  // 4 connected (6 in 3D) neighbors
  face_index_table0 = (int *)lcalloc(6*ynum*xnum,sizeof(int));
  face_index_table1 = (int *)lcalloc(6*ynum*xnum,sizeof(int));
  
  vertex = (tvertex_type *)lcalloc(MAXVERTICES,sizeof(tvertex_type));
  vertex_index_table = (int *)lcalloc(8*ynum*xnum,sizeof(int));
  
  // passing dir/COR- 
  read_images(argv[1]);
  
  make_surface();

  write_binary_surface2(ofpref);

  return 0;
}

static void
read_images(char *fpref)
{
  mri = MRIread(fpref);
  imnr0 = mri->imnr0;
  imnr1 = mri->imnr1;
  xnum = mri->width;
  ynum = mri->height;

  numimg = imnr1-imnr0+1;

  imin = 0;
  imax = ynum;
  jmin = 0;
  jmax = xnum;

  // change to UCHAR if not
  if (mri->type!=MRI_UCHAR)
  {
    MRI *mri_tmp ;
    
    type_changed = 1;
    printf("changing type of input volume to 8 bits/voxel...\n") ;
    mri_tmp = MRIchangeType(mri, MRI_UCHAR, 0.0, 0.999, FALSE) ;
    MRIfree(&mri) ; 
    mri = mri_tmp ;
  }
  else
    type_changed = 0 ;
  

/* Allocate memory */

  bufsize = ((unsigned long)xnum)*ynum;
  buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
}

static int face_index, vertex_index;

static void
add_face(int imnr, int i, int j, int f, int prev_flag)
{
  int pack = f*ynum*xnum+i*xnum+j;

  if (face_index >= MAXFACES-1)
    ErrorExit(ERROR_NOMEMORY, "%s: max faces %d exceeded", 
              Progname,MAXFACES) ;
  if (prev_flag)
    face_index_table0[pack] = face_index;
  else
    face_index_table1[pack] = face_index;
  face[face_index].imnr = imnr; // z
  face[face_index].i = i;       // y
  face[face_index].j = j;       // x
  face[face_index].f = f;
  face[face_index].num = 0;
  face_index++;
}

static int
add_vertex(int imnr, int i, int j)
{
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
facep(int im0, int i0, int j0, int im1, int i1, int j1)
{
  return (im0>=0&&im0<numimg&&i0>=imin&&i0<imax&&j0>=jmin&&j0<jmax&&
          im1>=0&&im1<numimg&&i1>=imin&&i1<imax&&j1>=jmin&&j1<jmax&&
          MRIvox(mri, j0, i0, im0) != MRIvox(mri, j1, i1, im1) &&
          ((MRIvox(mri, j0, i0, im0)==value || MRIvox(mri, j1, i1, im1)==value) || all_flag));
}

static void
check_face(int im0, int i0, int j0, int im1, int i1,int j1, 
	   int f, int n, int v_ind, int prev_flag)
{
  int f_pack = f*ynum*xnum+i0*xnum+j0; // f= 0, 1, 2, 3, 4, 5
  int f_ind;

  if ((im0>=0&&im0<numimg&&i0>=imin&&i0<imax&&j0>=jmin&&j0<jmax&&
       im1>=0&&im1<numimg&&i1>=imin&&i1<imax&&j1>=jmin&&j1<jmax))
  {
    if ((all_flag && ((MRIvox(mri, j0, i0, im0) != MRIvox(mri, j1, i1, im1)!=value))) ||
        (((MRIvox(mri, j0, i0, im0)==value) && (MRIvox(mri, j1, i1, im1)!=value))))
    {
      if (n==0)
      {
        add_face(im0,i0,j0,f,prev_flag);
      }
      if (prev_flag)
        f_ind = face_index_table0[f_pack];
      else
        f_ind = face_index_table1[f_pack];
      face[f_ind].v[n] = v_ind;
      if (vertex[v_ind].num<9)
        vertex[v_ind].f[vertex[v_ind].num++] = f_ind;
    }
  }
}

static void
make_surface(void)
{
  int imnr,i,j,f_pack,v_ind,f;

  face_index = 0;
  vertex_index = 0;

  for (imnr=0;imnr<=numimg;imnr++)
  {
    if ((vertex_index || face_index) && !(imnr % 10))
      printf("slice %d: %d vertices, %d faces\n",imnr,vertex_index,face_index);
    for (i=0;i<=ynum;i++)
    for (j=0;j<=xnum;j++)
    {
      if (facep(imnr,i-1,j-1,imnr-1,i-1,j-1) ||
          facep(imnr,i-1,j,imnr-1,i-1,j) ||
          facep(imnr,i,j,imnr-1,i,j) ||
          facep(imnr,i,j-1,imnr-1,i,j-1) ||
          facep(imnr-1,i,j-1,imnr-1,i-1,j-1) ||
          facep(imnr-1,i,j,imnr-1,i-1,j) ||
          facep(imnr,i,j,imnr,i-1,j) ||
          facep(imnr,i,j-1,imnr,i-1,j-1) ||
          facep(imnr-1,i-1,j,imnr-1,i-1,j-1) ||
          facep(imnr-1,i,j,imnr-1,i,j-1) ||
          facep(imnr,i,j,imnr,i,j-1) ||
          facep(imnr,i-1,j,imnr,i-1,j-1))
      {
        v_ind = add_vertex(imnr,i,j);
        check_face(imnr  ,i-1,j-1,imnr-1,i-1,j-1,0,2,v_ind,0);
        check_face(imnr  ,i-1,j  ,imnr-1,i-1,j  ,0,3,v_ind,0);
        check_face(imnr  ,i  ,j  ,imnr-1,i  ,j  ,0,0,v_ind,0);
        check_face(imnr  ,i  ,j-1,imnr-1,i  ,j-1,0,1,v_ind,0);
        check_face(imnr-1,i  ,j-1,imnr-1,i-1,j-1,2,2,v_ind,1);
        check_face(imnr-1,i  ,j  ,imnr-1,i-1,j  ,2,1,v_ind,1);
        check_face(imnr  ,i  ,j  ,imnr  ,i-1,j  ,2,0,v_ind,0);
        check_face(imnr  ,i  ,j-1,imnr  ,i-1,j-1,2,3,v_ind,0);
        check_face(imnr-1,i-1,j  ,imnr-1,i-1,j-1,4,2,v_ind,1);
        check_face(imnr-1,i  ,j  ,imnr-1,i  ,j-1,4,3,v_ind,1);
        check_face(imnr  ,i  ,j  ,imnr  ,i  ,j-1,4,0,v_ind,0);
        check_face(imnr  ,i-1,j  ,imnr  ,i-1,j-1,4,1,v_ind,0);

        check_face(imnr-1,i-1,j-1,imnr  ,i-1,j-1,1,2,v_ind,1);
        check_face(imnr-1,i-1,j  ,imnr  ,i-1,j  ,1,1,v_ind,1);
        check_face(imnr-1,i  ,j  ,imnr  ,i  ,j  ,1,0,v_ind,1);
        check_face(imnr-1,i  ,j-1,imnr  ,i  ,j-1,1,3,v_ind,1);
        check_face(imnr-1,i-1,j-1,imnr-1,i  ,j-1,3,2,v_ind,1);
        check_face(imnr-1,i-1,j  ,imnr-1,i  ,j  ,3,3,v_ind,1);
        check_face(imnr  ,i-1,j  ,imnr  ,i  ,j  ,3,0,v_ind,0);
        check_face(imnr  ,i-1,j-1,imnr  ,i  ,j-1,3,1,v_ind,0);
        check_face(imnr-1,i-1,j-1,imnr-1,i-1,j  ,5,2,v_ind,1);
        check_face(imnr-1,i  ,j-1,imnr-1,i  ,j  ,5,1,v_ind,1);
        check_face(imnr  ,i  ,j-1,imnr  ,i  ,j  ,5,0,v_ind,0);
        check_face(imnr  ,i-1,j-1,imnr  ,i-1,j  ,5,3,v_ind,0);
      }
    }
    for (i=0;i<xnum;i++)
    for (j=0;j<ynum;j++)
    for (f=0;f<6;f++)
    {
      f_pack = f*ynum*xnum+i*xnum+j;
      face_index_table0[f_pack] = face_index_table1[f_pack];
    }
  }
}

#define V4_LOAD(v, x, y, z, r)  (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z, VECTOR_ELT(v,4)=r) ;

static void
write_binary_surface2(char *fname)
{
  int k,n;
  double x,y,z;
  FILE *fp;
  MATRIX *m;
  VECTOR *vw, *vv;

  int useRealRAS = compatibility ? 0 : 1;
  if (useRealRAS == 1)
    printf("using real RAS to save vertex points...\n");
  else
    printf("using the conformed surface RAS to save vertex points...\n");

  fp = fopen(fname,"w");
  if (fp == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not write a file %s", 
              Progname,fname) ;
  else
    fprintf(stdout, "writing %s\n", fname);

  fwrite3(-1,fp);
  fwrite3(vertex_index,fp);
  fwrite3(face_index,fp);

  // matrix is the same all the time so cache it
  if (useRealRAS==1)
    m = extract_i_to_r(mri);
  else
    m = surfaceRASFromVoxel_(mri);

  vv = VectorAlloc(4, MATRIX_REAL);
  vw = VectorAlloc(4, MATRIX_REAL);
  for (k=0;k<vertex_index;k++)
  {

    V4_LOAD(vv, vertex[k].j-0.5, vertex[k].i-0.5, vertex[k].imnr-0.5, 1);
    MatrixMultiply(m, vv, vw);
    // we are doing the same thing as the following, but we save time in
    // calculating the matrix at every point
    // if (useRealRAS == 1)  // use the physical RAS as the vertex point
    //   MRIvoxelToWorld(mri, vertex[k].j-0.5, vertex[k].i-0.5, vertex[k].imnr-0.5, &x, &y, &z);
    // else
    //   MRIvoxelToSurfaceRAS(mri, vertex[k].j-0.5, vertex[k].i-0.5, vertex[k].imnr-0.5, &x, &y, &z);
    x = V3_X(vw);
    y = V3_Y(vw);
    z = V3_Z(vw);

    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
  }
  MatrixFree(&m);
  VectorFree(&vv);
  VectorFree(&vw);

  for (k=0;k<face_index;k++)
  {
    for (n=0;n<4;n++)
      fwrite3(face[k].v[n],fp);
  }
  // record whether use the physical RAS or not
  fwriteInt(TAG_USEREALRAS, fp); // first tag
  fwriteInt(useRealRAS, fp);     // its value
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

  if (!stricmp(option, "an option"))
  {
  }
  else if (!strcasecmp(option, "-version"))
  {
    fprintf(stderr, "Version: %s\n", MRI_TESSELLATE_VERSION);
    exit(0);
  }
  else switch (toupper(*option))
  {
  case 'A':
    all_flag = 1 ;
    printf("tessellating the surface of all voxels with different labels\n") ;
    break ;
  case 'N':
    compatibility = 0;
    printf("surface saved uses coordinates in the real RAS where c_(r,a,s) != 0\n");
    break;
  case '?':
  case 'U':
    printf("Usage: %s <input volume> <label> <output surface>\n",Progname);
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
