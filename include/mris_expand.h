#ifndef _MRIS_EXPAND_H
#define _MRIS_EXPAND_H

#define MAXVERTICES 10000
#define MAXFACES 10000

#define TRUE 1
#define FALSE 0

#define SQR(x) ((x)*(x))

/* FUNCTION PROTOTYPES */

static void write_geometry(char *fname) ;
static void read_geometry(char *fname);
static void compute_normals(void);
static void normal_face(int f,float *n);
static void expand_geometry(float mm);

/* TYPE DEFINITIONS */

typedef struct ss_vertex_type_
{
  float x,y,z;
  float nx,ny,nz;
  float xb,yb,zb;
  float nxb,nyb,nzb;
  float ox,oy,oz;
  float dx,dy,dz;
  float mx,my,mz;
  float nc,onc,snc;
  float thickness;
  int vnum;
  int v[10];
  int fnum;
  int f[10];
} ss_vertex_type;

#endif
