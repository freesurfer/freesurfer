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


#ifndef TOPO_PARMS
#define TOPO_PARMS

typedef struct
{
  int   vno1, vno2 ;
  float len ;
  short used ;
}
EDGE ;

typedef struct
{
  int **faces ;   /* array of pointer to list of ambiguous faces */
  int *nfaces ;   /* array of ints specifying how many faces are ambigous */
}
FACE_DEFECT_LIST, FDL ;

typedef struct
{
  float  nx, ny, nz ;    /* average normal in the defect */
  float  cx, cy, cz ;    /* centroid of the defect */
  float  area ;          /* total (original) area of the defect */
  int    *vertices ;     /* vertices in the defect */
  char   *status ;       /* keep or discard */
  int    nvertices ;     /* # of vertices in the defect */
  int    *border ;
  int    nborder ;
  int    *chull ;        /* vertices in the convex hull */
  int    nchull ;        /* # of vertices in the convex hull */

  EDGE *edges;    /* the set of edges constituting the border of the defect */
  int nedges;     /* inside is left of vn1->vn2 - outside is right */

  int    *vertex_trans;
  int defect_number;

  int intersect; /* do final surface self-intersect ? */
  int optimal_mapping; /* type of mapping */
  float fitness; /* fitness of the retessellated patch */

  /* intensity information */
  float white_mean,white_sigma,white_mean_ll;
  float gray_mean,gray_sigma,gray_mean_ll;

  float k1_mean,k1_var;
  float k2_mean,k2_var;

  /* fitness information */
  /* initial fitness */
  float initial_face_ll,initial_vertex_ll,initial_curv_ll,initial_qcurv_ll,initial_mri_ll,initial_unmri_ll;
  /* final fitness */
  float final_face_ll,final_vertex_ll,final_curv_ll,final_qcurv_ll,final_mri_ll,final_unmri_ll;

}
DEFECT ;

#define MAX_DEFECTS             25000

typedef struct
{
  int    ndefects ;
  DEFECT defects[MAX_DEFECTS] ;
}
DEFECT_LIST, DL ;

/* structure which contains the information about a specific retessellation */
typedef struct
{
  int *vertices; /* list of used vertices in the defect (first inside ones) */
  int nvertices; /* nvertices in the list */
  int ninside;   /* # of inside vertices (not on the border) */
  int ndiscarded; /* number of discarded vertices */

  int *faces;    /* list of used faces in the defect */
  int nfaces;

  int *edges;    /* list of edges in the defect */
  int nedges;

  /* fitness information */
  float face_ll,vertex_ll,curv_ll,qcurv_ll,mri_ll,unmri_ll;
  //the next fitness information takes into account the local area
  float fll,vll,cll,qcll,mrill,unmrill;

}
TESSELLATED_PATCH, TP ;




#endif
