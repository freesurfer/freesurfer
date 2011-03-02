/**
 * @file  fmarchmesh.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:43 $
 *    $Revision: 1.3 $
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


/* Program: fmarchmesh.c
 * Compute geodesic distance on a surface mesh using the fast-marching method
 * Author: Xiao Han
 */

/* The following function uses HEAP, which uses the LIST data structure */
#include "macros.h"
#include "mrisurf.h"
#include "fmarchmesh.h"
#define DEBUG 0

float ReCompute(int vIDc, int vIDa, int vIDb, MRI_SURFACE *mesh, float *T, unsigned char *label);
float ComputeTInAcute(float,float,float,float,float);

float *FastMarchMesh(MRI_SURFACE *mesh, int *contour, int numinitvert, float thred)
{
  /* Compute the geodesic distance of surface vertices to the initial contour*/
  /* vector 'contour' contains the indices of the vertice at the initial
     contour, i.e., vertices with zero distance. numinitvert is the size of
     contour.
     To save time, the distance computation stops when vertices with distance
     less than thred are all computed
  */
  /* The program assumes that the neighborhood structure of the surface mesh
     is already built, e.g., using function mrisFindNeighbors, which seems
     to be always performed during surface-file-reading
     Thus, the total number of neighors of vertex i is mesh->vertices[i].vnum,
     and the j-th neighbor of i is mesh->vertices[i].v[j].
     It also assumes the neighbors are stored in correct order, i.e.,
     two consecutive neighbors j and j+1 form a face with vertex i.
     Maybe a separate function need be written to make sure the neighborhood
     construction is correct!
     I think I may just use every face sharing the vertex directly, since the
     current program assumes the mesh is closed and simple, which may not be
     true. I think I need to use both, i.e., use face when need to compute
     distance from each shared face, and use v->v when need to enumerate
     every neighbors!
  */

  float *T;  /* final distances of each vertex to the contour */
  float tempvalue, newvalue;
  unsigned char *label; /* "Alive" => 1; "Narrow Band" => 2; "Far Away" => 3; */
  int *BackPointer;  /* backpointer to the narrowband heap */
  Xheap H;           /* Narrrow band heap */
  XheapElement he;
  //  int heapindex;

  FACE *face;

  int i, j, k, index, VN, VId, neighVId, n0, n1;

  VN = mesh->nvertices; /* total number of surface vertices */

  T = (float *) malloc(sizeof(float)*VN);

  for (i=0; i<VN; i++)
  {
    T[i] = INFINITY;  /*All distance initialized to a large value */
  }

  if (numinitvert <= 0)
  {
    printf("Warning, the initial contour is empty, no geodesic distance can be computed\n");

    return T;
  }

  label = (unsigned char *) malloc(sizeof(unsigned char)*VN);
  BackPointer = (int *) malloc(sizeof(int)*VN);

  /* Initialize heap structure */
  H = xhInitEmpty();

  /* Initialization for marching inwards and outwards simultaneously */
  for (i=0; i<VN; i++)
  {
    label[i] = (unsigned char) FAWAY;  /*All points are labelled as FAR AWAY */
  }

  for (i=0; i<numinitvert; i++)
  {
    label[contour[i]] = (unsigned char) ALIVE; /* all vertices on contour are alive */
    T[contour[i]] = 0; /* all vertices on contour have distance zero */
  }

  /* Initialize NarrowBand Heap */
  for (i=0; i<VN; i++)
  {
    if (label[i] == (unsigned char)ALIVE)
    {

      /* put its neighbors into the NarrowBand */
      if (DEBUG)
        printf("mesh->vertices[i].vnum = %d\n", mesh->vertices[i].vnum);
      for (j=0; j< mesh->vertices[i].vnum; j++)
      {
        neighVId = mesh->vertices[i].v[j];
        if (DEBUG) printf("neighVId = %d\n", neighVId);
        if (label[neighVId] == (unsigned char)FAWAY)
        {
          /* compute the distance to this point and add it to the heap */
          label[neighVId] = (unsigned char)NBAND;
          /* Check all possible triangles surrounding vertex neighVId1 */
          /* Note: Only ALIVE points contribute to the distance computation */
          newvalue = INFINITY;
          /* Compute distance from each face the vertex shares */
          /* k is the index of faces shared by the vertex neighVId */
          if (DEBUG)
            printf("mesh->vertices[neighVId].num = %d\n", mesh->vertices[neighVId].num);
          for (k=0; k<mesh->vertices[neighVId].num; k++)
          {
            index = mesh->vertices[neighVId].n[k]; /* # of vertex with ID
                                                           neighVId in the k-th face
                                                           that it is in */
            if (DEBUG)
              printf("index = %d\n", index);
            face = &mesh->faces[mesh->vertices[neighVId].f[k]]; /* ptr to the
                                                            k-th face */
            n0 = face->v[(index + VERTICES_PER_FACE -1)%VERTICES_PER_FACE];
            n1 = face->v[(index + 1)%VERTICES_PER_FACE];
            if (DEBUG) printf("no = %d, n1 = %d\n", n0, n1);
            tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label);
            if (tempvalue < newvalue)  newvalue = tempvalue;
          } /* end for-k */

          T[neighVId] = newvalue;
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);

        } /* end if */
      } /* end for */
    } /* end if */
  } /* end for */

  /* End of Initialization */

  /* Begin Fast Marching to get the unsigned distance function inwords
   * and outwards simultaneously
   * since points inside and outside the contour won't interfere with
   * each other
   */
  while (!xhIsEmpty(H))
  { /* There are still points not yet accepted */
    he = xhRemove(H); /* Label the point with smallest value among all
                                         NarrowBand points as ALIVE */
    VId = he.id;

    T[VId] = he.value;
    label[VId] = (unsigned char)ALIVE;

    if (he.value > thred) break;

    /* Update its neighbor */
    /* Put FARAWAY neighbors into NarrowBand, Recompute values at
     * NarrowBand neighbors,
     * Keep ALIVE (Accepted) neighbor unchanged
     */
    for (i=0; i<mesh->vertices[VId].vnum; i++)
    {
      neighVId = mesh->vertices[VId].v[i];

      /* Don't change ALIVE neighbors */
      if (label[neighVId] != (unsigned char)ALIVE)
      {
        newvalue = INFINITY;

        /* Compute distance from each face the vertex shares */
        /* j is the index of faces shared by the vertex neighVId */
        for (j=0; j<mesh->vertices[neighVId].num; j++)
        {
          face = &mesh->faces[mesh->vertices[neighVId].f[j]]; /* ptr to the
                                                  j-th face */

          index = mesh->vertices[neighVId].n[j]; /* index of the vertex
                                                   (neighVId) in the j-th
                                                   face that it belongs to */
          /* Take the ID of the remaining two vertices of the j-th face */
          n0 = face->v[(index + VERTICES_PER_FACE -1)%VERTICES_PER_FACE];
          n1 = face->v[(index + 1)%VERTICES_PER_FACE];
          tempvalue = ReCompute(neighVId, n0, n1, mesh, T, label);
          if (tempvalue < newvalue)  newvalue = tempvalue;
        } /* end for-j */

        /* If it was a FARAWAY point, add it to the NarrowBand Heap;
         * otherwise, just update its value
         * using the backpointer
         */
        if (label[neighVId] == (unsigned char)NBAND)
          xhChangeValue(BackPointer[neighVId],newvalue, H);
        else
        {
          xhInsert(newvalue, neighVId, &(BackPointer[neighVId]), H);
          label[neighVId] = (unsigned char)NBAND;
        }

      } /* end if ALIVE */
    } /* end updating neighbors */

  } /* end of marching loop */

  free(label);
  free(BackPointer);
  xhDestroy(H);

  return T;
}


/* vIDa, vIDb & vIDc are the three vertices making up the triangle */
/* mesh is the mesh to compute the distances on */
/* T is the current transit time values */
/* label is the label map */
/* compute the new T for vIDc based on the current T values of vIDa & vIDb */
/* Unfold the mesh around edge AB if C is an obtuse angle */
/* see Kimmel and Sethian for derivation of these formulae */
float ReCompute(int vIDc, int vIDa, int vIDb, MRI_SURFACE *mesh, float *T,
                unsigned char *label)
{
  float a, b, c;
  VERTEX *Va, *Vb, *Vc;
  unsigned char la, lb;
  float Ta, Tb;
  float t1, t2, t;

  int tmpvIDa;
  float tmpa;
  float tmpTa;
  unsigned char tmpla;

  /* Define some auxillary variables */
  int i, index, n0, n1, UNotID, UID, P1ID, P2ID;
  FACE *face;

  /* SQUARES of lengthes (Captital Letters) */
  float P1A, P1B, P1C, P2A, P2B, P2C, UA, UB, UC, P1P2, P1U, P2U;
  float cos12A, sin12A, cos12B, sin12B, cos12C, sin12C, cos12U, sin12U;
  float AC, BC, AB;

  VERTEX *P1, *P2, *U;
  MyVector vP1U, vP2U;

  int iters;

  la = label[vIDa];
  lb = label[vIDb];

  if ( (la != (unsigned char) ALIVE) && (lb != (unsigned char) ALIVE) )
    return INFINITY;

  Ta = T[vIDa];
  Tb = T[vIDb];

  Va = &mesh->vertices[vIDa];
  Vb = &mesh->vertices[vIDb];
  Vc = &mesh->vertices[vIDc];

  a = sqrt((Vb->x - Vc->x)*(Vb->x - Vc->x) + (Vb->y - Vc->y)*(Vb->y - Vc->y)
           + (Vb->z - Vc->z)*(Vb->z - Vc->z));
  b = sqrt((Va->x - Vc->x)*(Va->x - Vc->x) + (Va->y - Vc->y)*(Va->y - Vc->y)
           + (Va->z - Vc->z)*(Va->z - Vc->z));

  if (la != ALIVE) Ta = INFINITY;
  if (lb != ALIVE) Tb = INFINITY;

  if (Ta > Tb)
  {
    tmpTa = Ta;
    Ta = Tb;
    Tb = tmpTa;

    tmpla = la;
    la = lb;
    lb = tmpla;

    tmpa = a;
    a = b;
    b = tmpa;

    tmpvIDa = vIDa;
    vIDa = vIDb;
    vIDb = tmpvIDa;

    Va = &mesh->vertices[vIDa];
    Vb = &mesh->vertices[vIDb];
  }

  c = sqrt((Va->x - Vb->x)*(Va->x - Vb->x) + (Va->y - Vb->y)*(Va->y - Vb->y)
           + (Va->z - Vb->z)*(Va->z - Vb->z));

  AC = b*b;
  BC = a*a;
  AB = c*c;
  if (AB < AC+BC)
  {/* Acute Triangle */
    if (la == ALIVE && lb == ALIVE)
    {
      t = ComputeTInAcute(Ta, Tb, a,  b,  c);
      return t;
    }
    else  return (MIN(b+Ta, a+Tb)); /* Infinity + sth = Infinity */
  }

  /* Otherwise, perform unfolding */

  /*Initialization */
  P1ID = vIDa;
  P2ID = vIDb;
  UNotID = vIDc;

  P1A = 0;
  P1B = AB;
  P1C = AC;
  P2B = 0;
  P2A = P1B;
  P2C = BC;
  P1P2 = P1B;

  cos12A = 1;
  sin12A = 0;
  cos12B = 1;
  sin12B = 0;
  cos12C = (P1P2 + P2C - P1C)/(2*sqrt(P1P2*P2C));
  sin12C = (1-cos12C*cos12C); /* Notice: Square of sine */

  /* Now iteratively unfolding */
  iters = 0;
  while (iters <10)
  {
    /* Find the newly unfolded vertex ID */
    /* The following code assumes the neighbors form a circle, thus need be
     modified. But how?  */
    UID = -1;
    for (i=0; i<mesh->vertices[P1ID].num; i++)
    {
      index = mesh->vertices[P1ID].n[i]; /* # of vertex with ID
                                  P1ID in the i-th face
                                  that it is in */
      face = &mesh->faces[mesh->vertices[P1ID].f[i]]; /* ptr to the
                                 i-th face */
      n0 = face->v[(index + VERTICES_PER_FACE -1)%VERTICES_PER_FACE];
      n1 = face->v[(index + 1)%VERTICES_PER_FACE];
      if (n0 == P2ID && n1 != UNotID)
      {
        UID = n1;
        break;
      }
      else if (n1 == P2ID && n0 != UNotID)
      {
        UID = n0;
        break;
      }
    } /* end for-i */

    if (UID < 0) /* i.e., a valid UID cannot be found */
      break;

    P1 = &mesh->vertices[P1ID];
    P2 = &mesh->vertices[P2ID];
    U = &mesh->vertices[UID];

    vP1U.x = P1->x - U->x;
    vP1U.y = P1->y - U->y;
    vP1U.z = P1->z - U->z;
    vP2U.x = P2->x - U->x;
    vP2U.y = P2->y - U->y;
    vP2U.z = P2->z - U->z;

    P1U = (vP1U.x * vP1U.x + vP1U.y * vP1U.y + vP1U.z * vP1U.z);
    P2U = (vP2U.x * vP2U.x + vP2U.y * vP2U.y + vP2U.z * vP2U.z);

    cos12U = (P1P2 + P2U - P1U)/(2*sqrt(P1P2*P2U));
    sin12U = (1-cos12U*cos12U); /* Notice: Square of sine */

    /* Now compute three lengthes (squared) */
    UA = P2U + P2A - 2*sqrt(P2U*P2A)*(cos12A*cos12U - sqrt(sin12A*sin12U));
    UB = P2U + P2B - 2*sqrt(P2U*P2B)*(cos12B*cos12U - sqrt(sin12B*sin12U));
    UC = P2U + P2C - 2*sqrt(P2U*P2C)*(cos12C*cos12U - sqrt(sin12C*sin12U));

    /* Now Judge Which Side to continue unfolding */
    if (UA > (UC + AC))
    {/* Unfold along P1U */
      UNotID = P2ID;
      P2ID = UID;
      P1P2 = P1U;
      P2A = UA;
      P2B = UB;
      P2C = UC;
    }
    else if (UB > (UC + BC))
    { /* Unfold along P2U */
      UNotID = P1ID;
      P1ID = UID;
      P1P2 = P2U;
      P1A = UA;
      P1B = UB;
      P1C = UC;
    }
    else
    { /* Stop Unfolding and compute T*/
      /* Compute the actual lengthes */
      UC = sqrt(UC);
      UA = sqrt(UA);
      UB = sqrt(UB);
      if (label[UID] == (unsigned char) ALIVE)
      {
        if (la == ALIVE)
        {
          t1 = ComputeTInAcute(Ta, T[UID], UC, b, UA);
        }
        else t1 = INFINITY;
        if (lb == ALIVE)
        {
          t2 = ComputeTInAcute(Tb, T[UID], UC, a, UB);
        }
        else t2 = INFINITY;
        return MIN(t1,t2);
      }
      else return (MIN(b+Ta,a+Tb));
    }

    /* Update angles */
    cos12A = (P1P2 + P2A - P1A)/(2*sqrt(P1P2*P2A));
    if (P2B != 0) cos12B = (P1P2 + P2B - P1B)/(2*sqrt(P1P2*P2B));
    cos12C = (P1P2 + P2C - P1C)/(2*sqrt(P1P2*P2C));

    sin12A = 1 - cos12A*cos12A;
    sin12B = 1 - cos12B*cos12B;
    sin12C = 1 - cos12C*cos12C;

    iters++;
  }/* End of while loop */

  return (MIN(b+Ta, a+Tb));
}

float ComputeTInAcute(float Ta, float Tb, float a, float b, float c)
{
  float t1, t2, t, CD, costheta;
  float aa,bb,cc,u, tmp;

  costheta = (a*a+b*b-c*c)/(2*a*b);

  u = Tb - Ta;

  aa = a*a + b*b -2*a*b*costheta;
  bb = 2*b*u*(a*costheta-b);
  cc = b*b*(u*u-a*a*(1-costheta*costheta));

  tmp = bb*bb - 4*aa*cc;

  if (tmp < 0)   return (MIN(b+Ta,a+Tb));
  tmp = sqrt(tmp);

  t1 = (-bb + tmp)/(2*aa + 1e-15);
  t2 = (-bb - tmp)/(2*aa + 1e-15);
  t = MAX(t1,t2);
  CD = (b*(t-u))/t;

  if ( (u<t) && (a*costheta<CD) && (CD<(a/costheta)) )
    return (t+Ta);
  else
    return (MIN(b+Ta,a+Tb));

}
