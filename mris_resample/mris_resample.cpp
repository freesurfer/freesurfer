/**
 * @brief resample one surface to another, given their spherical regs
 *
 * This is a stripped-down version of mris_indirect_morph
 * 
 * Given two surfaces and their respective spherical registrations
 * (.sphere.reg files), this binary will resample the atlas onto the subject.
 *
 * The mechanism is the following:
 *  1. find for each vertex the closest vertex on the other surface.
 *  2. find the "closest" face
 *  3. do linear interpolation in spherical coordinates and use this
 *     interpolation on the original vertices
 *
 *  The additional functionality this binary provides is the
 *  sampling of an annotation file.
 *  The annotation uses the closest vertex principle derived
 *  in step 1. of the previous algorithm.
 */
/*
 * Original Author: Gheorghe Postelnicu, 2006
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 * onto the subject
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

// STL includes
#include <iostream>
#include <string>
#include <stdio.h>

// utilities
#include <ANN/ANN.h>

// FS includes
#include "argparse.h"
 
#include "error.h"
#include "mrisurf.h"


#include "mris_resample.help.xml.h"

const char* Progname;

// local fct declarations

struct IoParams;

void
process_input_data( const IoParams& params,
                    MRIS** mrisAtlasReg,
                    MRIS** mrisSubject );

void
resampleSurface( MRIS* mrisAtlasReg,
                 MRIS* mrisSubject,
                 bool passAnnotation = false);

double sqr(double x)
{
  return x*x;
}

double
v_to_f_distance(VERTEX *P0,
                MRI_SURFACE *mri_surf,
                int face_number,
                int debug,
                double *sopt,
                double *topt);

struct IoParams
{
  std::string strAtlasReg;

  std::string strSubjectReg;
  std::string strSubjectSurf;
  std::string strSubjectAnnot;
  bool annotation;

  std::string strOutput;
  std::string strOutAnnotation;

  void parse(int ac, char** av);
};

int
main(int ac, char** av)
{
  // parse cmd-line
  IoParams params;
  try
  {
    params.parse(ac,av);
  }
  catch (const std::exception& e)
  {
    std::cerr << "Exception caught while parsing cmd-line \n"
    << e.what() << std::endl;
    exit(1);
  }

  // input files
  MRIS* mrisAtlasReg;

  MRIS* mrisSubject;

  try
  {
    process_input_data( params,
                        &mrisAtlasReg,
                        &mrisSubject );
  }
  catch (const char* msg)
  {
    std::cerr << " Error processing input data "
    << msg << std::endl;
  }

  // perform resampling
  resampleSurface( mrisAtlasReg,
                   mrisSubject,
                   params.annotation );

  MRISrestoreVertexPositions( mrisAtlasReg, ORIGINAL_VERTICES);

  // copy the color table used for the initial volume
  mrisAtlasReg->ct = mrisSubject->ct;
  MRISwrite( mrisAtlasReg,
             const_cast<char*>(params.strOutput.c_str()));
  if ( params.annotation )
  {
    std::cout << " trying to write annotation\n";
    MRISwriteAnnotation
      ( mrisAtlasReg,
        const_cast<char*>( params.strOutAnnotation.c_str() ) );
  }
  return 0;
}

void
IoParams::parse(int ac, char** av)
{
  ArgumentParser parser;
  // required
  parser.addArgument("--atlas_reg", 1, String, true);
  parser.addArgument("--subject_reg", 1, String, true);
  parser.addArgument("--subject_surf", 1, String, true);
  parser.addArgument("--out", 1, String, true);
  // optional
  parser.addArgument("--annot_in", 1, String);
  parser.addArgument("--annot_out", 1, String);
  // help text
  parser.addHelp(mris_resample_help_xml, mris_resample_help_xml_len);
  parser.parse(ac, av);

  strAtlasReg = parser.retrieve<std::string>("atlas_reg");
  strSubjectReg = parser.retrieve<std::string>("subject_reg");
  strSubjectSurf = parser.retrieve<std::string>("subject_surf");
  strOutput = parser.retrieve<std::string>("out");

  std::string strAnnotationIn = parser.retrieve<std::string>("annot_in");
  std::string strAnnotationOut = parser.retrieve<std::string>("annot_out");

  annotation = false;

  if (strAnnotationOut.empty() && !strAnnotationIn.empty()) fs::fatal() << "missing --annot_out flag";
  if (!strAnnotationOut.empty() && strAnnotationIn.empty()) fs::fatal() << "missing --annot_in flag";

  if (!strAnnotationOut.empty() && !strAnnotationIn.empty()) {
    annotation = true;
    strSubjectAnnot = strAnnotationIn;
    strOutAnnotation = strAnnotationOut;
  }
}


void
process_input_data( const IoParams& params,
                    MRIS** mrisAtlasReg,
                    MRIS** mrisSubject )
{
  MRIS* mris;

  mris = MRISread
         ( const_cast<char*>(params.strAtlasReg.c_str()) );
  if (!mris) throw " Failed reading atlas registration file";
  *mrisAtlasReg = mris;

  mris = MRISread
         ( const_cast<char*>(params.strSubjectReg.c_str()));
  if (!mris) throw  " Failed reading subject spherical registration file";
  *mrisSubject = mris;

  MRISsaveVertexPositions(mris, TMP_VERTICES);

  if ( MRISreadVertexPositions
       ( mris, const_cast<char*>(params.strSubjectSurf.c_str()) )
       != NO_ERROR )
    throw " Failed reading subject orig surface file";

  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);

  if ( !params.strSubjectAnnot.empty() )
  {
    std::cout << " before reading annotation\n"
    << params.strSubjectAnnot << std::endl;
    if ( MRISreadAnnotation
         ( mris,
           const_cast<char*>(params.strSubjectAnnot.c_str()) )
         != NO_ERROR )
      throw " Failed reading annotation file for subject";
  }
}



void
resampleSurface( MRIS* mrisAtlasReg,
                 MRIS* mrisSubject,
                 bool passAnnotation )
{
  cheapAssert(mrisAtlasReg->origxyz_status == mrisSubject->origxyz_status);

  /* this function assumes the mapping function is stored
     as the "INFLATED vertex" of mris_template */

  int index, k, facenumber, closestface=-1;
  double tmps, tmpt; /* triangle parametrization parameters */
  double value, distance;
  double Radius, length, scale;
  double sumx, sumy, sumz, sumweight, weight;
  VERTEX *vertex, *V1, *V2, *V3;
  FACE *face;
  ANNpointArray QueryPt;

  VERTEX *pvtx;
  int numVertices;

  ANNpointArray pa = annAllocPts(mrisSubject->nvertices, 3);

  // populate the ANN point array
  pvtx = &( mrisSubject->vertices[0] );
  numVertices = mrisSubject->nvertices;
  for (index = 0;
       index < numVertices;
       ++pvtx, ++index)
  {
    pa[index][0] = pvtx->x;
    pa[index][1] = pvtx->y;
    pa[index][2] = pvtx->z;
  }

  // construct and initialize the tree
  ANNkd_tree *annkdTree = new ANNkd_tree(pa, mrisSubject->nvertices, 3);
  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];


  QueryPt = annAllocPts(1,3);

  bool normflag = false;
  if (normflag)
  {
    /* Compute the radius of the template mapping */
    vertex = &mrisSubject->vertices[0];
    Radius = sqrt(vertex->origx *vertex->origx + vertex->origy*vertex->origy
                  + vertex->origz*vertex->origz);
    printf("Radius = %g\n", Radius);
  }

  numVertices = mrisAtlasReg->nvertices;
  for (int index = 0;
       index < numVertices;
       ++index)
  {
    vertex = &mrisAtlasReg->vertices[index];
    
    QueryPt[0][0] = vertex->x;
    QueryPt[0][1] = vertex->y;
    QueryPt[0][2] = vertex->z;

    annkdTree->annkSearch( // search
      QueryPt[0],       // query point
      1,   // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,  // distance (returned)
      0);   // error bound

    /* annIndex gives the closest vertex on mris_template
       to the vertex #index of mris */
    /* Now need to find the closest face in order to
       perform linear interpolation */
    distance = 1000.0;


    if ( passAnnotation )
      vertex->annotation = mrisSubject->vertices[annIndex[0]].annotation;

    VERTEX_TOPOLOGY const * const vt = &mrisSubject->vertices_topology[annIndex[0]];
    
    for (k=0; k < vt->num; k++)
    {

      facenumber = vt->f[k]; /* index of the k-th face */
      if (facenumber < 0 || facenumber >= mrisSubject->nfaces) continue;
      value = v_to_f_distance(vertex, 
                              mrisSubject, 
                              facenumber, 
                              0, 
                              &tmps, 
                              &tmpt);

      if (distance > value)
      {
        distance = value;
        closestface = facenumber;
      }
    }

    //    distance = sqrt(distance); /* distance is no use here */
    /* Now interpolate the functions values to point (sopt,topt) */
    /* Linear interpolation is quite simple, just equal to
    sopt*vertex3 + topt*vertex2 + (1-s-t)*vertex1. check the function
    v_to_f_distance() to see how V1, V2, V3 is defined
    */
    /* Use the orig_vertex to store the interpolated cooridnate function */
    face = &mrisSubject->faces[closestface];
    V1 = &mrisSubject->vertices[face->v[0]];
    V2 = &mrisSubject->vertices[face->v[1]];
    V3 = &mrisSubject->vertices[face->v[2]];
    sumx = 0.0;
    sumy = 0.0;
    sumz = 0.0;
    sumweight = 0.0;
    weight = 1.0/(1e-20 + sqr(V1->x - vertex->x) +
                  sqr(V1->y - vertex->y) + sqr(V1->z - vertex->z) );
    sumx += weight*V1->origx;
    sumy += weight*V1->origy;
    sumz += weight*V1->origz;
    sumweight += weight;
    weight = 1.0/(1e-20 + sqr(V2->x - vertex->x) +
                  sqr(V2->y - vertex->y) + sqr(V2->z - vertex->z) );
    sumx += weight*V2->origx;
    sumy += weight*V2->origy;
    sumz += weight*V2->origz;
    sumweight += weight;
    weight = 1.0/(1e-20 + sqr(V3->x - vertex->x) +
                  sqr(V3->y - vertex->y) + sqr(V3->z - vertex->z) );
    sumx += weight*V3->origx;
    sumy += weight*V3->origy;
    sumz += weight*V3->origz;
    sumweight += weight;

    float origx = sumx /(sumweight + 1e-30);
    float origy = sumy /(sumweight + 1e-30);
    float origz = sumz /(sumweight + 1e-30);

    if (normflag)
    {
      length = sqrt(origx*origx + origy*origy + origz*origz);
      scale = Radius/(length + 1e-20);
      origx *= scale;
      origy *= scale;
      origz *= scale;
    }

    MRISsetOriginalXYZ(
      mrisAtlasReg, index,
      origx, origy, origz);
  }

  if (annkdTree) delete annkdTree;
  if (annIndex) delete[] annIndex;
  if (annDist) delete[] annDist;
  if (QueryPt) delete QueryPt;

  return;

}

double v_to_f_distance(VERTEX *P0,
                       MRI_SURFACE *mri_surf,
                       int face_number,
                       int debug,
                       double *sopt,
                       double *topt)
{
  /* Compute the distance from point P0 to a face of the surface mesh */
  /* (sopt, topt) determines the closest point inside the triangle from P0 */

  double a, b, c, d, e, f, det, s, t, invDet;
  double numer, denom, tmp0, tmp1;
  VERTEX *V1, *V2, *V3;
  FACE *face;

  struct { float x,y,z; } E0, E1, D;

  face = &mri_surf->faces[face_number];
  V1 = &mri_surf->vertices[face->v[0]];
  V2 = &mri_surf->vertices[face->v[1]];
  V3 = &mri_surf->vertices[face->v[2]];

  E0.x = V2->x - V1->x;
  E0.y = V2->y - V1->y;
  E0.z = V2->z - V1->z;
  E1.x = V3->x - V1->x;
  E1.y = V3->y - V1->y;
  E1.z = V3->z - V1->z;
  D.x = V1->x - P0->x;
  D.y = V1->y - P0->y;
  D.z = V1->z - P0->z;

  a = E0.x *E0.x + E0.y * E0.y + E0.z *E0.z;
  b = E0.x *E1.x + E0.y * E1.y + E0.z *E1.z;
  c = E1.x *E1.x + E1.y * E1.y + E1.z *E1.z;
  d = E0.x *D.x + E0.y * D.y + E0.z *D.z;
  e = E1.x *D.x + E1.y * D.y + E1.z *D.z;
  f = D.x *D.x + D.y * D.y + D.z *D.z;

  det = a*c - b*b;
  s = b*e - c*d;
  t = b*d - a*e;

  if (debug)
    std::cout << " det = " << det << std::endl;

  if (s + t <= det)
  {
    if (s < 0)
    {
      if (t<0)
      {
        /* Region 4 */
        if ( d < 0)
        { /* Minimum on edge t = 0 */
          s = (-d >= a ? 1 : -d/a);
          t = 0;
        }
        else if (e < 0)
        { /* Minimum on edge s = 0 */
          t = (-e >= c ? 1 : -e/c);
          s = 0;
        }
        else
        {
          s = 0;
          t = 0;
        }
        if (debug)
          printf("region 4,  d= %g, e = %g, s =%g, t =%g\n", d, e, s, t);
      }
      else
      {
        /* Region 3 */
        s = 0;
        t = ( e >= 0 ? 0 : (-e >= c ? 1 : (-e/c)));
        if (debug) printf("region 3, s =%g, t =%g\n", s, t);
      }
    }
    else if (t < 0)
    {
      /* Region 5 */
      t = 0;
      s = (d >= 0 ? 0 :(-d >= a ? 1 : (-d/a)));
      if (debug) printf("region 5, s =%g, t =%g\n", s, t);
    }
    else
    {
      /* Region 0 */
      invDet = 1/det;
      s *= invDet;
      t *= invDet;
      if (debug) printf("region 0, s =%g, t =%g\n", s, t);
    }

  }
  else
  {
    if ( s < 0 )
    {
      /* Region 2 */
      tmp0 = b + d;
      tmp1 = c + e;
      if (tmp1 > tmp0)
      {
        numer = tmp1 - tmp0;
        denom = a - b - b +c;
        s = (numer >= denom ? 1 : numer/denom);
        t = 1-s;

      }
      else
      {
        s = 0;
        /* t = (e >= 0 ? 0 : (-e >= c ? 0 > = c + e = tmp1) */
        t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
      }
      if (debug) printf("region 2, s =%g, t =%g\n", s, t);
    }
    else if ( t < 0)
    {
      /* Region 6 */
      tmp0 = b + e;
      tmp1 = a + d;
      if (tmp1 > tmp0)
      { /* Minimum at line s + t = 1 */
        numer = tmp1 - tmp0; /* Positive */
        denom = a + c - b -b;
        t = (numer >= denom ? 1 : (numer/denom));
        s = 1 - t;
      }
      else
      { /* Minimum at line t = 0 */
        s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d/a));
        t = 0;
      }
      if (debug) printf("region 6, s =%g, t =%g\n", s, t);
    }
    else
    {
      /* Region 1 */
      numer = c + e - b - d;
      if (numer <= 0)
      {
        s = 0;
      }
      else
      {
        denom = a + c - b - b; /* denom is positive */
        s = (numer >= denom ? 1 : (numer/denom));
      }
      t = 1-s;
      if (debug) printf("region 1, s =%g, t =%g\n", s, t);
    }
  }

  if ( s < 0  || s > 1 || t < 0 || t > 1)
  {
    printf("Error in computing s and t \n");
  }

  *sopt = s;
  *topt = t;
  /* return (sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f)); */
  /* The square-root will be taken later to save time */
  return (a* sqr(s) + 2*b*s*t + c* sqr(t) + 2*d*s + 2*e*t + f);

}
