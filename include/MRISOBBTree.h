/**
 * @brief Constructs an Oriented Bounding Box datastructure from a MRIS surface
 *
 * Goal is to do a fast Point Inclusion Test
 */
/*
 * Original Author: Krish Subramaniam
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

// include guard
#ifndef MRISOBBTree_h
#define MRISOBBTree_h

#include "mri.h"
#include "mrisurf.h"
#include "diag.h"
#include "error.h"

#include <list>
#include <iostream>
#include <vector>
#include <stack>
#include <limits>

#define export 	// These headers use a deprecated / obsolete "export template" feature 
		// and gcc 5 and higher emit error messages
		// https://stackoverflow.com/questions/5416872/using-export-keyword-with-templates
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include "utilsmath.h"
#undef export

typedef Math::Point<double> Pointf;

//! Class which represents a single OBB Node. 
/** Has a corner and the axis of the OBB. Since it's a hierarchical 
 * data-structure, it has child pointers */
class OBBNode
{
  public:
    double corner[3];
    double axes[3][3];           // the axes defining OBB - ordered from long axis to short
    OBBNode *left;               // both left and right are NULL implies leaf node
    OBBNode *right;              // 

    int *faceidlist;             // list of faces in the OBBnode
    Pointf *normals;             // corresponding normals to the face
    int numfaces;                // number of faces in the OBBnode ( length of above )

    OBBNode()
    {
      this->left = NULL;
      this->right = NULL;
      this->faceidlist = NULL;
      this->normals = NULL;
    }
    ~OBBNode()
    {
      if ( faceidlist )
        delete [] faceidlist;
      if ( normals )
        delete [] normals;
    }
};

//! Class which represents the entire OBB tree
/** and also has methods to test for OBB intersection with a line 
 * and point inclusion methods */
class MRISOBBTree
{
  private:
    int level;                // used for bookkeeping
    int DeepestLevel;         // the deepest level of the tree
    int maxlevel;             // the maximum level allowed during the creation of OBBTree
    int OBBcount;             // number of OBB nodes in the trees
    int *InsertedPoints;      
    double Tolerance;         // the epslion. usually it's zero
    typedef std::vector<Pointf> PointsList;
    PointsList pointslist;    // the list of verts temporarily used while constructing the tree 

    MRIS *mris;               // the surface
    OBBNode *tree;            // Root node
    Math::AABB *aabb;         // the Axis Aligned bounding box for the surface

  public:
    MRISOBBTree(MRIS *_mris): mris(_mris)
    {
      this->level = 4; //start with an arbit level more than deepest level ( for recursive building )
      this->DeepestLevel=0;
      this->maxlevel = 12; // just an intelligent guess
      this->OBBcount = 0;
      this->Tolerance = 0.0;
      // init lists
      this->InsertedPoints = new int[mris->nvertices];
      for (int i=0; i<mris->nvertices; i++)
        this->InsertedPoints[i] = 0;
      // init the pointers to NULL
      aabb = NULL;
      tree = NULL;
    }

    ~MRISOBBTree()
    {
      if ( this->tree )                // delete the tree
      {
        this->FreeTree(this->tree);
      }
      if ( !this->pointslist.empty() ) // the temp list of verts  used while computing OBB
        this->pointslist.clear();
      if ( this->InsertedPoints )      // the list of inserted points in every OBB
      {
        delete [] this->InsertedPoints;
        this->InsertedPoints = NULL;
      }
        
      if ( aabb )
        delete aabb;                   // the AABB
      MRISfree(&mris);                 // Finally, the surface
    }

    //! Recursively free the tree
    void FreeTree(OBBNode *pOBB)
    {
      if ( pOBB != NULL )
      {
        FreeTree(pOBB->left);
        FreeTree(pOBB->right);
        free(pOBB); 
      }
    }

    //! Setter for the maximum level of tree
    void SetMaxLevel(int ml){ this->maxlevel = ml; }
    //! Getter for the maximum level of tree
    int GetMaxLevel(){ return this->maxlevel; }
    
    //! Getter for nodes count
    int GetNodesCount(){ return this->OBBcount; }

    //! Getter for the deepest level
    int GetDeepestLevel(){ return this->DeepestLevel; }

    //! Get root node pointer. Check for NULL when using this
    OBBNode* GetRootNode() 
    {
      if ( this->tree )
        return this->tree;
      else
        return NULL;
    }

    //! The main method which constructs the OBB Tree for the mris surface
    void ConstructTree()
    {
      //create an axis aligned bounding box used to trivially reject points in inclusion tests
      this->aabb=new Math::AABB(this->mris, 0.0); 

      // init the face id list with the face numbers sequentially
      std::vector<int> faceidlistvector;
      faceidlistvector.reserve(mris->nfaces);
      for ( int i=0; i < mris->nfaces; i++)
        faceidlistvector.push_back(i); 

      /*if ( this->tree )
      {
         this->FreeTree(this->tree);
         delete this->tree;
      }*/
      
      // root node
      this->tree = new OBBNode;
      this->DeepestLevel = 0;
      // Recursive building of the tree
      this->BuildTreeRecursively( faceidlistvector, this->tree, 0);
      this->level = this->DeepestLevel;

      //std::cerr << "OBBCount " << this->OBBcount << std::endl;
      //std::cerr << "Level " << this->level << std::endl;
      
      if ( this->InsertedPoints )
      {
        delete [] this->InsertedPoints;
        this->InsertedPoints = NULL;
      }

      this->pointslist.clear();
    }
   
    //! Build an OBB tree
    /*!
     \param faceids - the list of faces to consider
     \param pOBB - the current node
     \param _level - the current level
     When called, all the faces, the root node and 0 ( topmost level ) are given as the arguments
     */
    void BuildTreeRecursively(std::vector<int> faceids, OBBNode* pOBB, int _level)
    {
      int numfaces = faceids.size();
      
      if ( _level > this->DeepestLevel )
        this->DeepestLevel = _level;

      // Compute OBB 
      this->ComputeOBB(faceids, pOBB);

      if ( _level < this->maxlevel )
      {
        std::vector<int> L_faceids; L_faceids.reserve(numfaces/2);
        std::vector<int> R_faceids; R_faceids.reserve(numfaces/2);
        double n[3], p[3], c[3], x[3], val, ratio, bestratio;
        int negative, positive, splitacceptable, splitplane;
        int foundbestsplit, bestplane=0;
        int numinLnode, numinRnode, i;
        FACE *face;
        VERTEX *vertex;


        // loop over the 3 split planes to find the best one for division
        for (i=0; i<3; i++)
          p[i] = pOBB->corner[i] + pOBB->axes[0][i]/2.0 + pOBB->axes[1][i]/2.0 + pOBB->axes[2][i]/2.0;
        bestratio       = 1.0;          // worst case ratio
        foundbestsplit  = 0;            // not yet found
        splitplane      = 0;
        splitacceptable = 0;
        while ( !splitacceptable && splitplane < 3)
        {
          for (i=0; i<3; i++)
          {
            n[i] = pOBB->axes[splitplane][i];
          }
          Math::Normalize(n);

          //traverse the faces, and populate the appropriate child list as needed
          std::vector<int>::iterator fv;
          for (fv=faceids.begin(); fv != faceids.end(); ++fv)
          {
            face = &mris->faces[*fv];
            c[0] = c[1] = c[2] = 0.0;
            negative = positive = 0;
            for (i=0; i < 3; i++)
            {
              vertex = &mris->vertices[face->v[i]];
              x[0]   = vertex->x;
              x[1]   = vertex->y;
              x[2]   = vertex->z;
              val    = n[0] * (x[0]-p[0]) + n[1]*(x[1]-p[1]) + n[2]*(x[2]-p[2]);
              c[0]  += x[0];
              c[1]  += x[1];
              c[2]  += x[2];
              if (val < 0.0) negative = 1;
              else positive = 1;
            }

            if ( negative && positive )
            {
              // straddle case, so use centroid to decide
              c[0] = c[0] / 3.0;
              c[1] = c[1] / 3.0;
              c[2] = c[2] / 3.0;
              if ( n[0]*(c[0]-p[0]) + n[1]*(c[1]-p[1]) + n[2]*(c[2]-p[2]) < 0.0 )
                L_faceids.push_back(*fv);
              else
                R_faceids.push_back(*fv);
            }
            else
            {
              if (negative)
                L_faceids.push_back(*fv);
              else
                R_faceids.push_back(*fv);
            }
          }//for each face

          //evaluate the split
          numinLnode = L_faceids.size();
          numinRnode = R_faceids.size();
          ratio      = fabs(((double)numinRnode - numinLnode) / numfaces );
          // did we find acceptable split plane?
          if ( ratio < 0.6 || foundbestsplit ) // yes, so accept this
          {
            splitacceptable = 1;
          }
          else
          { // try another split
            L_faceids.clear(); L_faceids.reserve ( numfaces/2 );
            R_faceids.clear(); R_faceids.reserve ( numfaces/2 );
            if ( ratio < bestratio )
            {
              bestratio = ratio;
              bestplane = splitplane;
            }
            ++splitplane;
            if ( splitplane == 3 && bestratio < 0.95 )
            { // simply choose the best plane of the 3 split planes
              splitplane = bestplane;
              foundbestsplit = 1;
            }
          } // try next split
        } // for each split
        
        if ( splitacceptable )
        {
          OBBNode *LNode = new OBBNode;
          OBBNode *RNode = new OBBNode;
          pOBB->left    = LNode;
          pOBB->right   = RNode;
          
          faceids.clear(); //don't need them anymore

          // recurse
          this->BuildTreeRecursively(L_faceids, LNode, _level+1 );
          this->BuildTreeRecursively(R_faceids, RNode, _level+1 );
        }
        else
        {
          // free up objects
          L_faceids.clear();
          R_faceids.clear();
        }
      }

      FACE *face;
      VERTEX *v;
      double v1[3], v2[3], v3[3];
      pOBB->numfaces   = faceids.size();
      pOBB->faceidlist = new int[numfaces];
      pOBB->normals    = new Pointf[numfaces];
      for ( int i=0; i < numfaces; i++)
      {
       pOBB->faceidlist[i] = faceids[i];
       // Precalculate normals and store in every node
       face  = &mris->faces[faceids[i]];
       v = &mris->vertices[face->v[0]];
       v1[0] = v->x; v1[1] = v->y; v1[2] = v->z;
       v = &mris->vertices[face->v[1]];
       v2[0] = v->x; v2[1] = v->y; v2[2] = v->z;
       v = &mris->vertices[face->v[2]];
       v3[0] = v->x; v3[1] = v->y; v3[2] = v->z;
       Math::ComputeNormal(v1, v2, v3, pOBB->normals[i].v);
      }
       faceids.clear();
    } 

    //! Given the list of faces, compute OBB bounds for them and fill the pOBB object
    /*!
     \param faceids - the list of faces
     \param pOBB - the node which will be filled with values after the OBB is computed
     */
    void ComputeOBB( std::vector<int> faceids, OBBNode *pOBB)
    {
      int i, j, k; //counter vars
      double tri_mass, tot_mass, t;
      double mean[3], a0[3], a1[3], a2[3], c[3], p[3], q[3], r[3];
      double dp0[3], dp1[3], xp[3], evalues[3], tMin[3], tMax[3];
      double max[3], mid[3], min[3]; //store the axis of max, mid and min direction resp.
      vnl_matrix<double> a(3,3);
      VERTEX *v, *vertex;
      Pointf vp;
      FACE *face = NULL;

      tot_mass = 0.0;
      this->OBBcount++;
      this->pointslist.clear();
      this->pointslist.reserve(mris->nvertices);

      //Compute mean & moments 
      mean[0]=mean[1]=mean[2]=0.0;
      tot_mass = 0.0;
      for (i=0; i<3; i++)
        a0[i] = a1[i] = a2[i] = 0.0;
 
      // for each face
      std::vector<int>::iterator fv;
      for (fv=faceids.begin(); fv != faceids.end(); ++fv ) 
      {
        face  = &mris->faces[*fv];
        v = &mris->vertices[face->v[0]];
        p[0]  = v->x; p[1]  = v->y; p[2] = v->z;
        v = &mris->vertices[face->v[1]];
        q[0]  = v->x; q[1]  = v->y; q[2] = v->z;
        v= &mris->vertices[face->v[2]];
        r[0]  = v->x; r[1]  = v->y; r[2] = v->z;
        for (k=0; k<3; k++)
        {
          // edge vectors 
          dp0[k] = q[k] - p[k];
          dp1[k] = r[k] - p[k];
          // centroids
          c[k]   = (p[k] + q[k] + r[k])/3.0;
        }
        // cross product
        Math::Cross( dp0, dp1, xp);
        tri_mass = 0.5 * Math::Norm(xp);
        tot_mass += tri_mass;
        // mean
        for(k=0; k<3; k++)
          mean[k] += tri_mass * c[k];
        // diagonal terms
        a0[0] += tri_mass*(9*c[0]*c[0] + p[0]*p[0] + q[0]*q[0] + r[0]*r[0])/12;
        a1[1] += tri_mass*(9*c[1]*c[1] + p[1]*p[1] + q[1]*q[1] + r[1]*r[1])/12;
        a2[2] += tri_mass*(9*c[2]*c[2] + p[2]*p[2] + q[2]*q[2] + r[2]*r[2])/12;
        // off-diagonal terms
        a0[1] += tri_mass*(9*c[0]*c[1] + p[0]*p[1] + q[0]*q[1] + r[0]*r[1])/12;
        a0[2] += tri_mass*(9*c[0]*c[2] + p[0]*p[2] + q[0]*q[2] + r[0]*r[2])/12;
        a1[2] += tri_mass*(9*c[1]*c[2] + p[1]*p[2] + q[1]*q[2] + r[1]*r[2])/12;
        // Insert the points
        for ( j=0; j<3; j++)
        {
          if ( this->InsertedPoints[face->v[j]] != this->OBBcount )
          {
            this->InsertedPoints[face->v[j]] = this->OBBcount;
            vertex  = &mris->vertices[face->v[j]];
            vp.v[0] = vertex->x;
            vp.v[1] = vertex->y;
            vp.v[2] = vertex->z;
            this->pointslist.push_back(vp);
          }
        } //end insert points
      } // for each face

      // check area of face
      if (iszero(tot_mass)) {
        printf("error: triangle face defined by vertices %d, %d, and %d has no area!\n", face->v[0], face->v[2], face->v[2]);
        exit(1);
      }
      
      //normalize data
      for (i=0; i<3; i++)
        mean[i] = mean[i]/ tot_mass;
      //make matrix symmetric
      a1[0] = a0[1];
      a2[0] = a0[2];
      a2[1] = a1[2];
      // a0, a1 and a2 are actually column vectors of a
      for (i=0; i<3; i++)
      {
        a[0][i] = a0[i];
        a[1][i] = a1[i];
        a[2][i] = a2[i];
      }
      // get covariance 
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          a[i][j] = a[i][j]/tot_mass - mean[i]*mean[j];
      // find the eigenvectors -- these are got back in max, mid and min arrays
      vnl_symmetric_eigensystem<double> eigenSystem(a);
      Math::GetSortedEigenVectors ( eigenSystem, evalues, max, mid, min);
      for (i=0; i<3; i++)
      {
        a[0][i] = mean[i] + max[i];
        a[1][i] = mean[i] + mid[i];
        a[2][i] = mean[i] + min[i];
      }

      // create OBBs by projecting points onto eigenvectors
      // with the help of "t", the OBB is made to fit very tightly with the mesh
      tMin[0] =  tMin[1] = tMin[2] = std::numeric_limits<float>::max() ;
      tMax[0] =  tMax[1] = tMax[2] = -std::numeric_limits<float>::max() ;
      PointsList::iterator iv;
      for (iv=pointslist.begin(); iv != pointslist.end(); ++iv ) 
      {
        for (i=0; i<3; i++)
        {
          t = Math::DistanceToLine(iv->v, mean, a[i]);
          if (t < tMin[i])
            tMin[i] = t;
          if (t > tMax[i])
            tMax[i] = t;
        }
      }
      for (i=0; i<3; i++)
      {
        pOBB->corner[i]  = mean[i] +tMin[0]*max[i] + tMin[1]*mid[i] + tMin[2]*min[i];
        pOBB->axes[0][i] = (tMax[0] - tMin[0]) * max[i];
        pOBB->axes[1][i] = (tMax[1] - tMin[1]) * mid[i];
        pOBB->axes[2][i] = (tMax[2] - tMin[2]) * min[i];
      }
    
    } // end ComputeOBB method 

    //! Check whether the given point is inside or outside
    /*!
     \param x, y and z - the test point 
     \return 1 if the point lies inside the surface
     \return -1 if point lies outside the surface
     \return 0 if no intersection is found
     */
    int PointInclusionTest(double x, double y, double z)
    {
    
      // First check whether the point can be rejected trivially by just checking its AABB
      // AABB test is very fast and a large number of pts can be rejected this way.
      if ( this->aabb->CheckOutside(x, y, z) )
        return -1;

      // the idea is to construct a ray that is guaranteed to hit one of the faces
      // and use that as pointinclusion test
      VERTEX *v;
      FACE *face;
      double pt1[3], pt2[3], pt3[3];
      for (int fno=0; fno < mris->nfaces; fno++)
      {
        face = &mris->faces[fno];
        v = &mris->vertices[face->v[0]];
        pt1[0] = v->x; pt1[1] = v->y; pt1[2] = v->z;
        v = &mris->vertices[face->v[1]];
        pt2[0] = v->x; pt2[1] = v->y; pt2[2] = v->z;
        v = &mris->vertices[face->v[2]];
        pt3[0] = v->x; pt3[1] = v->y; pt3[2] = v->z;

        double pt[3];
        pt[0] = (pt1[0] + pt2[0] + pt3[0])/3;
        pt[1] = (pt1[1] + pt2[1] + pt3[1])/3;
        pt[2] = (pt1[2] + pt2[2] + pt3[2])/3;
        // the following line is guaranteed to pass through the face
        pt[0] += pt[0] - x;
        pt[1] += pt[1] - y;
        pt[2] += pt[2] - z;
        //calculate the vector from the calculated point to the test point
        double x_pt[3];
        x_pt[0] = pt[0] - x;
        x_pt[1] = pt[1] - y;
        x_pt[2] = pt[2] - z;

        //the face shouldn't be parallel to this line
        double n[3];
        //n[0] = face->nx; n[1] = face->ny; n[2] = face->nz;
        Math::ComputeNormal(pt1, pt2, pt3, n);
        double dotprod = Math::Dot(n, x_pt);
        if ( dotprod < 0 )
          dotprod = -dotprod;
        if (dotprod >= this->Tolerance + 1e-6 )
          return this->DoesLineIntersectOBBTree(x, y, z, pt);
        //better luck with next face?
      }
      return 0;
    }

    //! Check  if the line segment connecting pts (x,y,z) and pt[3] intersects the OBB Tree
    /*!
     \param x, y and z - starting pt of the segment
     \param pt[3] - end point of the segment
     \return 1 if the segment intersects with the tree
     \return -1 if segment doesn't intersect
     */
    int DoesLineIntersectOBBTree(double x, double y, double z, double pt[3])
    {
      int returnval = 0, faceno;
      OBBNode **obbstack;
      OBBNode* node;
      FACE* face;
      VERTEX *v;
      double v1[3], v2[3], v3[3];
      int sense, listsize=0, listmaxsize=10;
      double *distancelist = new double[listmaxsize]; // stores the parameter t from the point to faces
      int *senselist       = new int[listmaxsize];    // stores the sense ( 1 if entering or -1 if exiting )
      double distance;

      // reserve spaces
      obbstack = new OBBNode *[this->maxlevel+1];
      //compute vector from pt to (x,y,z)
      double v12[3], p1[3];
      p1[0] = x;
      p1[1] = y;
      p1[2] = z;
      v12[0] = pt[0] - x;
      v12[1] = pt[1] - y;
      v12[2] = pt[2] - z;

      // preorder traversal of the OBBTree
      // simulate recursion
      obbstack[0] = this->tree;
      int depth   = 1;

      while ( depth > 0 )
      {
        --depth;
        node = obbstack[depth];
        // Check if the line drawn from current point to the definitely inside point 
        // intersects the node
        if ( this->DoesOBBoverlapLine(node, p1, pt) )
        {
          // leaf node, so get the faces and test for inclusion
          if ( node->left == NULL && node->right == NULL ) 
          {
            for (int fv=0; fv < node->numfaces; ++fv)
            {
              faceno = node->faceidlist[fv]; 
              face  = &mris->faces[faceno];
              v = &mris->vertices[face->v[0]];
              v1[0] = v->x; v1[1]  = v->y; v1[2] = v->z;
              v = &mris->vertices[face->v[1]];
              v2[0] = v->x; v2[1]  = v->y; v2[2] = v->z;
              v = &mris->vertices[face->v[2]];
              v3[0] = v->x; v3[1]  = v->y; v3[2] = v->z;

              if ( LineIntersectsTriangle(p1, pt, v1, v2, v3, node->normals[fv].v, distance, sense) <= 0 )
                continue; // line doesn't intersect the triangle

              // line intersects the triangle!
              if ( listsize >= listmaxsize)
              {
                listmaxsize            *= 2;
                double *tmpdistancelist = new double[listmaxsize];
                int *tmpsenselist       = new int[listmaxsize];
                for ( int k=0; k < listsize; k++)
                {
                  tmpdistancelist[k] = distancelist[k];
                  tmpsenselist[k] = senselist[k];
                }
                delete [] distancelist;
                distancelist = tmpdistancelist;
                delete [] senselist;
                senselist = tmpsenselist;
              }
              distancelist[listsize] = distance;
              senselist[listsize]    = sense;
              ++listsize;
            }
          }
          else
          {
            obbstack[depth] = node->left; 
            obbstack[depth+1] = node->right;
            depth += 2;
          }
        }
      } // end traversal of the tree

      if ( listsize != 0  )
      {
        double ptol = this->Tolerance / Math::Norm(v12);
        double lastdistance = 0.0;
        double lastsense = 0;
        int listremainder = listsize;
        while ( listremainder )
        {
          int minIdx = 0;
          for ( int j =1 ; j < listremainder; j++)
          {
            if ( senselist[j] != lastsense && distancelist[j] < distancelist[minIdx] )
              minIdx = j;
          }
          distance = distancelist[minIdx];
          sense    = senselist[minIdx];
          listremainder--;
          distancelist[minIdx] = distancelist[listremainder];
          senselist[minIdx]    = senselist[listremainder];

          if ( distance > lastdistance - ptol && sense != lastsense)
          {
            //return value according to the sense of the first intersection
            if ( returnval == 0 )
              returnval = sense;
            lastdistance = distance;
            lastsense    = sense;
          }
        }
      } // if distancelist is not empty

      delete [] distancelist;
      delete [] senselist;
      delete [] obbstack;
      // 1 if point is inside, 0 if point is outside
      return returnval;
    }
    
    //! Check if a line intersects an OBBNode. Returns 1 if even at-least some portion
    //  of the line lies within the OBB
    /*!
     \param pOBB - the OBB to test
     \param p0 - starting point
     \param p1 - ending point of the line
     \returns 1 if the line intersects and 0 otherwise
     */
    bool DoesOBBoverlapLine(OBBNode* pOBB, double p0[3], double p1[3])
    {
      double rangeAmin, rangeAmax, rangePmin, rangePmax, dotP;
      //double  eps =  this->Tolerance;
      int i;
  
      for( i=0; i<3; i++)
      {
        rangeAmin = Math::Dot(pOBB->corner, pOBB->axes[i]);
        rangeAmax = rangeAmin + Math::Dot(pOBB->axes[i], pOBB->axes[i]);
        rangePmin   = Math::Dot(p0, pOBB->axes[i]);
        rangePmax   = rangePmin;
        dotP        = Math::Dot(p1, pOBB->axes[i]);
        if ( dotP < rangePmin ) 
          rangePmin = dotP;
        else
          rangePmax = dotP;

        if ( (rangeAmax < rangePmin) || (rangePmax < rangeAmin) )
          return 0; // doesn't overlap
      }
      return 1; // OBB overlaps line
    }

    //! Check whether the line p1-p2 intersects the triangle formed by v1, v2 and v3
    //! parametric coordinate 't' of the line is returned
    //! sense of intersection is also returned ( +1 if entering, -1 if exiting 
    //! according to the normal of the triangle)
    //! Algorithm based on the vtkOBBtreeLineIntersectsTriangle method
    /*!
     \param p1, p2 - the line p1-p2
     \param pt1, pt2, pt3 - the triangle 
     \param normal - precalculated normal of the triangle which faces "outside"
     \returns t - the parametric coordinate from which distance can be calculated
     \return sense - the sense ( see above
     */
    inline int LineIntersectsTriangle(double p1[3], double p2[3], 
                                double pt1[3], double pt2[3], double pt3[3], 
                                double normal[3],
                                double &t, int &sense)
    {

      //double normal[3];
      Math::ComputeNormal(pt1, pt2, pt3, normal);

      //vector from p1 to p2
      double v12[3];
      v12[0] = p2[0] - p1[0];
      v12[1] = p2[1] - p1[1];
      v12[2] = p2[2] - p1[2];

      //vector from p1 to triangle
      double v1t[3];
      v1t[0] = pt1[0] - p1[0];
      v1t[1] = pt1[1] - p1[1];
      v1t[2] = pt1[2] - p1[2];

      double num = Math::Dot(normal, v1t);
      double den = Math::Dot(normal, v12);

      if (den == 0)
        return 0; // triangle and the line are really parallel

      double fabsden = den;
      sense          = -1;
      if ( fabsden < 0.0 )
      {
        sense   = 1;
        fabsden = -fabsden;
      }
      if ( fabsden > 1e-6 +  this->Tolerance )
      {
        t = num / den ;
        if ( t < 0.0 || t > 1.0 )
          return 0;
        double point[3];
        point[0] = p1[0] + t * v12[0];
        point[1] = p1[1] + t * v12[1];
        point[2] = p1[2] + t * v12[2];
        
        // find axis permutation to do the rest of the math in 2D
        int xi = 0, yi = 1, zi = 2;
        if ( normal[0]*normal[0] < normal[1]*normal[1] )
        {
          xi = 1; yi = 2; zi = 0;
        }
        if ( normal[xi]*normal[xi] < normal[2]*normal[2])
        {
          xi = 2; yi = 0; zi = 1;
        }
        //calculate vector from the triangle corner to the point
        double u0 = point[yi] - pt1[yi];
        double v0 = point[zi] - pt1[zi];
        // calculate edge vectors for triangle
        double u1 = pt2[yi] - pt1[yi];
        double v1 = pt2[zi] - pt1[zi];
        double u2 = pt3[yi] - pt1[yi];
        double v2 = pt3[zi] - pt1[zi];
        // area of projected triangle (multiplied by 2) via cross product
        double area = (v2*u1 - u2*v1);
        // sub-areas that must sum to less than the total area
        double alpha = (v2*u0 - u2*v0);
        double beta = (v0*u1 - u0*v1);
        double gamma = area - alpha - beta;
        if (area < 0)
        {
          area = -area;
          alpha = -alpha;
          beta = -beta;
          gamma = -gamma;
        }

        // inside of polygon
        if (alpha > 0  &&  beta > 0  &&  gamma > 0)
          return 1;
      }

      // TODO -- if ( Tolerance != 0 ) part 
      return 0;
    }

  /*
  // implements Moller Trumbore Ray Triangle test ( 1997 )
  // doesn't seem to work right now
  inline int LineIntersectsTriangle(double p1[3], double p2[3], 
                                double pt1[3], double pt2[3], double pt3[3], 
                                double normal[3],
                                double &t, int &sense)
  {
    // find sense
    double dir[3];
    dir[0] = p2[0] - p1[0];
    dir[1] = p2[1] - p1[1];
    dir[2] = p2[2] - p1[2];

    double den = Math::Dot(normal, dir);
    if (den == 0)
      return 0; // triangle and the line are really parallel
    sense  = -1;
    if ( den < 0.0 )
        sense = 1;

    //rest of the algorithm
    //find vectors
    double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
    double det, inv_det, u, v;
    // find vectors for two edges sharing pt1
    edge1[0] = pt2[0] - pt1[0]; edge1[1] = pt2[1] - pt1[1]; edge1[2] = pt2[2] - pt1[2];
    edge2[0] = pt3[0] - pt1[0]; edge2[1] = pt3[1] - pt1[1]; edge2[2] = pt3[2] - pt1[2];

    // begin calculating det - also used to calculate U parameter 
    Math::Cross(dir, edge2, pvec);

    // if det is very close to 0, ray lies in the plane of the triangle
    det = Math::Dot(edge1, pvec);

    if ( det > -1e-9 && det < 1e-9 ) 
      return 0;
    inv_det = 1.0 / det;

    // dist from pt1 to ray origin
    tvec[0] = p1[0] - pt1[0]; tvec[1] = p1[1] - pt1[1]; tvec[2] = p1[2] - pt1[2];

    // calculate U param and test bounds
    u = Math::Dot(tvec, pvec) * inv_det;
    if ( u < 0.0 || u > det)
      return 0;

    //prepare to test V param
    Math::Cross(tvec, edge1, qvec);

    //calculate V param and test bounds
    v = Math::Dot(dir, qvec) * inv_det;
    if ( v < 0.0 || u + v > det)
      return 0;

    t = Math::Dot(edge2, qvec) * inv_det;

    return 1;
  }
*/
   
};

#endif
