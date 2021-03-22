#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/**
 * @brief surface topology correction utilities
 *
 * "Automated Manifold Surgery: Constructing Geometrically Accurate and
 * Topologically Correct Models of the Human Cerebral Cortex",
 * Fischl, B., Liu, A. and Dale, A.M.
 * (2001) IEEE Transactions on Medical Imaging, 20(1):70-80.
 */
/*
 * Original Author: Florence Segonne
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

#include "mris_topology.h"
#include "patchdisk.h"
#include "utils.h"

#include "compilerdefs.h"

#define __PRINT_MODE 0
#define WHICH_OUTPUT stderr

bool doesMRISselfIntersect(MRIS *mris_work,TOPOFIX_PARMS &parms);

//check the new vertices : val2, val2bak... marked2=-1 ?

bool MRIScorrectDefect(MRIS *mris, 
                                  int defect_number,
                                  TOPOFIX_PARMS &parms) 
{
  int euler = MRISgetEuler(mris, defect_number);
  fprintf(WHICH_OUTPUT,
          "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n"
          "Correcting Topology of defect %d with "
          "euler number %d (%d loops) \n",
          parms.defect_number,euler,(1-euler)/2);

  if (euler == 1) {
    fprintf(WHICH_OUTPUT,
            "   Nothing to correct for defect %d!!\n",
            parms.defect_number);
    return true;
  }

  MRISinitDefectParameters(mris,&parms);

#if __PRINT_MODE
  fprintf(WHICH_OUTPUT,"   extracting MRIP\n");
#endif

  MRIP *mrip = MRIPextractFromMRIS(mris,defect_number);
  parms.mrip=mrip;

  if (mrip==NULL) {
    TOPOFIXfreeDP(&parms);
    return false;
  }

#if __PRINT_MODE
  euler = MRISgetEuler(mrip->mris);
  fprintf(WHICH_OUTPUT,"BEFORE TOP euler is %d\n",euler);
#endif

  //making sure of a few important things before the topology correction
  MRISclearFaceMarks(mris);
  MRISclearMarks(mris);
  // MRIScomputeNormals(mris);
  // MRIScomputeMetricProperties(mris);
  MRIScomputeTriangleProperties(mris);

  bool correct = MRIScorrectPatchTopology(mrip->mris,parms);
  if (correct == false) {
    fprintf(WHICH_OUTPUT,"PBM : Could not correct topology\n");
    MRIPfree(&mrip);
    TOPOFIXfreeDP(&parms);
    return false;
  }

#if __PRINT_MODE
  fprintf(WHICH_OUTPUT,"transferring corrections\n");
#endif

  MRISaddMRIP(mris,mrip);

#if __PRINT_MODE
  euler = MRISgetEuler(mris);
  fprintf(WHICH_OUTPUT,  "Now, euler is %d\n",euler);
#endif

  MRIPfree(&mrip);
  TOPOFIXfreeDP(&parms);

  // fprintf(WHICH_OUTPUT,"\n");
  //printing out results
  if (parms.verbose) {
    fprintf(WHICH_OUTPUT,"INITIAL FITNESS :  %3.5f \n",parms.initial_fitness);
    fprintf(WHICH_OUTPUT,"FINAL FITNESS   :  %3.5f \n",parms.fitness);
    fprintf(WHICH_OUTPUT,"# generated patches : %d (%d self-intersecting)\n",
            parms.ngeneratedpatches,parms.nselfintersectingpatches);
  }

  return true;
}

bool MRIScorrectPatchTopology(MRIS* &mris,TOPOFIX_PARMS &parms) 
{
#if __PRINT_MODE
  fprintf(WHICH_OUTPUT,"Correcting Topology\n");
#endif

  //first initialize the right parameters for this defect
  //what is the maximum face generating a cutting loop
  parms.max_face = mris->nfaces;
  parms.nattempts = nint(parms.nattempts_percent*mris->nfaces+1);
  parms.nminimal_attempts=
    nint(parms.minimal_loop_percent*mris->nfaces+parms.nminattempts);

  MRIScomputeNormals(mris);
  MRIScomputeTriangleProperties(mris);
  MRISsaveVertexPositions(mris,ORIGINAL_VERTICES);
  MRISinitDefectPatch(mris,&parms);
  MRIScomputeInitialFitness(mris,&parms);
  parms.fitness = -1.0;
  parms.ngeneratedpatches=0;
  parms.nselfintersectingpatches=0;
  int nloops = (1-MRISgetEuler(mris))/2;
  //we limit the total number of attempts to 2000!
  parms.nattempts = __MIN(2000/nloops,parms.nattempts);
  //parms.nminimal_attempts = 
  // __MIN(__MAX(1000/nloops,500),parms.nminimal_attempts);


  for (int n = 0 ; n < nloops ; n++) {
    if (parms.verbose>=VERBOSE_MODE_MEDIUM) {
      fprintf(WHICH_OUTPUT,
              "      max face = %d(%d) - loop = %d (%d)  - ntries = [%d,%d]\n",
              parms.max_face,mris->nfaces,n+1,nloops, 
              parms.nattempts,parms.nminimal_attempts);
    }
    bool correct = MRISincreaseEuler(mris,parms);
    if (correct == false) return false;
  }
  return true;
}

#define WS 0


bool MRISincreaseEuler(MRIS* &mris,TOPOFIX_PARMS &parms) 
{
  int nattempts=parms.nattempts;
  int final_euler = -1;

  static int s_nbr = 0;

  char fname[STRLEN] ;

  int nintersections=0,npatches=0;

  MRIS *best_mris = NULL;
  double best_fitness = -1;

#if __PRINT_MODE
  fprintf(WHICH_OUTPUT,"increasing euler number\n");
#endif

  for (int n = 0 ; n < nattempts ; n++) {
    if (parms.verbose>=VERBOSE_MODE_HIGH)
      fprintf(WHICH_OUTPUT, 
              "\r      %3.2f%% [%5d,%5d]",
              n*100.0/nattempts,n,nattempts);

    //generate a surface
    Surface *surface = MRIStoSurface(mris);
    surface->disk=(PatchDisk*)parms.patchdisk;

    //compute the euler number of the surface (mandatory for CutPatch!)
    int euler = surface->GetEuler();
    int init_euler = euler;

#if __PRINT_MODE
    fprintf(WHICH_OUTPUT,
            "BEFORE %d (%d: %d, %d)\n",
            euler,n,surface->nvertices,surface->nfaces);
#endif

    //increase the euler number by 2
    int correct;
    if (n == 0 ) //the first time, try to cut with a very small patch
      correct = surface->CutPatch(-2,parms.max_face,parms.nminimal_attempts,1);
    else
      correct = surface->CutPatch(-1,
                                  parms.max_face,
                                  10,1);//always trying 10 times at least

    npatches++;

    if (correct < 0) {
      fprintf(WHICH_OUTPUT,
              "\r      PBM: Euler Number incorrect for surface (defect %d)\n",
              parms.defect_number);
      delete surface;
      continue;
    }
#if __PRINT_MODE
    fprintf(WHICH_OUTPUT,
            "AFTER %d (%d,%d)\n",
            surface->GetEuler(),surface->nvertices,surface->nfaces);
#endif


    //transfer data into MRIS structure
    MRIS *mris_work = SurfaceToMRIS(surface);
    delete surface;

    MRIScopyHeader(mris,mris_work);
    MRISsaveVertexPositions(mris_work, INFLATED_VERTICES); /*saving current 
                                                             vertex positions 
                                                             into inflated */
    parms.mris_defect=mris_work;

    euler = MRISgetEuler(mris_work);
    if (euler != init_euler+2) {
      fprintf(WHICH_OUTPUT,
              "\r      PBM: Euler Number incorrect for mris (defect %d)\n",
              parms.defect_number);
      MRISfree(&mris_work);
      continue;
    }

#if __PRINT_MODE
    fprintf(WHICH_OUTPUT,"Topology -> %d \n",euler);
#endif

    // we have a correct surface : evaluate if valid
    MRISinitDefectPatch(mris_work,&parms);

#if WS
    sprintf(fname,"./defect_%d_%d.asc",parms.defect_number,s_nbr++);
    //fprintf(stderr,"%s!\n",fname);
    MRISwrite(mris_work,fname);
#endif

    //first check if the patch self-intersects
    bool selfintersect = doesMRISselfIntersect(mris_work,parms);

		
		//new version of mris_topo_fixer: generating small corrections first!
		if(selfintersect){ //the previous surface self-intersected! -> generate another one
			nintersections++;
      if (parms.verbose>=VERBOSE_MODE_HIGH) fprintf(WHICH_OUTPUT,"\r      SELF-INTERSECTING PATCH\n");
			MRISfree(&mris_work);

			//generate a new surface
			Surface *surface = MRIStoSurface(mris);
			surface->disk=(PatchDisk*)parms.patchdisk;

			//compute the euler number of the surface (mandatory for CutPatch!)
			int euler = surface->GetEuler();
			int init_euler = euler;

#if __PRINT_MODE
			fprintf(WHICH_OUTPUT,
							"BEFORE %d (%d: %d, %d)\n",
							euler,n,surface->nvertices,surface->nfaces);
#endif

			//increase the euler number by 2
			int correct;
			if (n == 0 ){
        correct = surface->CutPatch(-2,parms.max_face,parms.nminimal_attempts);
      }
			else {
        //always trying 10 times at least
        correct = surface->CutPatch(-1,parms.max_face,10);
      }
			npatches++;

			if (correct < 0) {
				fprintf(
          WHICH_OUTPUT,
          "\r      PBM: Euler Number incorrect for surface (defect %d)\n",
          parms.defect_number);
				delete surface;
				continue;
			}
#if __PRINT_MODE
			fprintf(WHICH_OUTPUT,
							"AFTER %d (%d,%d)\n",
							surface->GetEuler(),surface->nvertices,surface->nfaces);
#endif

			//transfer data into MRIS structure
			mris_work = SurfaceToMRIS(surface);
			delete surface;

			MRIScopyHeader(mris,mris_work);
			MRISsaveVertexPositions(mris_work, INFLATED_VERTICES); /*saving current 
																														 vertex positions 
                                                             into inflated */
			parms.mris_defect=mris_work;

			euler = MRISgetEuler(mris_work);
			if (euler != init_euler+2) {
				fprintf(WHICH_OUTPUT,
								"\r      PBM: Euler Number incorrect for mris (defect %d)\n",
								parms.defect_number);
				MRISfree(&mris_work);
				continue;
			}

#if __PRINT_MODE
			fprintf(WHICH_OUTPUT,"Topology -> %d \n",euler);
#endif

			// we have a correct surface : evaluate if valid
			MRISinitDefectPatch(mris_work,&parms);

#if WS
			sprintf(fname,"./defect_%d_%d.asc",parms.defect_number,s_nbr++);
			//fprintf(stderr,"%s!\n",fname);
			MRISwrite(mris_work,fname);
#endif

			//check again if the patch self-intersects
			selfintersect = doesMRISselfIntersect(mris_work,parms);
		}
		// end of the new version of mris_topo_fixer
		////////////////////////////////////////////////////////////////////
    if (selfintersect) {
      nintersections++;
      if (parms.verbose>=VERBOSE_MODE_HIGH) 
        fprintf(WHICH_OUTPUT,"\r      SELF-INTERSECTING PATCH\n");
      if (parms.minimal_mode && best_mris==NULL) {
        nattempts = __MAX(50,nattempts);
      }
    } else {
      //compute associated fitness
      double fitness = MRIScomputeFitness(mris_work,&parms,1);
      //if(parms.verbose>=VERBOSE_MODE_MEDIUM) 
      //  fprintf(WHICH_OUTPUT,"\r      fitness (o)is %3.5f \n",fitness);
      //update if necessary
      if (best_mris == NULL || fitness > best_fitness) {
        if (best_mris) MRISfree(&best_mris);
        best_mris = MRISduplicateOver(mris_work,1);
        best_fitness = fitness;
        final_euler = euler;
        if (parms.verbose>=VERBOSE_MODE_MEDIUM) {
          fprintf(WHICH_OUTPUT,
                  "\r      BEST FITNESS (o)is %3.5f \n",best_fitness);
          MRISprintInfo(&parms);
        }
#if WS
        sprintf(fname,"./best_%d_%d.asc",parms.defect_number,s_nbr);
        //fprintf(stderr,"%s!\n",fname);
        MRISwrite(mris_work,fname);
#endif
        if (parms.write) {
          sprintf(fname,"./defect_%d_%d.asc",parms.defect_number,s_nbr++);
          //fprintf(stderr,"%s!\n",fname);
          MRISwrite(mris_work,fname);
        }
      }
    };

    if (!parms.smooth && !parms.match) {
      MRISfree(&mris_work);
      continue;
    }

    //now try to use local intensities to improve fitness
    //first mark vertices
    MRISmarkPatchVertices(mris_work,&parms,mris->nvertices);
    for (int n = 0 ; n < 3 ; n++) MRISexpandMarked(mris_work);
    MRISmarkBorderVertices(mris_work,&parms,0);
    MRISdefectMatch(mris_work,&parms);
    npatches++;

#if WS
    sprintf(fname,"./defect_%d_%d.asc",parms.defect_number,s_nbr++);
    //fprintf(stderr,"%s!\n",fname);
    MRISwrite(mris_work,fname);
    sprintf(fname,"./curv_%d_%d.asc",parms.defect_number,s_nbr);
    //fprintf(stderr,"%s!\n",fname);
    for (int p =0; p < mris_work->nvertices;p++)
      mris_work->vertices[p].curv=mris_work->vertices[p].marked;
    MRISwriteCurvature(mris_work,fname);

#endif

    //if(parms.verbose>=VERBOSE_MODE_MEDIUM) 
    //fprintf(WHICH_OUTPUT,"\r      fitness (m)is %3.5f \n",fitness);
    //check if the surface self-intersects
    selfintersect = doesMRISselfIntersect(mris_work,parms);
    if (selfintersect) {
      nintersections++;
      if (parms.verbose>=VERBOSE_MODE_HIGH) 
        fprintf(WHICH_OUTPUT,"\r      SELF-INTERSECTING PATCH\n");
      MRISfree(&mris_work);
      continue;
    }
    //compute the fitness
    double fitness = MRIScomputeFitness(mris_work,&parms,1);
    //update if necessary
    if (best_mris == NULL || fitness > best_fitness) {
      if (best_mris) MRISfree(&best_mris);
      best_mris = mris_work;
      best_fitness = fitness;
      final_euler=euler;
      if (parms.verbose>=VERBOSE_MODE_MEDIUM) {
        fprintf(WHICH_OUTPUT,
                "\r      BEST FITNESS (M) is %3.5f \n",best_fitness);
        MRISprintInfo(&parms);
      }
#if WS
      sprintf(fname,"./best_%d_%d.asc",parms.defect_number,s_nbr);
      //fprintf(stderr,"%s!\n",fname);
      MRISwrite(mris_work,fname);
#endif
      if (parms.write) {
        sprintf(fname,"./defect_%d_%d.asc" , parms.defect_number,s_nbr++);
        //fprintf(stderr,"%s!\n",fname);
        MRISwrite(mris_work,fname);
      }
    } else {
      MRISfree(&mris_work);
    }
  }
  if (parms.verbose>=VERBOSE_MODE_HIGH)
    fprintf(WHICH_OUTPUT,"\r                             \r");

  if (parms.verbose>=VERBOSE_MODE_MEDIUM)
    fprintf(WHICH_OUTPUT,
            "      %d patches have been generated - %d self-intersected\n",
            npatches,nintersections);

  //update the surface
  if (best_mris == NULL ) {
    fprintf(WHICH_OUTPUT,
            "PBM : Could Not Increase Euler Number for defect %d\n",
            parms.defect_number);
    return false;
  }

  if (final_euler==1) s_nbr=0;

  MRISfree(&mris);
  mris = best_mris;
  parms.fitness = best_fitness;
  parms.ngeneratedpatches += npatches;
  parms.nselfintersectingpatches += nintersections;

  return true;
}

bool doesMRISselfIntersect(MRIS *mris_work,TOPOFIX_PARMS &parms) 
{
  if (parms.no_self_intersections==false) return false;

  //checking self-intersections
  return IsMRISselfIntersecting(mris_work);
}

MRIS *MRISduplicateOver(MRIS *mris,int mode) 
{
  //clone the surface mris
  MRIS * mris_dst = MRISclone(mris);
  
  cheapAssert(mris_dst->origxyz_status == mris->origxyz_status);
  
  for (int n = 0 ; n < mris->nvertices ; n++) {
    VERTEX *vdst = &mris_dst->vertices[n];
    VERTEX *vsrc = &mris->vertices[n];
    vdst->marked2 = vsrc->marked2;
    vdst->marked  = 0 ;
    vdst->ripflag = 0 ;
    
    MRISsetOriginalXYZ(mris_dst, n,
      vsrc->origx,
      vsrc->origy,
      vsrc->origz);
    
    vdst->val = vsrc->val;
    vdst->val2 = vsrc->val2;
    vdst->val2bak = vsrc->val2bak;
  }

#if 0
  // ATH: I don't think it's necessary to reallocate vertices and faces here.
  // In fact, doing this ends up causing issues during patching, since too many
  // vertices get added.

  int n_extra_vertices,n_extra_faces;

  if (mode==0) {
    //count the number of loops to be corrected
    int nloops = (2-MRISgetEuler(mris))/2;

    //count the number of defective vertices
    n_extra_vertices = 0;
    for (int n = 0 ; n < mris->nvertices ; n++)
      if (mris->vertices[n].marked2) n_extra_vertices++;
    n_extra_vertices = 3*n_extra_vertices + nloops*MAX_EXTRA_VERTICES; //2

    //count the number of defective faces
    n_extra_faces = 0 ;
    for (int n = 0 ; n < mris->nfaces ; n++) {
      bool isface=true;
      for (int i = 0 ; i < 3 ; i++) {
        if (mris->vertices[mris->faces[n].v[i]].marked2==0) {
          isface = false;
          break;
        }
      }
      if (isface) n_extra_faces++;
    }
    n_extra_faces = 5*n_extra_faces + nloops*MAX_EXTRA_FACES; //4
  } else {
    n_extra_vertices = __MAX(mris->max_vertices - mris->nvertices, 0);
    n_extra_faces    = __MAX(mris->max_faces    - mris->nfaces,    0);
  }
  
  MRISreallocVerticesAndFaces(mris_dst, mris_dst->nvertices + n_extra_vertices, mris_dst->nfaces + n_extra_faces);
#endif

  return mris_dst;
}


int MRISgetEuler(MRIS *mris, int defect_number) 
{
  int *list_of_faces,nfaces;

  if (defect_number < 0) {
    //counting faces
    nfaces = mris->nfaces ;
    if (nfaces == 0) return 0;

    //allocate the list of faces
    list_of_faces = (int*)malloc(nfaces*sizeof(int));

    for (int n = 0 ; n < nfaces ; n++)
      list_of_faces[n]=n;
  } else {

    //counting faces
    nfaces = 0 ;
    for (int n = 0 ; n < mris->nfaces ; n++) {
      bool isface=true;
      for (int i = 0 ; i < 3 ; i++) {
        if (mris->vertices[mris->faces[n].v[i]].marked2!=defect_number) {
          isface=false;
          break;
        }
      }
      if (isface) nfaces++;
    }
    if (nfaces==0) return 0;

    //allocate the list of faces
    list_of_faces = (int*)malloc(nfaces*sizeof(int));

    //initialize the list of faces
    nfaces = 0 ;
    for (int n = 0 ; n < mris->nfaces ; n++) {
      bool isface = true;
      for (int i = 0 ; i < 3 ; i++) {
        if (mris->vertices[mris->faces[n].v[i]].marked2!=defect_number) {
          isface = false;
          break;
        }
      }
      if (isface) list_of_faces[nfaces++]=n;
    }
  }

  int euler = MRISgetEulerNumber(mris, list_of_faces,nfaces);

  delete [] list_of_faces;

  return euler;
}

int MRISgetEulerNumber(const MRIS *mris,  const int *list_of_faces, int nfs) 
{
  int nv,nf,ne;

  //we need to allocate the edge structure for the faces in the list
  //mark and count the vertices
  for (int n = 0 ; n < mris->nvertices ; n++) {
    mris->vertices         [n].marked = 0;
    mris->vertices_topology[n].e=0;
  }
  nv = 0 ;
  for (int n = 0 ; n < nfs ; n++)
    for (int i = 0 ; i < 3 ; i++) {
      VERTEX_TOPOLOGY * const vt=&mris->vertices_topology[mris->faces[list_of_faces[n]].v[i]];
      VERTEX          * const v =&mris->vertices         [mris->faces[list_of_faces[n]].v[i]];
      if (v->marked==0) {
        nv++;
        v->marked=1;
        if (vt->e) free(vt->e);
        vt->e = (int*)calloc(vt->vnum,sizeof(int));
      };
    }

  //mark and count the faces
  for (int n = 0 ; n < mris->nfaces ; n++)
    mris->faces[n].marked=0;
  nf=ne=0;
  for (int n = 0 ; n < nfs ; n++) {
    FACE *face= &mris->faces[list_of_faces[n]];
    if (face->marked==0) {
      nf++;
      face->marked=1;
      int vn0,vn1;
      for (int i = 0 ; i < 3 ; i++) {
        vn0 = face->v[i];
        if (i==2) vn1 = face->v[0];
        else vn1 = face->v[i+1];
        VERTEX_TOPOLOGY * const v0 = &mris->vertices_topology[vn0];
        VERTEX_TOPOLOGY * const v1 = &mris->vertices_topology[vn1];
        //edge vn0 <--> vn1 ?
        for (int p = 0 ; p < v0->vnum ; p++) {
          if (v0->v[p] == vn1 && v0->e[p]==0) {
            ne++;
            v0->e[p]=1;
            //mark the other edge
            for (int m = 0 ; m < v1->vnum ; m++) {
              if (v1->v[m]==vn0) {
                v1->e[m]=1;
                break;
              }
            }
            break;
          }
        }
      }
    }
  }

  //free the edge structure - reset marks to 0
  for (int n = 0 ; n < nfs ; n++) {
    FACE *face= &mris->faces[list_of_faces[n]];
    face->marked = 0;
    for (int i = 0 ; i < 3 ; i++) {
      VERTEX_TOPOLOGY * const vt=&mris->vertices_topology[face->v[i]];
      VERTEX          * const v =&mris->vertices         [face->v[i]];
      if (v->marked) {
        v->marked=0;
        if (vt->e) free(vt->e);
        vt->e=0;
      }
    }
  }

  // fprintf(WHICH_OUTPUT,"(%d,%d,%d)\n",nv,ne,nf);

  return (nv-ne+nf);
}

MRIP* MRIPextractFromMRIS(MRIS *mris, int defect_number) 
{
  int *list_of_faces, nfaces;

  /* count the number of faces and vertices */
  //counting the number of vertices
  int nvertices = 0 ;
  for (int n = 0 ; n < mris->nvertices ; n++)
    if (mris->vertices[n].marked2 == defect_number)
      nvertices++;
  if (nvertices==0) return NULL;
  //counting the number of faces
  nfaces = 0 ;
  for (int n = 0 ; n < mris->nfaces ; n++) {
    bool is_face = true;
    FACE *face=&mris->faces[n];
    face->marked=0;
    for (int i = 0 ; i < 3 ; i++)
      if (mris->vertices[face->v[i]].marked2 != defect_number) {
        is_face = false;
        break;
      };
    if (is_face) {
      nfaces++;
      face->marked=1;
    }
  }
  list_of_faces = new int[nfaces];
  nfaces=0;
  for (int n = 0 ; n < mris->nfaces ; n++)
    if (mris->faces[n].marked) list_of_faces[nfaces++]=n;

  /* first compute euler number of defect : if one return  0 */
  int euler = MRISgetEulerNumber(mris, list_of_faces, nfaces);

  //checking the value of the euler number
  if (euler==1) {
    delete [] list_of_faces;
    return NULL;
  }
  if ((euler+1)%2) {
    fprintf(WHICH_OUTPUT,"\n\n Surface non Valid %d\n\n",euler);
    delete [] list_of_faces;
    return NULL;
  }

  // number of loops to be cut
  int ncorrections = (1-euler)/2;

  /* allocate the patch with extra space */
  int max_vertices,max_faces;
  // _OverAlloc(2*loop.npoints+2*pdisk->disk.nvertices,2*(2*loop.npoints+pdisk->disk.nfaces+pdisk->ring.npoints));
  max_vertices = nvertices + ncorrections*(2*nvertices + MAX_EXTRA_VERTICES);
  max_faces = nfaces + ncorrections*(4*nvertices+ MAX_EXTRA_FACES);

  MRIP *mrip = MRIPalloc(max_vertices, max_faces);
  mrip->n_vertices = nvertices; // number of vertices before correction
  mrip->n_faces = nfaces; //number of faces before correction
  mrip->mris_source = mris;

  MRIS *mris_dst = mrip->mris;
  /* copy the necessary information */
  MRISreallocVerticesAndFaces(mris_dst, nvertices, nfaces);

  // the corresponding tables
  int *vt_to,*vt_from;
  vt_to = new int[max_vertices];
  for (int n = 0 ; n < max_vertices ; n++)
    vt_to[n]=-1;
  vt_from = new int[mris->nvertices];

  nvertices=0;
  for (int n = 0 ; n < mris->nvertices ; n++) {
    vt_from[n]=-1;
    if (mris->vertices[n].marked2 == defect_number) {
      VERTEX *vsrc = &mris->vertices[n];
      vt_from[n]=nvertices;
      vt_to[nvertices]=n;
      // copy the strict necessary
      MRISsetXYZ(mris_dst, nvertices, vsrc->x, vsrc->y, vsrc->z);
      nvertices++;
    }
  }
  mrip->vtrans_to = vt_to;
  mrip->vtrans_from = vt_from;

  //the corresponding tables
  int *ft_to,*ft_from;
  ft_to = new int[max_faces];
  for (int n = 0 ; n < max_faces ; n++)
    ft_to[n]=-1;
  ft_from = new int[mris->nfaces];

  nfaces=0;
  for (int n = 0 ; n < mris->nfaces ; n++) {
    ft_from[n]=-1;
    bool is_face = true;
    for (int i = 0 ; i < 3 ; i++)
      if (mris->vertices[mris->faces[n].v[i]].marked2 != defect_number) {
        is_face = false;
        break;
      };
    if (is_face) {
      FACE *fsrc = &mris->faces[n];
      FACE *f = &mris_dst->faces[nfaces];
      ft_from[n]=nfaces;
      ft_to[nfaces++]=n;
      for (int i = 0 ; i < 3 ; i++)
        f->v[i] = vt_from[fsrc->v[i]];
    }
  }
  mrip->ftrans_to = ft_to;
  mrip->ftrans_from = ft_from;

  //init the rest
  MRISinitSurface(mris_dst);

  return mrip;
}

void MRISinitSurface(MRIS *mris) 
{
  for (int n = 0 ; n < mris->nvertices ; n++) {
    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[n];
    VERTEX          * const v  = &mris->vertices         [n];
    vt->num = 0;
    v->marked = 0;
    if (vt->f) free(vt->f);
    if (vt->n) free(vt->n);
    if (vt->v) free(vt->v);
    vt->f = NULL;
    vt->n = NULL;
    vt->v = NULL;
  }

  // counting the number of faces per vertex
  for (int n = 0 ; n < mris->nfaces ; n++)
    for (int i = 0 ; i < 3 ; i++)
      mris->vertices_topology[mris->faces[n].v[i]].num++;

  // allocate the list of faces
  for (int n = 0 ; n < mris->nvertices ; n++) {
    VERTEX_TOPOLOGY * const v = &mris->vertices_topology[n];
    v->f = (int   *)calloc(mris->vertices_topology[n].num,sizeof(int));
    v->n = (uchar *)calloc(mris->vertices_topology[n].num,sizeof(uchar));
    v->num = 0;
  }

  // initialize the list of faces
  for (int n = 0 ; n < mris->nfaces ; n++) {
    for (int i = 0 ; i < 3 ; i++) {
      int vno = mris->faces[n].v[i];
      VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vno];
      v->f[v->num] = n;
      v->n[v->num++] = i;
    }
  }


  // counting the list of vertices
  for (int n = 0 ; n < mris->nvertices ; n++) {
    VERTEX_TOPOLOGY * const v = &mris->vertices_topology[n];
    clearVnum(mris,n);
    for (int p = 0 ; p < v->num ; p++) {
      FACE *face = &mris->faces[v->f[p]];
      for ( int i = 0 ; i < 3 ; i++) {
        int vn=face->v[i];
        if (vn==n) continue;
        if (mris->vertices[vn].marked) continue;
        mris->vertices[vn].marked=1;
        vnumAdd(mris,n,1);
      }
    }
    // allocate the list of vertices
    v->v = (int*)calloc(mris->vertices_topology[n].vnum,sizeof(int));
    clearVnum(mris,n);
    for (int p = 0 ; p < v->num ; p++) {
      FACE *face = &mris->faces[v->f[p]];
      for ( int i = 0 ; i < 3 ; i++) {
        int vn = face->v[i];
        if (vn == n) continue;
        if (mris->vertices[vn].marked == 0 ) continue;
        mris->vertices[vn].marked = 0;
        v->v[vnumAdd(mris,n,1)] = vn;
      }
    }
    v->v2num = v->vnum;
    v->v3num = v->vnum;
    v->vtotal = v->vnum;
  }

#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
#endif
  mrisCheckVertexFaceTopology(mris);
#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic pop
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic pop
#endif
}


MRIP *MRIPalloc(int nvertices, int nfaces) 
{
  MRIP *patch;

  /* allocate the patch */
  patch = (MRIP*)calloc(1,sizeof(MRIP));
  /* allocate the surface */
  MRIS *mris = MRISoverAlloc(nvertices,nfaces,nvertices,nfaces);

  patch->mris=mris;

  return patch;
}

void MRIPfree(MRIP **mrip) {
  MRIP *patch;

  patch = *mrip;
  *mrip = NULL;
  MRIS *mris=patch->mris;

  MRISfree(&mris);

  if (patch->vtrans_to) free(patch->vtrans_to);
  if (patch->ftrans_to) free(patch->ftrans_to);
  if (patch->vtrans_from) free(patch->vtrans_from);
  if (patch->ftrans_from) free(patch->ftrans_from);

  free(patch);

}

MRIP *MRIPclone(MRIP *src) 
{
  MRIP *dst;

  /* allocate the patch */
  dst = (MRIP*)calloc(1,sizeof(MRIP));
  /* clone the surface */
  dst->mris = MRISclone(src->mris);

  // do not clone the rest
  return dst;
}

Surface *MRIStoSurface(MRIS *mris) 
{
  //allocation of the surface
  Surface *surface=new Surface(mris->max_vertices,mris->max_faces);
  surface->nvertices = mris->nvertices;
  surface->nfaces = mris->nfaces;

  for (int n = 0 ; n < mris->nvertices;n++) {
    Vertex *vdst = &surface->vertices[n];
    VERTEX *vsrc = &mris->vertices[n];
    vdst->x=vsrc->x;
    vdst->y=vsrc->y;
    vdst->z=vsrc->z;
  }

  for (int n = 0 ; n < mris->nfaces ;n++) {
    Face *fdst = &surface->faces[n];
    FACE *fsrc = &mris->faces[n];
    fdst->v[0]=fsrc->v[0];
    fdst->v[1]=fsrc->v[1];
    fdst->v[2]=fsrc->v[2];
  }

  surface->InitSurface();

  return surface;
}

void MRIScopyHeader(MRIS *mris_src,MRIS *mris_dst) 
{
  mris_dst->type = mris_src->type;  // missing
  mris_dst->hemisphere = mris_src->hemisphere ;
  mris_dst->xctr = mris_src->xctr ;
  mris_dst->yctr = mris_src->yctr ;
  mris_dst->zctr = mris_src->zctr ;
  mris_dst->xlo = mris_src->xlo ;
  mris_dst->ylo = mris_src->ylo ;
  mris_dst->zlo = mris_src->zlo ;
  mris_dst->xhi = mris_src->xhi ;
  mris_dst->yhi = mris_src->yhi ;
  mris_dst->zhi = mris_src->zhi ;
  mris_dst->min_curv = mris_src->min_curv ;
  mris_dst->max_curv = mris_src->max_curv ;
  mris_dst->total_area = mris_src->total_area ;
  mris_dst->orig_area = mris_src->orig_area ;

  mris_dst->radius = mris_src->radius; // to be checked

  // just copy the pointer ///////////////////////////////////
#if 0
  mris_dst->linear_transform = mris_src->linear_transform ;
  mris_dst->inverse_linear_transform = mris_src->inverse_linear_transform ;
#endif
  mris_dst->lta = mris_src->lta;
  mris_dst->SRASToTalSRAS_ = mris_src->SRASToTalSRAS_;
  mris_dst->TalSRASToSRAS_ = mris_src->TalSRASToSRAS_;
  mris_dst->free_transform = 0 ;  // mark not to try to free them
  /////////////////////////////////////////////////////////////
  if (mris_src->v_frontal_pole)
    mris_dst->v_frontal_pole =
      &mris_dst->vertices[mris_src->v_frontal_pole - mris_src->vertices] ;
  if (mris_src->v_occipital_pole)
    mris_dst->v_occipital_pole =
      &mris_dst->vertices[mris_src->v_occipital_pole - mris_src->vertices] ;
  if (mris_src->v_temporal_pole)
    mris_dst->v_temporal_pole =
      &mris_dst->vertices[mris_src->v_temporal_pole - mris_src->vertices] ;

  // copy geometry info
  copyVolGeom(&mris_src->vg, &mris_dst->vg);
}

static MRIS* SurfaceToMRISwkr_new(Surface *surface) 
{
  MRIS *mris =
      MRISoverAlloc(surface->maxvertices,
                       surface->maxfaces,
                       surface->nvertices,
                       surface->nfaces);

  for (int n = 0 ; n < surface->nvertices ; n++) {
    Vertex const * const vsrc = &surface->vertices[n];
    
    MRISsetXYZ(mris,n,
      vsrc->x,
      vsrc->y,
      vsrc->z);
  }
  
  for (int n = 0 ; n < surface->nfaces ; n++) {
    Face const * const fsrc = &surface->faces[n];
    mrisAttachFaceToVertices(mris, n, fsrc->v[0], fsrc->v[1], fsrc->v[2]);
  }

  mrisCompleteTopology(mris);
  
  return mris;
}

static MRIS* SurfaceToMRISwkr_old(Surface *surface) 
{
  MRIS *mris =
      MRISoverAlloc(surface->maxvertices,
                       surface->maxfaces,
                       surface->nvertices,
                       surface->nfaces);

  //Vertices
  for (int n = 0 ; n < surface->nvertices ; n++) {
    Vertex const * const vsrc = &surface->vertices[n];
    
    VERTEX_TOPOLOGY * const vdstt = &mris->vertices_topology[n];
    VERTEX          * const vdst  = &mris->vertices         [n];
    
    MRISsetXYZ(mris,n,
      vsrc->x,
      vsrc->y,
      vsrc->z);
    
    vdst->ripflag = 0;
    vdst->marked  = 0;
    //vertices
    freeAndNULL(vdstt->v);
    vdstt->v = (int*)calloc(vsrc->vnum , sizeof(int));
    for (int p = 0 ; p < vsrc->vnum ; p++)
      vdstt->v[p]=vsrc->v[p];
      
    setVnum(mris,n,vsrc->vnum);
    vdstt->v2num  = vsrc->vnum;
    vdstt->v3num  = vsrc->vnum;
    vdstt->vtotal = vsrc->vnum;
    
    //faces
    if (vdstt->f) free(vdstt->f);
    vdstt->f=NULL;
    
    if (vdstt->n) free(vdstt->n);
    vdstt->n=NULL;
    
    vdstt->f = (int*)calloc(vsrc->fnum , sizeof(int));
    vdstt->n = (uchar*)calloc(vsrc->fnum , sizeof(uchar));
    for (int p = 0 ; p < vsrc->fnum ; p++) {
      vdstt->f[p]=vsrc->f[p];
      vdstt->n[p]=(uchar)vsrc->n[p];
    }
    vdstt->num=vsrc->fnum;
  }
  
  //Faces
  for (int n = 0 ; n < surface->nfaces ; n++) {
    FACE *fdst=&mris->faces[n];
    Face *fsrc = &surface->faces[n];
    fdst->v[0]=fsrc->v[0];
    fdst->v[1]=fsrc->v[1];
    fdst->v[2]=fsrc->v[2];
  }

#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
#endif
  mrisCheckVertexFaceTopology(mris);
#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic pop
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic pop
#endif
  
  return mris;
}

MRIS* SurfaceToMRIS(Surface *surface) {
  //
  // The possible difference in the order of faces around the vertices may cause slight differences
  // hence the old code is kept for comparison purposes
  //
  MRIS* mris =
    false 
    ? SurfaceToMRISwkr_old(surface) 
    : SurfaceToMRISwkr_new(surface);
    
  MRIScomputeNormals(mris);
  MRIScomputeTriangleProperties(mris);
  MRISsaveVertexPositions(mris,ORIGINAL_VERTICES);

  return mris;
}

bool MRISaddMRIP(MRIS *mris_dst, MRIP *mrip) 
{
  MRIS * const mris = mrip->mris;

  int const nvertices = mrip->n_vertices;
  int const nfaces    = mrip->n_faces;
  
  int * const vto     = mrip->vtrans_to;
  int * const fto     = mrip->ftrans_to;

  int currNumVertices = mris_dst->nvertices;
  int currNumFaces    = mris_dst->nfaces;
  int newNumVertices  = currNumVertices + mris->nvertices - nvertices;
  int newNumFaces     = currNumFaces    + mris->nfaces    - nfaces;

  MRISreallocVerticesAndFaces(mris_dst, newNumVertices, newNumFaces);

  cheapAssert(mris_dst->origxyz_status == mris->origxyz_status);
  

  //vertices
  for (int n = 0 ; n < mris->nvertices ; n++) {
    VERTEX const * const vsrc = &mris->vertices[n];
    
    if (n >= nvertices) vto[n] = currNumVertices++;
    int const vno_dst = vto[n];
    
    MRISsetXYZ(mris_dst, vno_dst,
      vsrc->x,
      vsrc->y,
      vsrc->z);
    
    MRISsetOriginalXYZ(mris_dst, vno_dst, 
      vsrc->x, 
      vsrc->y, 
      vsrc->z);
  }

  //faces
  for (int n = 0 ; n < mris->nfaces ; n++) {
    FACE const * const fsrc = &mris->faces[n];
    
    if (n >= nfaces) fto[n] = currNumFaces++;
    
    FACE * const fdst = &mris_dst->faces[fto[n]];
    
    //vertex indices
    fdst->v[0] = vto[fsrc->v[0]];
    fdst->v[1] = vto[fsrc->v[1]];
    fdst->v[2] = vto[fsrc->v[2]];
  }

  //stuff in vertices and faces
  MRISinitSurface(mris_dst);

  return true;
}
