#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <array>

#include "ventfix.h"
#include "colortab.h"
#include "mri_identify.h"
#include "diag.h"

// fix ventricles that are underlabled in FreeSurfer
// 1. binarize input asegps using threshmin as threshold
//    voxels not belonging to any clusters in asegps (voxels not labelled) will have value 1; 
//    otherwise, they will have value 0
// 2. use brainmask to cut out the binarize volume
// 3. create clusters from the binary volume (the areas not labelled in asegps).
// 4. output cluster list as MRI structure (ocn) with voxel values as assigned cluster numbers
// 5. call ExpandSegIndices() to find/label voxels in segid cluster but not labeled in asegps
MRI* VentFix::fixasegps(MRI *asegps, MRI *brainmask, char *segids, float threshmin, int niters, int nmax, int topo)
{
  // Prepare the output binary volume
  MRI *binVol = MRIallocSequence(asegps->width, asegps->height, asegps->depth, MRI_INT, 1);
  if (binVol == NULL) exit(1);
  MRIcopyHeader(asegps, binVol);

  // binarize input asegps
  // if a voxel value < threshmin (not belonging to any clusters in asegps), set its value to 1;
  // otherwise, set its value to 0
  MRIbinarize(asegps, binVol, threshmin, 1, 0);

  // apply brainmask to binary volume
  MRI *newbinVol = MRImask(binVol, brainmask, NULL, 0, 0);
  MRIfree(&binVol);
  binVol = newbinVol;

  // create clusters from the binary volume
  // these are areas not labelled in asegps
  int nClusters;
  VOLCLUSTER **ClusterList = clustGetClusters(binVol, 0, threshmin, 0, 0, 0, NULL, &nClusters, NULL);

  // output cluster list as MRI structure with voxel values as assigned cluster numbers
  MRI *ocnVol = clustClusterList2Vol(ClusterList, nClusters, binVol, 0, 0);

  // mgz doesn't have cluster information.
  // ??? do we need to do this if we are not outputing .lut ???
#if 1
  ocnVol->ct = CTABalloc(nClusters+1);
  strcpy(ocnVol->ct->entries[0]->name, "Unknown");

  int n;
  for (n = 0; n < nClusters; n++)
    sprintf(ocnVol->ct->entries[n+1]->name, "Cluster-%03d", n+1);
#endif

#if 0  // testing
  CTABwriteFileASCII(ocnVol->ct, "ocn.lut");
#endif

  // free binary volume
  MRIfree(&binVol);

  // free ClusterList
  clustFreeClusterList(&ClusterList, nClusters);


  int nc;
  fsPointSet centroid;
  MRI *segVol = NULL, *newsegVol = NULL;

  char *restsegids = (char*)segids;
  char *nextsegid = strtok_r(segids, ",", &restsegids);
  while (nextsegid != NULL)
  {
    if (segVol == NULL)
      segVol = asegps;
    else
      segVol = newsegVol;

    int segid = atoi(nextsegid);
    nextsegid = strtok_r(NULL, ",", &restsegids);
    printf("\nexpandSegIndices for segid %d\n", segid);
    newsegVol = ExpandSegIndices(segVol, segid, ocnVol, niters, nmax, topo, &nc, centroid, NULL);
    printf("segid %d niters %d nmax %d topo %d   nc %5d \n", segid, niters, nmax, topo, nc);

    // free MRI structures returned by previous ExpandSegIndices() call
    if (segVol != asegps)
      MRIfree(&segVol);
  }

  MRIfree(&ocnVol);

  return newsegVol;
}

//-----------------------------------------------------
// 1. copy MRI *seg (asegps) to MRI *newseg
// 2. find all voxels in asegps with segid, remember their (c, r, s) in a vector idx (first iteration)
// 3. go through the vector idx, 
//    for each voxel in the vector
//      look for its neighbor defined by topo
//      if the neighbor belongs to a cluster in MRI *ocn && 
//        it doesn't belong to any clusters in MRI *newseg (voxel = 0, this is the voxel underlabled when creating asegps), 
//          add the voxel (c, r, s) to be processed vector idxnew, assign it value segid in MRI *newseg
// 4. add members of vector idxnew to vector idx
// 5. repeat #3, #4 (???skip previous processed voxels???) until we find all neighbors, 
//    we have gone through niters iterations, or nmax voxels have been added to the cluster
MRI* VentFix::ExpandSegIndices(MRI *seg, int segid, MRI *ocn, int niters, int nmax, int topo, int *nexpansions, fsPointSet &centroid, MRI* newseg)
{
  if(seg == NULL) return(NULL);
  newseg = MRIcopy(seg,newseg);
  if(newseg == NULL) return(NULL);
  MRIcopyHeader(seg, newseg);
  MRIcopyPulseParameters(seg, newseg);

  std::vector<std::array<int,3>> idx;
  int ok = 1;
  int iter = 0;
  int nprev = 0;
  double csum=0, rsum=0, ssum=0;
  *nexpansions = 0;
  while(ok){
    int nchanges = 0;
    iter++;
    nprev= idx.size();
    if(iter == 1){
      // In the first iteration, just find the voxels in the segid and record their indices
      int c;
      for(c=0; c < seg->width; c++){
	int r, s;
	for(r=0; r < seg->height; r++){
	  for(s=0; s < seg->depth; s++){
	    float v = MRIgetVoxVal(seg,c,r,s,0);
	    if(v != segid) continue;
	    std::array<int,3> crs = {c,r,s};
	    idx.push_back(crs);
	  }
	}
      }
      nchanges = idx.size();
      printf("Iter %3d %6d %6d %6d\n",iter,(int)idx.size(),nchanges,*nexpansions);
      if(niters > 0 && iter == niters) break;
      continue;
    }
    // In the 2nd+ iterations, go through all the segid voxels found
    // in the previous iterations, and look at their nearest
    // neighbors. If the NN is 0 and in the ocn, then set to segid
    std::vector<std::array<int,3>> idxnew;
    for(int n=0; n < idx.size(); n++){
      int c0 = idx[n][0];
      int r0 = idx[n][1];
      int s0 = idx[n][2];

      // Go through the nearest neighbors
      int c,r,s,dc,dr,ds;
      for(dc = -1; dc <= 1; dc++){
	c = c0 + dc;
	if(c < 0 || c >= seg->width) continue;
	for(dr = -1; dr <= 1; dr++){
	  r = r0 + dr;
	  if(r < 0 || r >= seg->height) continue;
	  for(ds = -1; ds <= 1; ds++){
	    s = s0 + ds;
	    if(s < 0 || s >= seg->depth) continue;
	    int debug = 0;
	    if(c == 139 && r == 153 && s == 73) debug = 1;
	    // proceed or not if this voxel is in the topology contraint
	    // topo = 1 = face neighbors only
	    // topo = 2 = face neighbors + edge neighbors only
	    // topo = 3 = face neighbors + edge neighbors + corner neighbors
	    if((abs(dc)+abs(dr)+abs(ds)) > topo) {
	      if(debug) printf("#@# %d %d %d topo skip (%d %d %d)\n",c,r,s,c0,r0,s0);
	      fflush(stdout);
	      continue; 
	    }
	    if(ocn){
	      // Only add a new point if it is in the ocn
	      float m = MRIgetVoxVal(ocn,c,r,s,0);
	      // The "ocn" here is expected to be the OCN from the
	      // cluster analysis where the unfilled ventricle is
	      // expected to be its own cluster/island. Each cluster
	      // has its own number, sorted by cluster size. The first
	      // cluster is probably the entire background, so skip that
	      // by only considering cluster numbers > 1
	      if(debug) printf("#@# %d %d %d  m=%g\n",c,r,s,m);
	      fflush(stdout);
	      if(m < 1.5) continue; 
	    }
	    float v = MRIgetVoxVal(newseg,c,r,s,0);
	    // Only proceed if neighboring voxel=0
	    if(debug) printf("#@# %d %d %d  m=%g\n",c,r,s,v);
	    if(v != 0 || v == segid) continue;
	    if(debug) printf("#@# %d %d %d  adding\n",c,r,s);
	    std::array<int,3> crs = {c,r,s};
	    // Add this voxel to the list
	    idxnew.push_back(crs);
	    csum += c;
	    rsum += r;
	    ssum += s;
	    // Seg the voxel to be segid in the new output volume
	    MRIsetVoxVal(newseg,c,r,s,0,segid);
	  } //dc
	} //dr
      } //ds
    } // n
    // Count the number of changes
    // First, get a unique list of indices
    idx.insert(idx.end(),idxnew.begin(),idxnew.end());
    std::vector<std::array<int,3>>::iterator it;
    it = std::unique(idx.begin(), idx.end());
    idx.resize(std::distance(idx.begin(),it) );
    nchanges = idx.size()-nprev;
    *nexpansions += nchanges;
    printf("Iter %3d %6d %6d %6d\n",iter,(int)idx.size(),nchanges,*nexpansions);
    if(niters > 0 && iter == niters) break;
    if(nchanges == 0) break;
    if(nmax > 0 && *nexpansions > nmax) break;
  } // iters

  if(*nexpansions > 0){
    csum /= *nexpansions;
    rsum /= *nexpansions;
    ssum /= *nexpansions;
    
    MATRIX *crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[1][1] = csum;
    crs->rptr[2][1] = rsum;
    crs->rptr[3][1] = ssum;
    crs->rptr[4][1] = 1;
    MATRIX *vox2ras = MRIxfmCRS2XYZ(newseg, 0);
    MATRIX *ras = MatrixMultiply(vox2ras,crs,NULL);
    printf("Centriod CRS %6.4lf %6.4lf %6.4lf\n",csum,rsum,ssum);
    printf("Centriod RAS %6.4lf %6.4lf %6.4lf\n",
	   ras->rptr[1][1],ras->rptr[2][1],ras->rptr[3][1]);
    
    centroid.vox2ras = "scanner_ras";
    fsPointSet::Point centroidp;
    centroidp.x = ras->rptr[1][1];
    centroidp.y = ras->rptr[2][1];
    centroidp.z = ras->rptr[3][1];
    centroid.add(centroidp);
    
    MatrixFree(&crs);
    MatrixFree(&vox2ras);
    MatrixFree(&ras);
  }

  return(newseg);
}

