//
// for calculating geodesics on a polyhedral surface
//

#include <stdlib.h>
#include <algorithm>  
#include <iomanip>
#include <iostream>
#include <map>
#include <stack>
#include <vector>

#include "geodesics.h"

#include "macros.h"
#include "mrisurf.h"
#include "timer.h"

#ifdef _OPENMP
#include "romp_support.h"
#endif

// Vertex
struct Vertex
{
  float x, y;
  int id;
};

// StackItem
struct StackItem
{
  Vertex a, b, c;
  int idx;
  float mina, maxa;
};

// Triangle
struct Triangle
{
  float length[3];
  float angle[3];
  int vert[3];
  int neighbor[3];
  bool inChain;
};

static int getIndex(int *arr, int vid);
static float distanceBetween(int v1, int v2, MRIS *surf);
static int findNeighbor(int faceidx, int v1, int v2, MRIS *surf);
static Vertex extendedPoint(Vertex A, Vertex B, float dA, float dB, float dAB);
static std::pair< int, int > makeKey(int a, int b);
static void progressBar(float progress);

Geodesics *computeGeodesics(MRIS *surf, float maxdist)
{
  int msec;
  Timer mytimer;
  printf("computeGeodesics(): maxdist = %g, nvertices = %d\n", maxdist, surf->nvertices);
  fflush(stdout);

  // pre-compute and set-up required values to build triangle chain:
  FACE *face;
  Triangle *triangle;
  std::vector< Triangle > triangles(surf->nfaces);
  int idx1, idx2;
  for (int nf = 0; nf < surf->nfaces; nf++) {
    face = &surf->faces[nf];
    triangle = &triangles[nf];
    for (int ns = 0; ns < 3; ns++) {
      idx1 = (ns + 1) % 3;
      idx2 = (ns + 2) % 3;
      triangle->length[ns] = distanceBetween(face->v[idx1], face->v[idx2], surf);
      triangle->neighbor[ns] = findNeighbor(nf, face->v[idx1], face->v[idx2], surf);
      triangle->vert[ns] = face->v[ns];
      triangle->angle[ns] = face->angle[ns];
    }
    triangle->inChain = false;
  }

  msec = mytimer.milliseconds();
  printf("precompute t = %g min\n", msec / (1000.0 * 60));
  fflush(stdout);

  std::cout << "computing geodesics within distance of " << maxdist << " mm\n";
  fflush(stdout);
  std::vector< std::vector< int > > nearestverts(surf->nvertices);
  int idxlookup[] = {0, 2, 1, 0};  // fast lookup table to find remaining index
                                   // can be removed... there's an easier way
  std::map< std::pair< int, int >, float > pathmap;
  // private vars:
  std::vector< int > chain;
  std::stack< StackItem > stack;
  std::vector< bool > isnearest(surf->nvertices);
  std::map< std::pair< int, int >, float >::iterator edge;
  StackItem stackitem;
  Vertex A, B, C, D;
  int iA, iB, iC, iD;
  int current_idx;
  float min_angle, max_angle, current_angle, distance;

  // ------ STEP 1 ------
  // compute each geodesic using the LOS algorithm. each geodesics is added
  // to the pathmap, but this will not account for every path.
  for (int vertexID = 0; vertexID < surf->nvertices; vertexID++) {
    VERTEX_TOPOLOGY const * const basevertex = &surf->vertices_topology[vertexID];
    // begin chain with each face that neighbors the current base vertex:
    for (int i = 0; i < basevertex->num; i++) {
      // clear triangle chain and surrounding vertices index:
      for (unsigned int c = 0; c < chain.size(); c++) triangles[chain[c]].inChain = false;
      chain.clear();
      for (unsigned int c = 0; c < isnearest.size(); c++) isnearest[c] = false;
      // set up initial triangle in plane:
      current_idx = basevertex->f[i];
      triangle = &triangles[current_idx];
      // formally add to chain:
      chain.push_back(current_idx);
      triangle->inChain = true;
      // get vertex indices in relation to their
      // placement in the face's vertex array:
      iC = getIndex(triangle->vert, vertexID);
      iA = (iC + 1) % 3;
      iB = (iC + 2) % 3;
      // compute min and max fov angles:
      min_angle = 0.0;
      max_angle = triangle->angle[iC];
      // compute vertex A along x axis:
      A.x = triangle->length[iB];
      A.y = 0.0;
      A.id = triangle->vert[iA];
      // compute vertex B in positive y (no need for this to be pre-computed):
      B.x = triangle->length[iA] * cos(max_angle);
      B.y = triangle->length[iA] * sin(max_angle);
      B.id = triangle->vert[iB];
      // reset vertex C (base vertex which represents the origin):
      C.x = 0.0;
      C.y = 0.0;
      C.id = triangle->vert[iC];
      // formally consider the distance from C to A as a geodesic:
      nearestverts[vertexID].push_back(A.id);
      isnearest[A.id] = true;
      edge = pathmap.find(makeKey(vertexID, A.id));
      if (edge != pathmap.end())
        edge->second = triangle->length[iB];
      else
        pathmap.insert(make_pair(makeKey(vertexID, A.id), triangle->length[iB]));
      // formally consider the distance from C to B as a geodesic:
      nearestverts[vertexID].push_back(B.id);
      isnearest[B.id] = true;
      edge = pathmap.find(makeKey(vertexID, B.id));
      if (edge != pathmap.end())
        edge->second = triangle->length[iA];
      else
        pathmap.insert(make_pair(makeKey(vertexID, B.id), triangle->length[iA]));
      // chain initialiaztion complete. get next triangle
      // and begin building chain:
      current_idx = triangle->neighbor[iC];
      // ------ build triangle chain ------
      while (true) {
        // check if the current triangle is valid or if it
        // already exists in the chain:
        if ((current_idx < 0) || (triangles[current_idx].inChain)) {
          // move on to next base triangle if the stack is empty:
          if (stack.empty()) break;
          // if not, just revert to the last stack item:
          else {
            stackitem = stack.top();
            A = stackitem.a;
            B = stackitem.b;
            C = stackitem.c;
            min_angle = stackitem.mina;
            max_angle = stackitem.maxa;
            current_idx = stackitem.idx;
            triangle = &triangles[current_idx];
            // trim the chain back to the current triangle:
            while ((chain.back() != current_idx) && (chain.size() > 0)) {
              triangles[chain.back()].inChain = false;
              chain.pop_back();
            }
            stack.pop();
          }
        }
        // triangle is valid, so add it to the chain:
        else {
          triangle = &triangles[current_idx];
          chain.push_back(current_idx);
          triangle->inChain = true;
          // find appropriate vertex indices for new triangle:
          iA = getIndex(triangle->vert, A.id);  // this can be optimized
          iB = getIndex(triangle->vert, B.id);
          iD = idxlookup[iA + iB];
          // calculate the planar position of the extended vertex D:
          D = extendedPoint(A, B, triangle->length[iB], triangle->length[iA], triangle->length[iD]);
          D.id = triangle->vert[iD];
          // calculate the angle that the vector D makes with x-axis:
          current_angle = atan2(D.y, D.x);
          // now calculate the distance to the origin:
          distance = sqrt(D.x * D.x + D.y * D.y);
          if (distance > maxdist) {
            current_idx = -1;  // this forces the next triangle invalid
            continue;
          }
          // potentially use the geodesics matrix here...
          if (!isnearest[D.id]) {
            nearestverts[vertexID].push_back(D.id);
            isnearest[D.id] = true;
          }
          // check if angle is visible within the fov:
          if ((current_angle < min_angle)) {
            C = A;
            A = D;
          }
          else if ((current_angle > max_angle)) {
            C = B;
            B = D;
          }
          else if (((current_angle <= max_angle) && (current_angle >= min_angle))) {
            // add geodesic to path map if shorter than the previous distance:
            edge = pathmap.find(makeKey(vertexID, D.id));
            if (edge != pathmap.end()) {
              if ((distance < edge->second) || (edge->second < 0.0)) edge->second = distance;
            }
            else {
              pathmap.insert(make_pair(makeKey(vertexID, D.id), distance));
            }
            // push triangle to the stack:
            stackitem.a = A;
            stackitem.b = D;
            stackitem.c = B;
            stackitem.idx = current_idx;
            stackitem.mina = min_angle;
            stackitem.maxa = current_angle;
            stack.push(stackitem);
            C = A;
            A = D;
            min_angle = current_angle;
          }
          // this is used to find bugs within the surface (so far I've only
          // seen problems in the fsaverage surface)
          else {
            // std::cout << "nan\n";
            // fflush(stdout);
            current_idx = -1;
            continue;
          }
        }
        // get the next triangle and repeat:
        iC = getIndex(triangle->vert, C.id);
        current_idx = triangle->neighbor[iC];
      }
    }
    if (vertexID % 1000 == 0) progressBar((float)vertexID / surf->nvertices);
  }
  progressBar(1.0);
  std::cout << std::endl;
  msec = mytimer.milliseconds();
  printf("step 1 t = %g min\n", msec / (1000.0 * 60));
  fflush(stdout);

  std::cout << "computing shortest paths and non-geodesics\n";
  fflush(stdout);
  Geodesics *geo = (Geodesics *)calloc(surf->nvertices, sizeof(Geodesics));
  // ------ STEP 2 ------
  // compute the shortest paths between each vertex (within given limit)
  int vi, vj;
  for (int k = 0; k < surf->nvertices; k++) {
    // set up nearest vertices reference:
    for (unsigned int d = 0; d < nearestverts[k].size(); d++) {
      isnearest[nearestverts[k][d]] == true;
    }
    for (unsigned int i = 0; i < nearestverts[k].size(); i++) {
      vi = nearestverts[k][i];
      edge = pathmap.find(makeKey(k, vi));
      if ((vi == k) || (edge != pathmap.end())) continue;

      for (unsigned int j = 0; j < nearestverts[k].size(); j++) {
        vj = nearestverts[k][j];
        if ((vi == vj) || (vj == k)) continue;
        // retrieve distance from k to i:
        edge = pathmap.find(makeKey(k, vj));
        if (edge == pathmap.end()) continue;
        distance = edge->second;
        // retrieve distance from j to k:
        edge = pathmap.find(makeKey(vj, vi));
        if (edge == pathmap.end()) continue;
        distance += edge->second;
        // retrieve distance from i to j:
        edge = pathmap.find(makeKey(k, vi));
        if ((edge == pathmap.end()) && (distance < maxdist))
          pathmap.insert(make_pair(makeKey(k, vi), distance));
        else if (distance < edge->second)
          edge->second = distance;
      }
      // search for vertices that are within distance limits
      // but weren't discovered by the triangle chain
      edge = pathmap.find(makeKey(k, vi));
      if (edge != pathmap.end()) {
        for (int side = 0; side < surf->vertices_topology[vi].vnum; side++) {
          if (!isnearest[surf->vertices_topology[vi].v[side]]) {
            distance = edge->second + 0.5;
            if (distance < maxdist) {
              isnearest[surf->vertices_topology[vi].v[side]] = true;
              nearestverts[k].push_back(surf->vertices_topology[vi].v[side]);
            }
          }
        }
      }
    }
    geo[k].vnum = 0;
    // geo[k].v = (int*) calloc(nearestverts[k].size(), sizeof(int));
    // geo[k].dist = (float*) calloc(nearestverts[k].size(), sizeof(float));
    for (unsigned int d = 0; d < nearestverts[k].size(); d++) {
      edge = pathmap.find(makeKey(k, nearestverts[k][d]));
      if (edge != pathmap.end()) {
        geo[k].v[geo[k].vnum] = nearestverts[k][d];
        geo[k].dist[geo[k].vnum] = edge->second;
        geo[k].vnum += 1;
        if (geo[k].vnum > MAX_GEODESICS) {
          std::cerr << "error: too many neighbors, try a smaller max distance\n";
          fflush(stdout);
          exit(1);
        }
      }
      isnearest[nearestverts[k][d]] == false;
    }

    if (k % 100 == 0) {
      progressBar((float)k / surf->nvertices);
    }
  }
  progressBar(1.0);
  std::cout << std::endl;

  msec = mytimer.milliseconds();
  printf("t = %g min\n", msec / (1000.0 * 60));
  fflush(stdout);

  return geo;
}

void geodesicsWrite(Geodesics *geo, int nvertices, char *fname)
{
  int vtxno;
  FILE *fp;
  int msec;

  printf("geodesicsWrite(): uniquifying\n");
  Timer mytimer;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (vtxno = 0; vtxno < nvertices; vtxno++) geodesicsUniquify(&geo[vtxno]);
  msec = mytimer.milliseconds();
  printf(" uniquification took %g min\n", msec / (1000.0 * 60));

  fp = fopen(fname, "wb");
  fprintf(fp, "FreeSurferGeodesics\n");
  fprintf(fp, "%d\n", -1);
  fprintf(fp, "%d\n", nvertices);
  for (vtxno = 0; vtxno < nvertices; vtxno++) {
    fwrite(&geo[vtxno].vnum, sizeof(int), 1, fp);
    fwrite(geo[vtxno].v, sizeof(int), geo[vtxno].vnum, fp);
    fwrite(geo[vtxno].dist, sizeof(float), geo[vtxno].vnum, fp);
  }
  fclose(fp);
}

Geodesics *geodesicsRead(char *fname, int *pnvertices)
{
  int magic, nthvtx;
  char tmpstr[1000];
  FILE *fp;

  fp = fopen(fname, "rb");
  if (fscanf(fp, "%s", tmpstr) != 1) {
    printf("ERROR (%s): could not read file\n", fname);
  }
  if (strcmp(tmpstr, "FreeSurferGeodesics")) {
    fclose(fp);
    printf("ERROR: %s not a geodesics file\n", fname);
    return (NULL);
  }
  if (fscanf(fp, "%d", &magic) != 1) {
    printf("ERROR (%s): could not read file\n", fname);
  }
  if (magic != -1) {
    fclose(fp);
    printf("ERROR: %s wrong endian\n", fname);
    return (NULL);
  }
  if (fscanf(fp, "%d", pnvertices) != 1) {
    printf("ERROR (%s): could not read file\n", fname);
  }
  fgetc(fp);  // swallow the new line
  printf("    geodesicsRead(): %s nvertices = %d, magic = %d\n", fname, *pnvertices, magic);
  fflush(stdout);
  Geodesics *geo = (Geodesics *)calloc(*pnvertices, sizeof(Geodesics));
  for (nthvtx = 0; nthvtx < *pnvertices; nthvtx++) {
    if (nthvtx % (*pnvertices / 10) == 0) {
      printf("%2d%% ", (int)round(100 * (float)nthvtx / (*pnvertices)));
      fflush(stdout);
    }
    if(!fread(&geo[nthvtx].vnum, sizeof(int), 1, fp)){
      printf("ERROR: %s failed fread\n", fname);
    }
    if(!fread(geo[nthvtx].v, sizeof(int), geo[nthvtx].vnum, fp)){
      printf("ERROR: %s failed fread\n", fname);
    }
    if(!fread(geo[nthvtx].dist, sizeof(float), geo[nthvtx].vnum, fp)){
      printf("Error: %s failed fread\n", fname);
    }
  }
  printf("\n");
  fclose(fp);
  return (geo);
}

/*!
\fn int geodesicsUniquify(Geodesics *geod)
\brief Removes relicants from the v (and dist) lists; vnum is updated.
*/
int geodesicsUniquify(Geodesics *geod)
{
  int nthnbr, *vlist, nunique, k, *vuniq;
  float *dist;

  // make a copy of the original list
  vlist = (int *)calloc(sizeof(int), geod->vnum);
  memcpy(vlist, geod->v, geod->vnum * sizeof(int));
  dist = (float *)calloc(sizeof(float), geod->vnum);
  memcpy(dist, geod->dist, geod->vnum * sizeof(float));

  // get unique list of vertices
  vuniq = unqiue_int_list(geod->v, geod->vnum, &nunique);
  if (nunique == geod->vnum) {
    // nothing changed
    free(vuniq);
    free(vlist);
    free(dist);
    return (nunique);
  }

  // replace structure values with unqiue
  for (nthnbr = 0; nthnbr < nunique; nthnbr++) {
    for (k = 0; k < geod->vnum; k++) {
      if (vlist[k] == vuniq[nthnbr]) {
        geod->dist[nthnbr] = dist[k];
        geod->v[nthnbr] = vlist[k];
        break;
      }
    }
  }
  geod->vnum = nunique;

  free(vuniq);
  free(vlist);
  free(dist);
  return (nunique);
}

static int getIndex(int *arr, int vid)
{
  int idx = std::distance(arr, std::find(arr, arr + 3, vid));
  // this can be removed:
  if (idx > 2) {
    std::cerr << "ERROR: getIndex returned value of 3\n";
    exit(1);
  }
  return idx;
}

static float distanceBetween(int v1, int v2, MRIS *surf)
{
  VERTEX_TOPOLOGY const * const vt = &surf->vertices_topology[v1];
  VERTEX          const * const v  = &surf->vertices         [v1];
  int const *ns = vt->v;
  int idx = std::distance(ns, std::find(ns, ns + vt->vnum, v2));
  return v->dist[idx];
}

static int findNeighbor(int faceidx, int v1, int v2, MRIS *surf)
{
  VERTEX_TOPOLOGY const * const vt = &surf->vertices_topology[v1];
  for (int nf = 0; nf < vt->num; nf++) {
    int f = vt->f[nf];
    if (f == faceidx) continue;
    for (int nv = 0; nv < 3; nv++) {
      if (surf->faces[f].v[nv] == v2) return f;
    }
  }
  return -1;
}

static Vertex extendedPoint(Vertex A, Vertex B, float dA, float dB, float dAB)
{
  Vertex D;
  float a, h, px, py;
  a = (dA * dA - dB * dB + dAB * dAB) / (2 * dAB);
  h = sqrt(dA * dA - a * a);
  px = A.x + a * (B.x - A.x) / dAB;
  py = A.y + a * (B.y - A.y) / dAB;
  D.x = px + h * (B.y - A.y) / dAB;
  D.y = py - h * (B.x - A.x) / dAB;
  return D;
}

static std::pair< int, int > makeKey(int a, int b)
{
  if (a < b)
    return std::pair< int, int >(a, b);
  else
    return std::pair< int, int >(b, a);
}

static void progressBar(float progress)
{
  if (!isatty(fileno(stdout))) return;
  int barwidth = 25;
  std::cout << "  [";
  int pos = barwidth * progress;
  for (int i = 0; i < barwidth; i++) {
    if (i < pos)
      std::cout << "-";
    else
      std::cout << " ";
  }
  std::cout << "]  " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}


MRI *GeoSmooth(MRI *src, double fwhm, MRIS *surf, Geodesics *geod, MRI *volindex, MRI *out)
{
  int vtxno;
  double gvar, gstd, gf;

  if(out == NULL)  {
    out = MRIallocSequence(src->width, src->height, src->depth,
                            MRI_FLOAT, src->nframes);
    if (out==NULL){
      printf("ERROR: GeoSmooth: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(src,out);
    MRIcopyPulseParameters(src,out);
  }
  MRIScomputeMetricProperties(surf);

  gstd = fwhm/sqrt(log(256.0));
  gvar = gstd*gstd;
  // scale factor does not really matter but good check because kernel should sum to 1
  gf = sqrt(pow(2*M_PI,2)*gvar*gvar); // = 2*M_PI*gvar
  printf("fwhm %g gstd %g %g %g, nv=%d\n",fwhm,gstd,gvar,gf,surf->nvertices);

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) 
#endif
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    ROMP_PFLB_begin
    int nthnbr,nbrvtxno,frame;
    double ksum, *sum, d, vkern;
    Geodesics *vtxgeod;
    if(vtxno % 10000 == 0) printf("vtxno %d\n",vtxno);
    if(surf->vertices[vtxno].ripflag) continue;

    // Remove replicate neighbors that have the same volume vertex no
    if(volindex) vtxgeod = VtxVolPruneGeod(geod, vtxno, volindex);
    else         vtxgeod = &geod[vtxno];

    // Set up init using self
    //vkern = surf->vertices[vtxno].area/gf; // scale by the area
    vkern = 1/gf;
    ksum = vkern;
    sum = (double *) calloc(sizeof(double),src->nframes);
    for(frame = 0; frame < src->nframes; frame++)
      sum[frame] = (vkern*MRIgetVoxVal(src,vtxno,0,0,frame));

    for(nthnbr = 0 ; nthnbr < vtxgeod->vnum; nthnbr++) {
      nbrvtxno = vtxgeod->v[nthnbr];
      if(surf->vertices[nbrvtxno].ripflag) continue;
      d = vtxgeod->dist[nthnbr];
      //vkern = surf->vertices[nbrvtxno].area*exp(-(d*d)/(2*gvar))/gf; // scale by the area
      vkern = exp(-(d*d)/(2*gvar))/gf; 
      ksum += vkern;
      for(frame = 0; frame < src->nframes; frame++)
	sum[frame] += (vkern*MRIgetVoxVal(src,nbrvtxno,0,0,frame));
    }

    for(frame = 0; frame < src->nframes; frame++)
      MRIsetVoxVal(out,vtxno,0,0,frame,(sum[frame]/ksum));
    surf->vertices[vtxno].valbak = ksum;
    surf->vertices[vtxno].val2bak = vtxgeod->vnum;
    free(sum);
    if(volindex) free(vtxgeod);
    ROMP_PFLB_end
  } // vtxno
  ROMP_PF_end

  return(out);
}

int GeoCount(Geodesics *geod, int nvertices)
{
  int c=0,cmax=0,vtxno;

  for(vtxno = 0; vtxno < nvertices; vtxno++){
    c += geod[vtxno].vnum;
    if(cmax < geod[vtxno].vnum) cmax = geod[vtxno].vnum;
  }
  printf(" GeoCount %d %d\n",c,cmax);
  return(c);
}

int GeoDumpVertex(char *fname, Geodesics *geod, int vtxno)
{
  FILE *fp;
  int nthnbr;
  fp = fopen(fname,"w");
  for(nthnbr = 0; nthnbr < geod[vtxno].vnum; nthnbr++){
    fprintf(fp,"%5d %g\n",geod[vtxno].v[nthnbr],geod[vtxno].dist[nthnbr]);
  }
  fclose(fp);
  return(0);
}


int geodesicsWriteV2(Geodesics* geo, int nvertices, char* fname) 
{
  int vtxno,*vnum,*vlist,nth,nthnbr,nnbrstot;
  float *dist;
  FILE *fp;

  nnbrstot = 0;
  for(vtxno = 0; vtxno < nvertices; vtxno++) nnbrstot += geo[vtxno].vnum;
  printf(" GeoCount %d\n",nnbrstot);

  fp = fopen(fname, "wb");
  fprintf(fp,"FreeSurferGeodesics-V2\n");
  fprintf(fp,"%d\n",-1);
  fprintf(fp,"%d\n",nvertices);
  fprintf(fp,"%d\n",nnbrstot);

  // Pack the number of neighbors into an array, then fwrite
  vnum = (int *) calloc(sizeof(int),nvertices);
  for(vtxno = 0; vtxno < nvertices; vtxno++) vnum[vtxno] = geo[vtxno].vnum;
  fwrite(vnum,sizeof(int), nvertices, fp);
  free(vnum);

  // Pack the neighbor vertex numbers and dist into an arrays, then fwrite
  vlist = (int *)   calloc(sizeof(int),  nnbrstot);
  dist  = (float *) calloc(sizeof(float),nnbrstot);
  nth = 0;
  for(vtxno = 0; vtxno < nvertices; vtxno++){
    for(nthnbr = 0 ; nthnbr < geo[vtxno].vnum; nthnbr++){
      vlist[nth] = geo[vtxno].v[nthnbr];
      dist[nth]  = geo[vtxno].dist[nthnbr];
      nth ++;
    }
  }
  fwrite(vlist,sizeof(int), nnbrstot, fp);
  free(vlist);
  fwrite(dist,sizeof(float), nnbrstot, fp);
  free(dist);

  fclose(fp);
  return(0);
}

Geodesics* geodesicsReadV2(char* fname, int *pnvertices) 
{
  int magic;
  char tmpstr[1000];
  FILE *fp;
  int vtxno,*vnum,*vlist,nth,nthnbr,nnbrstot;
  float *dist;
  int msec; 

  fp = fopen(fname, "rb");
  fscanf(fp,"%s",tmpstr);
  if(strcmp(tmpstr,"FreeSurferGeodesics-V2")){
    fclose(fp);
    printf("   %s\n",tmpstr);
    printf("ERROR: %s not a geodesics file\n",fname);
    return(NULL);
  }
  fscanf(fp,"%d",&magic);
  if(magic != -1){
    fclose(fp);
    printf("ERROR: %s wrong endian\n",fname);
    return(NULL);
  }
  fscanf(fp,"%d",pnvertices);
  fscanf(fp,"%d",&nnbrstot);
  fgetc(fp); // swallow the new line
  printf("    geodesicsReadV2(): %s nvertices = %d, magic = %d\n",fname,*pnvertices,magic);fflush(stdout);

  printf("  alloc vnum \n");fflush(stdout);
  vnum = (int *) calloc(sizeof(int),*pnvertices);
  printf("  reading in vnum %d \n",*pnvertices);fflush(stdout);
  Timer mytimer;
  fread(vnum,sizeof(int), *pnvertices, fp);
  msec = mytimer.milliseconds() ;
  printf("  t = %g min\n",msec/(1000.0*60));

  mytimer.reset();
  printf("  allocing vlist and dlist %d\n",nnbrstot);fflush(stdout);
  vlist = (int *)   calloc(sizeof(int),  nnbrstot);
  dist  = (float *) calloc(sizeof(float),nnbrstot);
  printf("  reading in vlist %d\n",nnbrstot);fflush(stdout);
  fread(vlist,sizeof(int), nnbrstot, fp);
  msec = mytimer.milliseconds(); printf("  t = %g min\n",msec/(1000.0*60));
  printf("  reading in dist %d\n",nnbrstot);fflush(stdout);
  fread(dist, sizeof(float), nnbrstot, fp);
  msec = mytimer.milliseconds(); printf("  t = %g min\n",msec/(1000.0*60));
  printf("Alloc geo\n");
  Geodesics *geo = (Geodesics*) calloc(*pnvertices, sizeof(Geodesics));
  msec = mytimer.milliseconds(); printf("  t = %g min\n",msec/(1000.0*60));
  printf("Setting\n");
  nth = 0;
  for(vtxno = 0; vtxno < *pnvertices; vtxno++){
    geo[vtxno].vnum = vnum[vtxno];
    for(nthnbr = 0 ; nthnbr < geo[vtxno].vnum; nthnbr++){
      geo[vtxno].v[nthnbr]    = vlist[nth];
      geo[vtxno].dist[nthnbr] = dist[nth];
      nth ++;
    }
  }
  free(vlist);
  free(dist);
  free(vnum);
  msec = mytimer.milliseconds(); printf("  t = %g min\n",msec/(1000.0*60));

  fclose(fp);

  return(geo);
}

// distance along the sphere between two  vertices
double MRISsphereDist(MRIS *sphere, VERTEX *vtx1, VERTEX *vtx2)
{
  double r, dot, d;
  r = sqrt(vtx1->x * vtx1->x + vtx1->y * vtx1->y + vtx1->z * vtx1->z);
  dot = (vtx1->x * vtx2->x + vtx1->y * vtx2->y + vtx1->z * vtx2->z)/(r*r);
  d = acos(dot)*r;
  //printf("r = %g, dot = %g, d = %g\n",r,dot,d);
  return(d);
}

double geodesicsCheckSphereDist(MRIS *sphere, Geodesics *geod)
{
  int vtxno;
  double e;

  e = 0;
  for(vtxno = 0; vtxno < sphere->nvertices; vtxno++){
    int nbrvtxno,nthnbr;
    VERTEX *vtx1,*vtx2;
    double dA, dB;
    vtx1 = &sphere->vertices[vtxno];
    for(nthnbr = 0 ; nthnbr < geod[vtxno].vnum; nthnbr++) {
      nbrvtxno = geod[vtxno].v[nthnbr];
      dA = geod[vtxno].dist[nthnbr];
      vtx2 = &sphere->vertices[nbrvtxno];
      dB = MRISsphereDist(sphere, vtx1, vtx2);
      printf("%5d %5d %6.4f %6.4f  %8.6f %8.4f\n",vtxno,nbrvtxno,dA,dB,dA-dB,e);
      e += fabs(dA-dB);
    }
  }

  return(e);
}


int VtxVolIndexCompare(const void *a, const void *b)
{

  VTXVOLINDEX *vvi1;
  VTXVOLINDEX *vvi2;

  vvi1 = ((VTXVOLINDEX *) a);
  vvi2 = ((VTXVOLINDEX *) b);

  /* Sort by Volume Index */
  if(vvi1->volindex < vvi2->volindex) return(-1);
  else if(vvi1->volindex > vvi2->volindex) return(+1);
  else{
    if(vvi1->vtxno < vvi2->vtxno) return(-1);
    return(+1);
  }
}

int VtxVolIndexSort(VTXVOLINDEX *vvi, int nlist)
{
  qsort(vvi,nlist,sizeof(VTXVOLINDEX),VtxVolIndexCompare);
  return(0);
}

int VtxVolIndexPrint(VTXVOLINDEX *vvi, int nlist)
{
  int n;
  for(n=0; n < nlist; n++)
    printf("%3d %5d %7d  %6.4f\n",n,vvi[n].vtxno,vvi[n].volindex,vvi[n].dist);
  return(0);
}

int VtxVolIndexSortTest(int nlist)
{
  VTXVOLINDEX *vvi,*vvi2;
  int n,nunique;
  vvi = (VTXVOLINDEX *) calloc(sizeof(VTXVOLINDEX),nlist);
  for(n=0; n < nlist; n++){
    vvi[n].vtxno = round(100*drand48());
    vvi[n].volindex = round(100*drand48());
    vvi[n].dist = 10*drand48();
  }
  vvi[0].volindex = vvi[nlist-1].volindex;
  vvi[1].volindex = vvi[nlist-1].volindex;
  vvi[2].volindex = vvi[nlist-2].volindex;
  
  VtxVolIndexPrint(vvi, nlist);
  VtxVolIndexSort(vvi, nlist);
  printf("-----------------\n");
  VtxVolIndexPrint(vvi, nlist);

  vvi2 = VtxVolIndexUnique(vvi, nlist, &nunique);
  printf("%d -----------------\n", nunique);
  VtxVolIndexPrint(vvi2, nunique);

  return(0);
}

VTXVOLINDEX *VtxVolIndexUnique(VTXVOLINDEX *vvi, int nlist, int *nunique)
{
  int n,nthu,nhits;
  VTXVOLINDEX *uvvi;
  double distsum;

  VtxVolIndexSort(vvi, nlist);

  uvvi = (VTXVOLINDEX *) calloc(sizeof(VTXVOLINDEX),nlist);

  nthu = 0;
  n = 0;
  memcpy(&uvvi[nthu],&vvi[n],sizeof(VTXVOLINDEX));
  distsum = vvi[n].dist;
  nhits = 1;
  for(n=1; n<nlist; n++) {
    if(uvvi[nthu].volindex != vvi[n].volindex) {
      uvvi[nthu].dist = distsum/nhits;
      nthu ++;
      memcpy(&uvvi[nthu],&vvi[n],sizeof(VTXVOLINDEX));
      distsum = vvi[n].dist;
      nhits = 1;
    }
    else {
      distsum += vvi[n].dist;
      nhits ++;
    }
  }
  *nunique = nthu+1;
  return(uvvi);
}

VTXVOLINDEX *VtxVolIndexPack(Geodesics *geod, int vtxno, MRI *volindex)
{
  int nthnbr,nbrvtxno;
  VTXVOLINDEX *vvi;
  //Geodesics *geodvtx = &geod[vtxno];

  vvi = (VTXVOLINDEX *) calloc(sizeof(VTXVOLINDEX),geod[vtxno].vnum);

  for(nthnbr = 0 ; nthnbr < geod[vtxno].vnum; nthnbr++) {
    nbrvtxno = geod[vtxno].v[nthnbr];
    vvi[nthnbr].vtxno = nbrvtxno;
    vvi[nthnbr].dist = geod[vtxno].dist[nthnbr];
    vvi[nthnbr].volindex = MRIgetVoxVal(volindex,nbrvtxno,0,0,0);
  }
  //vvi[nthnbr].vtxno = vtxno;
  //vvi[nthnbr].dist = 0;


  return(vvi);
}

Geodesics *VtxVolPruneGeod(Geodesics *geod, int vtxno, MRI *volindex)
{
  VTXVOLINDEX *vvi, *uvvi; 
  Geodesics *ugeod;
  int nunique;
  int nthnbr;

  vvi = VtxVolIndexPack(geod, vtxno, volindex);
  //if(vtxno == 123093)  VtxVolIndexPrint(vvi, geod[vtxno].vnum);

  uvvi = VtxVolIndexUnique(vvi, geod[vtxno].vnum, &nunique);

  ugeod = (Geodesics *) calloc(sizeof(Geodesics),1);
  ugeod->vnum = nunique;
  for(nthnbr = 0 ; nthnbr < ugeod->vnum; nthnbr++) {
    ugeod->v[nthnbr] = uvvi[nthnbr].vtxno;
    ugeod->dist[nthnbr] = uvvi[nthnbr].dist;
  }

  free(vvi);
  free(uvvi);

  return(ugeod);
}

