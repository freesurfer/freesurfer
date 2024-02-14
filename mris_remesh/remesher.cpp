#ifdef NDEBUG
#undef NDEBUG  // enable assert
#endif

#include <cassert>
#include <queue>
#include <utility>

#include "remesher.h"

using namespace std;

#define unused(x) ((void)(x))


Remesher::Remesher(const MRIS *surf)
{
  points3d.resize(surf->nvertices);
  for (unsigned int i = 0; i < surf->nvertices; i++) {
    VERTEX *vertex = &surf->vertices[i];
    points3d[i][0] = vertex->x;
    points3d[i][1] = vertex->y;
    points3d[i][2] = vertex->z;
  }

  tria.resize(surf->nfaces);
  for (unsigned int i = 0; i < surf->nfaces; i++) {
    FACE *face = &surf->faces[i];
    tria[i] = { face->v[0], face->v[1], face->v[2] };
  }

  init();
}


MRIS * Remesher::toSurface()
{
  MRIS * surf = MRISalloc(points3d.size(), tria.size());

  for (int n = 0 ; n < surf->nvertices ; n++) {
    Vector vtx = points3d[n];
    MRISsetXYZ(surf, n, vtx[0], vtx[1], vtx[2]);
    MRISsetOriginalXYZ(surf, n, vtx[0], vtx[1], vtx[2]);
  }

  setFaceAttachmentDeferred(surf, true);
  for (int n = 0 ; n < surf->nfaces ; n++) {
    std::vector<int> tri = tria[n];
    mrisAttachFaceToVertices(surf, n, tri[0], tri[1], tri[2]);
  }
  setFaceAttachmentDeferred(surf, false);

  mrisCompleteTopology(surf);
  MRIScomputeMetricProperties(surf);

  return surf;
}


void Remesher::init()
{
  rmFreeVertices(false);  // do not call init form there, uses only tria
  createVtoT();           // create vertex to tria structure
  createEdges();          // create edge structures
  createVBoundary();      // create vertex boundary structures and test if manifold
  // reduceComponents();
  // orientTria();
  createQualities();
}


/** Remove free vertices (not used in any triangle)
  Uses only tria structure!
 If fix=true then init is called.*/
int Remesher::rmFreeVertices(bool fix)
{
  if (verbose > 0) cout << "Remesher::rmFreeVertices(bool fix) ... " << flush;
  
  bool hasridge = (onridge.size() == points3d.size());
  
  // keep vertices used in tria
  // mark others for deletion (-1)
  vector < int > newvidx(points3d.size(),-1);
  for (unsigned int t=0;t<tria.size();t++)
  for (unsigned int i=0;i<3;i++)
  {
    newvidx[tria[t][i]] = 0;
  }
  
  // create newvidx (new vertex numbers, -1 means delete)
  int vcount = 0;
  for (unsigned int v = 0;v<points3d.size(); v++)
    if (newvidx[v] ==0)
    {
      newvidx[v] = vcount;
      vcount++;
    }
    
  
  //delete vertices
  int delvertex = points3d.size() - vcount;
  if (delvertex > 0)
  {
    vector < Vector > new_points3d(vcount);
    int vvcount = 0;
    for (unsigned int v=0 ; v< points3d.size() ; v++)
    {
      if (newvidx[v] > -1)
      {
        new_points3d[vvcount] = points3d[v];
        vvcount++;
      }
    }
    points3d = new_points3d;
    
    // update tria
    for (unsigned int t =0;t<tria.size();t++)
    for (unsigned int i=0;i<3;i++)
    {
      tria[t][i] = newvidx[tria[t][i]];
      assert(tria[t][i] != -1);
    }
    

    if (hasridge)
    {
       vector <bool> newridge (new_points3d.size());
       for (unsigned int i = 0;i<newvidx.size();i++)
         if (newvidx[i] != -1) newridge[newvidx[i]] = onridge[i];
       onridge = newridge;
    }    
    
    
    
     if (verbose > 0)  cout << delvertex << " vertices deleted! " << endl;
    
    
        
  }
  else if (verbose > 0) cout << " no free vertices! " << endl;
  
  if (fix) init(); // needs to be called in all cases because trias could have been deleted earlier
  
  
  return (delvertex);
}


/** Creates vertex to tria structures:<br>
    v_to_t contains the tria for each vertex<br>
    V_to_lv contains the local index.<br>
    Uses only tria.
*/
void Remesher::createVtoT(void)
{
  v_to_t.clear();
  v_to_t.resize(points3d.size());
  v_to_lv.clear();
  v_to_lv.resize(points3d.size());
  
  for (unsigned int t = 0; t<tria.size(); t++)
  for (unsigned int j = 0; j<3; j++)
  {
    bool found = false;
    unsigned int i = 0;
    int vertex = tria[t][j];
    while (!found && i<v_to_t[vertex].size())
    {
      found = (v_to_t[vertex][i] == (int)t);
      i++;
    }
    if (!found)
    {
       v_to_t[vertex].push_back(t);
       v_to_lv[vertex].push_back(j);
    }
  }
}


/** Create edge datastructures. 
  Every edge has two vertices (stored in e_to_v).
  Some vertex combinations have a connecting edge (stored +1 in matrix v_to_e).
  Every edge is in at least one triangle, all triangles are stored in e_to_t.
  Every triangle has exactly three edges (stored in t_to_e).*/
void Remesher::createEdges()
{
  if (verbose > 0) cout << "Remesher::createEdges() ... " << flush;
  e_to_t.clear();
  t_to_e.resize(tria.size());
  e_to_v.clear();
  std::vector< int > temp_et;
  std::vector< int > temp_ev(2);
  pointcount = points3d.size();
  v_to_e.redim(pointcount,pointcount);
  v_to_e.fillzero();
  int ecount=0;
  // got through trias
  for (unsigned int t = 0; t<tria.size(); t++)
  {
  //cout << endl << " working on tria : " << t << endl;
    t_to_e[t].resize(tria[t].size());
    if (tria[t].size() != 3)
    {
      cerr << " tria[ " << t << " ] : " << tria[t] << endl;
    }
    assert(tria[t].size() == 3);
    // go through (usually 3) edges
    for (unsigned int e = 0; e<tria[t].size(); e++)
    {
      int v1 = tria[t][e];
      int v2 = tria[t][(e+1)%3];
      if (v1> v2) { int temp = v1;  v1 = v2; v2 = temp;  }
      
      if (v1 < 0 || v2 >= pointcount)
      {
        cerr << " tria[ "<< t <<" ]" << tria[t] << endl;
        cerr << " pointcount: " << pointcount << " v1: " << v1 << " v2: " << v2 << endl;
        assert (v1 > -1);
        assert (v2 < pointcount);
      }
      
      // check if edge exist
      int edge = (int)v_to_e(v1, v2) - 1;
    
      if (edge == -1) // new edge
      {
        // save new edge number in Matrix
        v_to_e(v1, v2) = ecount+1;
        v_to_e(v2, v1) = ecount+1;
        edge = ecount;
        
        // Save the two vertices in edge list (for vertices)
        temp_ev[0] = v1;
        temp_ev[1] = v2;
        e_to_v.push_back(temp_ev);
        
        // Save this triangle in the edge list (for triangles)
        temp_et.resize(1);
        temp_et[0] = t;
        e_to_t.push_back(temp_et);
        
        // increase edge count
        ecount++;
      }
      else
      {
        // make sure, this tria is not in edge list yet
        for (unsigned int tt = 0;tt<e_to_t[edge].size();tt++)
        {
          if (e_to_t[edge][tt] == (int) t)
          {
            cerr << endl <<  "Tria " << t << " found in e_to_t[ " << edge << " ] : " << e_to_t[edge] << endl;
            cerr << " Maybe tria has double vertex? tria : " << tria[t] << endl;
          }
          assert(e_to_t[edge][tt] != (int)t);
        }
        // Add this triangle to edge list for triangles
        e_to_t[edge].push_back(t);
      }
      
      // store edge in t_to_e list for this tria
      t_to_e[t][e] = edge;
  
    }
  }
  if (verbose > 0) cout << " DONE" << endl;
}


/** Creates boundary information for vertices.
  onboundary is an array specifying if specific vertex is on boundary.
 indexwoboundary is index for only the inner vertices.
 innercount is the number of inner points.
 It is also determined, if we have a manifold (each edge has 1 or two triangles).*/
void Remesher::createVBoundary()
{
  if (verbose > 0) cout << "Remesher::createVBoundary() ... " << flush;
  unsigned int esize = e_to_t.size();
  assert (esize > 0);
  assert (esize == e_to_v.size());
  
  onboundary.clear();
  bedgecount = 0;
  unsigned int psize = points3d.size();
  onboundary.resize(psize,false);
  int maxt = 0;
  int medg = -1;
  for (unsigned int e = 0; e< esize; e++)
  {
    int tnum = e_to_t[e].size();
    if (tnum > maxt)
    {
       maxt=tnum;
       medg=e;
    }
    if (tnum == 1)
    {
      onboundary[e_to_v[e][0]] = true;
      onboundary[e_to_v[e][1]] = true;
      bedgecount++;    
    }
  
  }
  
  indexwoboundary.resize(psize);
  innercount = 0;
  pointcount = psize;
  for (unsigned int v =0;v<psize; v++)
  {
    if (onboundary[v]) indexwoboundary[v]=-1;
    else
    {
      indexwoboundary[v] = innercount;
      innercount++;
    }  
  }
  
  ismanifold = true;
  if (maxt > 2)
  {
     std::cerr << "NOT Manifold: e.g. edge " << medg << " ( " << e_to_v[medg][0] <<" , " <<e_to_v[medg][1]<< " )  has " << maxt << " triangles: " << e_to_t[medg] << std::endl;
    ismanifold = false;
  }
  

  //if (ismanifold) createLoops();

  
  if (verbose > 0) cout << "   DONE" << endl;
}


/**
 * Computes tria quality (and worst, average).
 * If we have trias with qual < tmin, prints a warning to cerr.
**/
void Remesher::createQualities(double tmin) const
{
  worstqual = 1;
  averagequal = 0;
  unsigned int k;
  double t;  
  qualities.resize(tria.size());
  int badcount = 0;
  for (k= 0; k < tria.size(); k++)
  {    
    t = fabs(computeQuality(k));
    qualities[k] = t;
    if (t < tmin)
    {
      badcount ++;
//      cerr << endl;
//      cerr << "Warning (bad triangulation) "<< endl;
//      cerr << "    Triangle: "<< k << " Quality:"<< t<< endl;
    }
    if (t < worstqual) worstqual = t;
    averagequal += t;
  }
  averagequal = averagequal / tria.size();
  if (badcount > 0 && verbose > 0)
  {
      cerr << endl;
      cerr << "Warning (bad triangulation) "<< endl;
      cerr << " counted " <<badcount << " trias with Quality < " << tmin << endl;
  }
}


/** Remeshing by Botsch and Kobelt 2004
 l is the target edge length (paper: slightly less than average edge length)
 it the number of iterations (paper says 5) */
void Remesher::remeshBK(unsigned int it, double l, bool ridge)
{
  //l =  0.355856;
  
  //cout << " checkStructure: " << checkStructure() << endl;
  double qual0 = getAverageQuality();
  //  tangentialSmoothing(3,true);
  
  if (l < 0)
  {
    double f = 0.8;
    l = f * getAverageEdgeLength();
    cout << "   adjusted l: " << l << endl;
  }
  
  if (ridge) makeRidge();
  else onridge.clear();

  std::cout << "remeshing to edge length " << l << " with " << it << " iterations" << std::endl;

  // loop through iterations
  for (unsigned int ii = 0; ii<it; ii++)
  {
    // 1. Split edges at midpoint that are longer than 4/3 l
    while ( insertVerticesOnLongEdgesQueue(l*4.0/3.0) > 0) {};
    //insertVerticesOnLongEdges(l*4.0/3.0);
    init(); // v_to_t might be missing
    
    //int cs = checkStructure();
    //if (verbose > 0) cout << " checkStructure (1): " << cs << endl;
    //if (cs > 0 ) exit(cs);
    //  if (isnan(computeArea()))
    //  {cout <<"ERROR: Area is nan (1)!" << endl; exit(1);}
//exportOFF("test-i.off");
    assert(onridge.size() == points3d.size() || onridge.size() == 0);

    // 2. Collapse edges shorter than 4/5 l into midpoint
    //int ccc = 0;
    while ( contractShortEdgesQueue(l*4.0/5.0) > 0) {};
    //contractShortEdges(l*4.0/5.0);
    //if (verbose > 0) cout << " checkStructure (2): " << checkStructure() << endl;
    //      if (isnan(computeArea()))
    //  {cout <<"ERROR: Area is is nan (2)!" << endl; exit(1);}
//exportOFF("test-ii.off");
    assert(onridge.size() == points3d.size() || onridge.size() == 0);

    
    // 3. Flip edges to optimize valence (6 or 4 on boundary)
    //flipEdges();
    //if (verbose > 0)  cout << " checkStructure (3): " << checkStructure() << endl;
    //     if (isnan(computeArea()))
    // {cout <<"ERROR: Area is is nan (3)!" << endl; exit(1);}
//exportOFF("test-iii.off");

    //assert(onridge.size() == points3d.size() || onridge.size() == 0);

    // 4. Tangential Smoothing
    tangentialSmoothing(2,true);
    //     if (isnan(computeArea()))
    //   {cout <<"ERROR: Area is is nan (4)!" << endl; exit(1);}
    // if (verbose > 0)   cout << " checkStructure (4): " << checkStructure() << endl;
//exportOFF("test-iiii.off");
    //while ( contractShortEdgesQueue(l*4.0/5.0) > 0) {};
  }
 
  init();
  int cs = checkStructure();
  if (verbose > 0) cout << " checkStructure (final): " << cs << endl;
  if (cs > 0 ) exit(cs);
  
  double qual1 = getAverageQuality();
  cout << endl;
  cout << "avg qual before   : " << qual0   << "  after: " << qual1 << endl;
  cout << endl;
}


/**
  Remeshing by Botsch and Kobelt 2004
  vnum is the target vertex number
  it the number of iterations (paper says 5)
*/
void Remesher::remeshBKV(unsigned int it, int vnum, bool ridge)
{
  cout << " Remesher::remeshBKV( "<<it << " , " << vnum << " )" << endl;
  if (vnum < 0) {
    cerr << "Error: specify number of target vertices!" <<endl;
    exit(1);
  }
  double qual0 = getAverageQuality();
  
  int vcount = getPointCount();
  double avel = getAverageEdgeLength();
  double area = computeArea();
  
  double s0 = sqrt( 2 * area / (sqrt(3)* vcount));
  double st = sqrt( 2 * area / (sqrt(3)*vnum));
  
  cout << "vcount = " << vcount << "  vnum: " << vnum << endl;
  cout << "area   = " << area   << "  avel: " << avel << endl;
  cout << "s0     = " << s0     << "  st:   " << st << endl;
  if (ridge) makeRidge();
  else onridge.clear();
  //int maxnup = -1;
  int maxndown = -1;
  
  // loop through iterations
  for (unsigned int ii = 0; ii<it; ii++)
  {
    if (it > 1 && verbose > 0) cout << endl << ">>>  working on iteration " << ii+1 << endl << endl;
    // get to target vertices in it equal steps
    // compute target num of vertices and target edge length for this step
    double ii1 = ii+1.0;
    int targetn = ((it-ii1)/it)*vcount + (ii1/it) *vnum;
    double targetl = ((it-ii1)/it)*avel + (ii1/it) *st;
    
   
    cout << endl << "Targetn: " << targetn << endl << "points: " << getPointCount() << endl;
    cout << "Targetl: " << targetl << "  avg l: " << getAverageEdgeLength() << endl;
    //cout << " Maxnup: " << maxnup << endl;
   //exit(1);
    int n = 100;
    while (n>0) {n = insertVerticesOnLongEdgesQueue(targetl*4.0/3.0);};
    init(); // v_to_t might be missing  
    maxndown = getPointCount() - targetn;
    // if we want to move further up (more vertices) only go down very little:
    if (maxndown <0) maxndown = (int)(0.1*fabs(maxndown));
    // else go down slightly more than targetn:
    else maxndown = (int)(1.1 * maxndown);
    // only in last step go exactly to target n
    if (ii == it-1)
    {
      maxndown = getPointCount() - vnum;
      if (maxndown<0) { cerr<<"Error, too small before last shrinking step: " << getPointCount() << endl; exit(1);}
      // if this ever happens, we could insert a couple more and then shrink back.
    }
    cout << endl << " points: " << getPointCount() << " Maxndown:" << maxndown << endl;
    
    n=1; 
    while ( n>0 && maxndown!=0) { n= contractShortEdgesQueue(targetl*4.0/5.0,maxndown);maxndown-=n;};
    tangentialSmoothing(2,true);
    cout << " Points: " << getPointCount() << endl << "  avg edge: " << getAverageEdgeLength() << endl;
    
  }
  assert(maxndown == 0);
  cout << "Final Points: " << getPointCount() << endl << "final avg edge: " << getAverageEdgeLength() << endl;
  double qual1 = getAverageQuality();
  cout << endl;
  cout << "avg qual before   : " << qual0   << "  after: " << qual1 << endl;
  cout << endl;
}



int Remesher::checkStructure()
{
  // test if sizes are correct:
  int test = 1;
  if (tria.size() != t_to_e.size()) return test;
  test = 2;
  if (e_to_v.size() != e_to_t.size()) return test;
  test = 3;
  if (v_to_e.getNoRow() != (int)points3d.size()) return test;
  test = 4;
  if (onboundary.size() != points3d.size()) return test;
  test = 5;
  if (pointcount != (int)points3d.size()) return test;
  test = 6;
  if ((int)indexwoboundary.size() != pointcount) return test;
  test = 7;
  for (unsigned int t = 0;t<tria.size();t++)
  {
    if (tria[t].size() != 3) return test;
    if (t_to_e[t].size() != 3) return test+1;
  }
  test = 9;
  for (unsigned int e=0;e<e_to_v.size();e++)
    if (e_to_v[e].size() != 2) return test;

  // check if e_to_v is also in v_to_e
  test = 10;
  for (unsigned int e = 0;e<e_to_v.size(); e++)
  {
    int v0 = e_to_v[e][0];
    int v1 = e_to_v[e][1];
    if (v0 == -1 || v1 == -1) return test;
    if ( (int)v_to_e(v0, v1) -1 != (int)e)
    {
      cout << " v_to_e not correct (1)" << endl;
      return test;
    }
    if ( (int)v_to_e(v1, v0) -1 != (int)e)
    {
      cout << " v_to_e not correct (2)" << endl;
      return test+1;
    }
  }
  test = 12;
  // check if all edges from v_to_e exist in e_to_v
  int maxe = (int)v_to_e.getMax() -1;
  int mine = (int)v_to_e.getMin() -1;
  if (mine < -1 || maxe >= (int)e_to_v.size()) return test;

  // check if all points and edges in trias exist
  test = 13;
  for (unsigned int t = 0;t<tria.size(); t++)
  for (unsigned int j = 0;j<3;j++)
  {
    if (tria[t][j] < 0 || tria[t][j] >= pointcount) return test;
    if (t_to_e[t][j] < 0 || t_to_e[t][j] >= (int)e_to_t.size()) return test+1;
  }
  
  // check if all e_to_t exist
  test = 15;
  for (unsigned int e = 0;e<e_to_t.size();e++)
    if (e_to_t.size() ==0) return test;
  
  // check if vertices in tria are consistent with edges
  test = 16;
  for (unsigned int t = 0;t<tria.size(); t++)
  for (unsigned int j = 0;j<3;j++)
  {
    int v0 = tria[t][j];
    int v1 = tria[t][(j+1)%3];
    if (v0 == -1 || v1 == -1) return test;
    int e = t_to_e[t][j];
    // make sure v_to_e contains this edge
    if ((int)v_to_e(v0, v1) != e+1) return test;
    if ((int)v_to_e(v1, v0) != e+1) return test+1;
    // make sure edge contains these vertices
    if (e_to_v[e][0] != v0 && e_to_v[e][1] != v0 ) return test+2;
    if (e_to_v[e][0] != v1 && e_to_v[e][1] != v1 ) return test+3;
    // make sure this tria is in e_to_t list:
    unsigned int lt;
    for ( lt=0;lt < e_to_t[e].size(); lt++)
      if (e_to_t[e][lt] == (int)t) break;
    if (lt == e_to_t[e].size()) return test+4;
  }
  
  // make sure every tria is in correct e_to_t list
  test = 21;
  for (unsigned int t = 0;t<tria.size(); t++)
  for (unsigned int j = 0;j<3;j++)
  {
    int e = t_to_e[t][j];
    // make sure this tria is in e_to_t list:
    unsigned int lt;
    for ( lt=0;lt < e_to_t[e].size(); lt++)
      if (e_to_t[e][lt] == (int)t) break;
    if (lt == e_to_t[e].size()) return test;    
  }
  
  // and every tria in t_to_e has really that edge
  test = 22;
  for (unsigned int e = 0;e<e_to_t.size();e++)
  for (unsigned int lt = 0;lt<e_to_t[e].size();lt++)
  {
    unsigned int t = e_to_t[e][lt];
    unsigned int j;
    for (j = 0;j<3;j++)
      if (t_to_e[t][j] == (int)e) break;
    if (j == 3) return test;
  }

  // check onboundary
  test = 23;
  for (unsigned int e = 0;e<e_to_t.size();e++)
    if (e_to_t[e].size() ==1)
      if (!onboundary[e_to_v[e][0]] || !onboundary[e_to_v[e][1]]) return test;
  // also innercount
  test = 24;
  int icount =0 ;
  for (unsigned int p = 0;p<points3d.size();p++)
  {
    if (onboundary[p] && (indexwoboundary[p] != -1)) return test;
    if (!onboundary[p] &&  (indexwoboundary[p] < 0 || indexwoboundary[p] >= innercount) )return test+1;
    if (!onboundary[p]) icount++;
  }
  test = 26;
  if (icount != innercount) return test;

  // everything is ok
  return 0;
}


double Remesher::computeArea(unsigned int tidx) const
{
  assert(tidx < tria.size());
  unsigned int i = tidx;
  Vector cr = cross((points3d[tria[i][1]]-points3d[tria[i][0]]),(points3d[tria[i][2]]-points3d[tria[i][0]]));
  return cr.norm() / 2.0;
}


double Remesher::computeArea() const
// Computes Surface Area
{
  //if (verbose > 0) cout << "Remesher::computeArea " << flush;
  double area=0;
  Vector cr;
  for (int i=0; i<(int)tria.size();i++)
  {
    cr = cross((points3d[tria[i][1]]-points3d[tria[i][0]]),(points3d[tria[i][2]]-points3d[tria[i][0]]));
    area += cr.norm();
  }
  //if (verbose > 0) cout << " A = " << area/2.0 << endl;
  return area/2.0;
}


/** Contracts edges shorter than l.
 Default for l is AvgEdgeLength * 3.6 / 5.0.
 Runs through all existing (old) edges, skips contraction for boundary edges.
 Finally removes all collapsed triangles.
 */
int Remesher::contractShortEdges(double l, int maxn)
{
  if (l == -1) l = 0.9 * getAverageEdgeLength() * 4.0/5.0;
  if (verbose > 0) cout << " Remesher::contractShortEdges( "<<l<<" ) " << endl;

  createElementN(true); // need to check if normals flip

  int en = 0;
  //int es = e_to_v.size();
  int v1,v2;
  int bskip = 0;
  //vector < int > rmt;
  for (int i = 0;i<(int)e_to_v.size();i++) // only run trough all old edges (ignore newly created edges in this run)
  {
    if (e_to_t[i].size() == 0) continue; // edge was deleted
    v1 = e_to_v[i][0];
    v2 = e_to_v[i][1];
    //cout << " checking edge: " << i << " v1 : " << v1 << " v2: " << v2 << endl;
    if ((points3d[v1] - points3d[v2]).norm() < l)
    {
      //cout << " found short edge: " << i << endl;
      //cout << " e_to_v[i] " << e_to_v[i] << endl;
      if (e_to_t[i].size() < 2 || onboundary[v1] || onboundary[v2])
      {
        //cout << " skipping boundary ( e_to_t[i].size: " << e_to_t[i].size() << " ) " << endl;  
        bskip++;
        continue;
      }
//       if ( onridge.size() == points3d.size() && ( onridge[v1] || onridge[v2] )) 
//       {
//         cout << " skipping ridge ( e_to_t[i].size: " << e_to_t[i].size() << " ) " << endl;  
//         continue;
//       }
         

      // store tria idxs for deletion later
      //for (unsigned int j = 0; j<e_to_t[i].size();j++)
      //  rmt.push_back(e_to_t[i][j]);
        
      if (contractEdge(i,false)) en++;
//       else
//       {
//          Remesher T = computePatchT(e_to_t[i][0],5);
//          T.exportOFF("test-patch.off");
//          exit(1);
//       }
      
//if (rmt.size() % 5000 == 0) 
//   if (i == 364)
//   {
//      Remesher T(*this);
//      T.exportOFF("test-ii-sub.off");
//      //getchar();
//   }
// 
//       int cs = checkStructure();
//       if (cs > 0 || ! ismanifold )
//        {
//          cout << " check: " << cs << endl << endl << endl;    
//          exit(cs);         
//       }
        if (en == maxn) break;
    }
  
  }
  // remove all collapsed trias:
  rmTrias(true);

  if (verbose > 0) cout << "  contractShortEdges: " << en << " edges contracted ( "<< bskip <<" boundary edges skipped)" << endl;
  
  return en;

}


int Remesher::contractShortEdgesQueue(double l, int maxn)
{
  if (l == -1) l = 0.9 * getAverageEdgeLength() * 4.0/5.0;
  if (verbose > 0) cout << " Remesher::contractShortEdgesQueue( " << l << " ) " << endl;
  createElementN(true); // need to check if normals flip
  int en = 0;
  //int es = e_to_v.size();
  int v1,v2;
  int bskip = 0;
  //vector < int > rmt;
  double len = 0.0;
  
  std::priority_queue<std::pair<double, int>, std::vector<std::pair<double,int> >, std::greater<std::pair<double,int> >  > q;
  
  for (int i = 0;i<(int)e_to_v.size();i++) // only run trough all old edges (ignore newly created edges in this run)
  {
    if (e_to_t[i].size() == 0) continue; // edge was deleted
    v1 = e_to_v[i][0];
    v2 = e_to_v[i][1];
    if (e_to_t[i].size() < 2 || onboundary[v1] || onboundary[v2])
    {
      //cout << " skipping boundary ( e_to_t[i].size: " << e_to_t[i].size() << " ) " << endl;  
      bskip++;
      continue;
    }
    len = (points3d[v1] - points3d[v2]).norm();
    q.push(std::pair<double,int>(len,i));
  }
  
  while (!q.empty() && q.top().first < l && maxn != en)
  {
    int i = q.top().second;
    len = q.top().first;
    q.pop();
    if (e_to_t[i].size() == 0) continue; // edge was deleted
    v1 = e_to_v[i][0];
    v2 = e_to_v[i][1];
    double ll = (points3d[v1] - points3d[v2]).norm();
    if (ll > len) q.push(std::pair<double,int>(ll,i));
    else if (contractEdge(i,false)) en++;
  }
  // remove all collapsed trias:
  rmTrias(true);
  if (verbose > 0) cout << "  contractShortEdgesQueue: " << en << " edges contracted ( "<< bskip <<" boundary edges skipped)" << endl;
  
  return en;
}


double Remesher::getAverageEdgeLength() const
{
  double l = 0;
  for (unsigned int ei = 0; ei<e_to_v.size(); ei++) {
    l += (points3d[e_to_v[ei][0]] - points3d[e_to_v[ei][1]]).norm();
  }

  l /= e_to_v.size();

  return l;
}


/** Inserts a vertex on the midpoint of all edges longer than l.
 Returns number of inserted vertices. */
int Remesher::insertVerticesOnLongEdges(double l, int maxn)
{
  if (l == -1) l = 0.9 * getAverageEdgeLength() * 4.0/3.0;
  if (verbose > 0) cout << "Remesher::insertVerticesOnLongEdges( "<<l<< " )" << endl;

  int es = e_to_v.size();
  int vnum = 0;
  int v1,v2;
  for (int i = 0;i<es;i++) // only run trough all old edges (ignore newly created edges in this run)
  {
    v1 = e_to_v[i][0];
    v2 = e_to_v[i][1];
    if ((points3d[v1] - points3d[v2]).norm() > l)
    {
      vnum++;
      insertVertexOnEdge(Vector(v1,v2,0.5));
      // cout << " check: " << checkStructure() << endl;
    }
    if (vnum == maxn) break;  
  }  
  if (verbose > 0) cout << "  " << vnum << " vertices inserted! " << endl;
  return vnum;
}


int Remesher::insertVerticesOnLongEdgesQueue(double l, int maxn)
{
  if (l == -1) l = 0.9 * getAverageEdgeLength() * 4.0/3.0;
  if (verbose > 0) cout << "Remesher::insertVerticesOnLongEdgesQueue( " << l << " )" << endl;
  int es = e_to_v.size();
  int vnum = 0;
  int v1,v2;
  double len = 0.0;
  std::priority_queue<std::pair<double, int> > q;
 
  for (int i = 0;i<es;i++) // only run trough all old edges (ignore newly created edges in this run)
  {
    v1 = e_to_v[i][0];
    v2 = e_to_v[i][1];
    len = (points3d[v1] - points3d[v2]).norm();
    q.push(std::pair<double,int>(len,i));
  }
  
  while (!q.empty() && q.top().first > l && maxn != vnum)
  {
    int i = q.top().second;
    len = q.top().first;
    q.pop();
    if (e_to_t[i].size() == 0) continue; // edge was deleted
    v1 = e_to_v[i][0];
    v2 = e_to_v[i][1];
    double ll = (points3d[v1] - points3d[v2]).norm();
    if (ll < len) q.push(std::pair<double,int>(ll,i));
    else 
    {
      insertVertexOnEdge(Vector(v1,v2,0.5));
      vnum++;
    }
  }
  if (verbose > 0) cout << "  " << vnum << " vertices inserted! " << endl;
  return vnum;
}


/** Places vertices on ridge at edges with large angles.
 angle can be 0..PI (although pi/6 to pi/3 makes more sense) */
int Remesher::makeRidge(double angle)
{
  if (angle < 0) angle = M_PI/3.0; // 60 degree
  if (!ismanifold) return 0;
  onridge.clear();
  onridge.resize(points3d.size(),false);
  double c = cos(angle);
  cout << " Remesher::makeRidge( " << angle << " )" <<endl;
  createElementN(true);
  int ret = 0;
  for (unsigned int e = 0;e<e_to_t.size();e++)
  {
     if (e_to_t[e].size() < 2) continue;
     assert (e_to_t[e].size() == 2);
     if (elementn[e_to_t[e][0]] * elementn[e_to_t[e][1]] < c)
     {
        onridge [e_to_v[e][0]] = true;
        onridge [e_to_v[e][1]] = true;
        ret++;
     }
  }
  
  return ret;
}


/** Smooth mesh vertices.
  gravity: vertices with larger area attract stronger.
  tangent: switch off attraction towards tangent plane.
  If tangent = true averages will be moved back into tangent plane (using damping).*/
void Remesher::tangentialSmoothing(int it, bool gravity, bool tangent)
{

  if (verbose > 0) cout << "Remesher::tangentialSmoothing( " <<  it << " , bool gravity)" << endl;
  vector < double > A;
  
  for (int run = 0 ; run < it; run++)
  {  
  if (gravity)
  {
    A.clear();
    A.resize(points3d.size(),0.0);
    double ta;
    for (int t = 0;t<(int)tria.size();t++)
    {
      ta = computeArea(t)/3.0;
      for (int i = 0;i<3;i++)
        A[tria[t][i]]+= ta;
    }
  }
  else
  {
    A.clear();
    A.resize(points3d.size(),1.0);
  }

  createVertexNSimple(true);
  
  Vector gi;
  double d;
  for (unsigned int i = 0;i<points3d.size();i++)
  {
    if (onboundary[i]) continue;
    if (onridge.size() == points3d.size() && onridge[i]) continue;
    
    const vector<int> & N = get1Ring(i,false);
    if (N.empty()) continue; // if nonmanifold
    gi = Vector(0,0,0);
    d = 0;
    for (unsigned int j = 0;j<N.size();j++)
    {
      d += A[N[j]];
      gi += A[N[j]] * points3d[N[j]];
    }
    gi = (1.0/d) * gi;

    if (tangent)
    {
      // update into tangent plane
      double lambda = 0.99; // damping factor
      Matrix I(1,0,0,0,1,0,0,0,1);
      Matrix M;
      for (int ii = 0; ii<3;ii++)
      for (int jj = 0; jj<3;jj++)
        M(ii, jj) = vertexn[i][ii] * vertexn[i][jj];
      points3d[i] = points3d[i] + lambda * ((I - M ) * (gi -points3d[i]));
    }
    else
    {
      points3d[i] = gi;
    }
    
  }
  }
}


/** Returns indices of 1 ring vertices.
 If closed, last index will equal first,
 else it will start at a boundary and to the other.
 If not ordered and vertex is at a non-manifold edge, return empty set.
 If ordered, we want to have manifold!*/
std::vector < int > Remesher::get1Ring(unsigned int vn, bool ordered) const
{
  
  assert(vn < points3d.size());
  
  //cout << "Remesher::get1Ring( "<<vn<<" )" << endl;
  
  // get reference to row (neigbors of vertex vn)  
  const std::list < SparseMatrix::entry > & row = v_to_e.getRow(vn);
  vector < int > ret;             // return vector
  
  if (!ordered)
  {
    list<SparseMatrix::entry>::const_iterator iter;
     for (iter=row.begin(); iter != row.end(); iter++)
     {
       ret.push_back((*iter).pos);
      int e = (int)iter->value -1;
      if (e_to_t[e].size() > 2) {ret.clear(); return ret;}
     }
     //cout << " neigbors: " << ret << endl;  
  }
 else 
 {
  if (!ismanifold) 
  {
    cout << "Remesher::get1Ring: Ordered only possible for Manifolds!" << endl;
    return vector < int >() ;
  }
  // set up variables and
  // push first vertex of 1ring (arbitrary):
  int ec = (int)(*row.begin()).value -1; // current edge
  int vc = (*row.begin()).pos;    // current vertex (neighbor of vn)
  int tc = e_to_t[ec][0];         // current tria
  int tother = -1;                // other tria at current edge
  bool hitborder = false;         // signal if we hit a boundary
  if (e_to_t[ec].size()>1) tother =  e_to_t[ec][1];
  else hitborder = true; // we have a border at the start!
  
  ret.push_back(vc);
  
  // move around:
  int j,temp;
  do
  {
    //cout << "vc: " << vc <<  "  ec: " << ec<< "  tc: " << tc << endl;
    //cout << "tria: " << tria[tc] << endl;
    //cout << "edges: " << t_to_e[tc] << endl;
    // find next vertex in current triangle:
    for (j = 0; j< 3;j++)
      if (tria[tc][j] != (int)vn && tria[tc][j]!=vc) break;
    assert (j < 3);
    vc = tria[tc][j];
    ret.push_back(vc);
    
    // find current edge and then next triangle:
    ec = (int)(v_to_e.getVal(vc,vn))-1;
    
    if (e_to_t[ec].size()<2)
    {
      // we hit a boundary:
      //cout << " hit a boundary! " << endl;
      // if this happens for the second time, we are done:
      if (hitborder) break;
      
      hitborder = true;
      
      // reverse order of ret
      int rs = ret.size();
      for (j=0;j<rs/2;j++)
      {
        temp = ret[j];
        ret[j] = ret[rs-1-j];
        ret[rs-1-j] = temp;
      }
      // setup values to move further in the other direction
      tc = tother;
      vc = ret[rs-1];
    
    }
    // here we have another neighbor triangle:
    else if (e_to_t[ec][0] == tc) tc = e_to_t[ec][1];
    else
    {
      assert(e_to_t[ec][1] == tc);
      tc = e_to_t[ec][0];
    }
  
  } while(ret[0] != ret.back());
  
  
  //cout << "ret: " << ret << endl;
  //cout << "rowsize: " << row.size() << endl;
 }

  
  int rs = ret.size();
  if (ret[0] == ret[rs-1]) rs--;
  assert(rs == (int)row.size());
  
  return ret;
}


/** Simply average tria normals at vertex (no angle normalization)
 to create vertex normals.*/
void Remesher::createVertexNSimple(bool force) const
{


  //cout << "Tria3d.createVertexN()" << endl;
  if (vertexn.size() == points3d.size()  && ! force) return;
  // erstellt zuerst Element Normalen fuer alle Elemente
  createElementN(force);
  int numpoint = points3d.size();
  vertexn.clear();
  Vector empty(0,0,0);
  vertexn.resize(numpoint,empty);

  int idx, elem, knoten;
  // laufe ueber alle Dreiecke und ueber jede Ecke
  for (elem =0; elem < (int)tria.size(); elem++)
  for (knoten=0; knoten < 3; knoten++)
  {
    //cout << "Element: "<< elem << " Knoten: " << knoten << endl;
    idx = tria[elem][knoten];
    
//     if (idx == 178)
//     {
//       cout << "Element: "<< elem << " Knoten: " << knoten << endl;
//       cout << "   tnormal: " << elementn[elem] << "  angle: " << getAngle(elem,knoten) << endl;
//     
//     
//     }
    vertexn[idx] += elementn[elem];
    
  }
  
  // check length of vertexn and normalize
  double length;
  for (idx=0;idx<(int)vertexn.size();idx++)
  {
    
      length = vertexn[idx].norm();
      if (length < 0.0000001)
      {
        cout << "Remesher.createVertexN "<< idx<<": normal length is zero!" << endl;
        // get unordered list of triangles at this vertex
        //vector < int > tn =  getv_to_t(idx,false); 
        const vector < int >& tn =  v_to_t[idx]; 
        cout << "  tn : " << tn << endl;
        for (unsigned int ii=0;ii<tn.size(); ii++)
        {
          cout << " tria " << tn[ii] << " : " << tria[tn[ii]] << " tnormal: " << elementn[tn[ii]] << endl;
        }
        //exit(1);
      }
      vertexn[idx] *= (1.0/length);
      //cout << "Vertex normal: " << vertexn[idx] << endl;
      //cout << "Tangent U = " << getTangentU(vertexn[idx]) << endl ;
      //cout << "Tangent V = " << getTangentV(vertexn[idx]) << endl;
      //cout << "U x V = " << cross(getTangentU(vertexn[idx]), getTangentV(vertexn[idx])) << endl;;
      
  }

  //cout << " points3d: " << endl << points3d << endl;
  //cout << " Vertex Normals : " << endl << vertexn << endl;
  
}


/** The new vertex will be inserted at w*v1 + (1-w)*v2
 where pos is a vector int vertex1, int vertex2, double w.
 Datastructures are updated, but not v_to_t
 */
int Remesher::insertVertexOnEdge(const Vector& pos)
{
  int v1 = (int) pos[0];
  int v2 = (int) pos[1];
  double w = pos[2];
  
  // swap so that v1 < v2
  if (v1 > v2)
  {
    int ti = v1;
    v1 = v2;
    v2 = ti;
    w = 1.0 - w;
  }
  
  assert(v1 < (int)points3d.size());
  assert(v2 < (int)points3d.size());
  int e = (int)v_to_e(v1, v2) - 1;
  assert (e > -1);
  assert ((e_to_v[e][0] == v1 && e_to_v[e][1] == v2) || (e_to_v[e][0] == v2 && e_to_v[e][1] == v1));
  
//cout << " working on edge: " << e << "  v1: " << v1 << " v2: " << v2 << endl;  
  
  // insert new point
  //Vector nv = points3d[v1]*w + points3d[v2]*(1.0-w);
  
////  cout << " new vertex: " << nv << endl;
//  if (isnan(nv[0])) exit(1);
  
//  int ni = points3d.size();
//  points3d.push_back(nv);
  int ni = createVertexOnEdge(pos);
  pointcount++;
  if (e_to_t[e].size() == 1)
  {
    onboundary.push_back(true);
    indexwoboundary.push_back(-1);
  }
  else
  {
    onboundary.push_back(false);
    indexwoboundary.push_back(innercount);
    innercount++;
  }
  if (onridge.size() == points3d.size() -1)
  {
     if (onridge[v1] && onridge[v2])
        onridge.push_back(true);
     else onridge.push_back(false);
     // not really great as it is not clear that this edge is on a ridge,
     // if the two endpoints are on ridge, edge could still be inside
  }  
  
  // insert new edge (ni to v2) (e_to_v and v_to_e)
  vector < int > newedge(2); newedge[0] = ni; newedge[1] = v2;
  int en = e_to_v.size();
  e_to_v.push_back(newedge);
  v_to_e.redim(ni+1,ni+1);
  v_to_e(v2, ni) = en+1;
  v_to_e(ni, v2) = en+1;
  //   and adjust old edge (v1 to ni)
  if (e_to_v[e][0] == v2) e_to_v[e][0] = ni;
  else
  {
    assert(e_to_v[e][1] == v2);
    e_to_v[e][1] = ni;
  }
  v_to_e(v1, ni) = e+1;
  v_to_e(ni, v1) = e+1;
  v_to_e.erase(v1,v2); // delete old edge
  v_to_e.erase(v2,v1);
  // create empty e_to_t[en]
  e_to_t.resize(e_to_t.size()+1);

  // create new triangles and fix structure
  vector < int > etote;
  for (int ti = 0 ; ti<(int)e_to_t[e].size(); ti++)
  {
    int t = e_to_t[e][ti];
    
    //cout << "  " << ti<<":   working on tria " << t << endl;
  
    // compute local index of other point in this tria:
    int vother = 0;
    for (vother =0; vother< 3; vother++) if (tria[t][vother] != v1 && tria[t][vother] != v2)break;
    assert(vother < 3);

    // keep track of the direction of the edge
    bool reverse = false;
    //cout << " v1: " << v1 << "  v2: " << v2 << endl;
    //cout << " tria[t]: " << tria[t] << endl;
    if (tria[t][(vother+1)%3] == v2)
    {
      assert (tria[t][(vother+2)%3] == v1);
      reverse = true;
    }
    else
    {
      assert(tria[t][(vother+1)%3] == v1);
      assert(tria[t][(vother+2)%3] == v2);
    }
    
    // create new edge within tria (update e_to_v and v_to_e)
    newedge[0] = tria[t][vother]; 
    newedge[1] = ni;
    int eint = e_to_v.size();
    e_to_v.push_back(newedge);
    v_to_e(newedge[0], ni) = eint+1;
    v_to_e(ni, newedge[0]) = eint+1;

    // add new tria  
    vector < int > t2(3);
    t2[0] = ni;
    t2[1] = tria[t][vother];
    t2[2] = tria[t][(vother+1)%3];
    int t2i = tria.size();
    tria.push_back(t2);
    
    // also add new t_to_e (misuse variable t2)
    t2[0] = eint;
    t2[1] = t_to_e[t][vother];
    if (reverse) t2[2] = en;
    else t2[2] = e;
    t_to_e.push_back(t2);
    
    // adjust e_to_t at edge now pointing to t2i (instead of t)
    int tte = t_to_e[t][vother];
    int cc = 0;
    for (unsigned int tt = 0; tt< e_to_t[tte].size(); tt++)
      if (e_to_t[tte][tt] == t)
      { e_to_t[tte][tt] = t2i;cc++;}
    assert(cc ==1 );


    // update also e_to_t[eint] (misuse variable newedge)
    newedge[0] = t; newedge[1] = t2i;
    e_to_t.push_back(newedge);

    // adjust old tria t to be the other half
    tria[t][(vother+1)%3] = ni;
    
    // adjust t_to_e[t]
    t_to_e[t][vother] = eint;
    if (reverse) t_to_e[t][(vother+1)%3] = e;
    else t_to_e[t][(vother+1)%3] = en;
    // the third edge stays the same
    
     
    // adjust e_to_t for the split edge e (and en)
    // e_to_t[en] gets a new tria
    // e_to_t[e] has to be changed later (we are working with it)
    // store tria and add later    
    if (reverse)
    {
      e_to_t[en].push_back(t2i);
      etote.push_back(t);
    }
    else
    {
      e_to_t[en].push_back(t);
      etote.push_back(t2i);
    }
    
        
  }
  
  // now get rid of old e_to_t[e] and replace with etote
  e_to_t[e] = etote;
  return pointcount-1;
}


/** Contracts edge eidx to its midpoint.
 e_to_t, t_to_e, trias and e_to_v are kept intact.
 cleanup true means that the empty trias are deleted from the tria vector (using rmTria).*/
bool Remesher::contractEdge(int eidx, bool cleanup)
{
  //cout << endl << "Remesher::contractEdge( " << eidx << " )" << endl;

  if (e_to_t[eidx].size() ==0)
  {
    if (verbose > 2)
      cout << " edge " << eidx << " already deleted! " << endl;
    return false; // if this edge has allready been voided
    
  }
  
  assert(ismanifold);  
  //assert(bedgecount==0);

  // find two vertices at this edge
  int v0 = e_to_v[eidx][0];
  int v1 = e_to_v[eidx][1];
  //cout << " v0 : " << v0 << "  v1: " << v1 << endl;
  // find two triangles at this edge
  assert (e_to_t[eidx].size() ==2); // boundary edges not allowed
  int t0 = e_to_t[eidx][0];
  int t1 = e_to_t[eidx][1];
  // find two opposite vertices in triangles
  int vt0 = tria[t0][0];
  if (vt0 == v0 || vt0 == v1) vt0 = tria[t0][1];
  if (vt0 == v0 || vt0 == v1) vt0 = tria[t0][2];
  assert (vt0 != v0 && vt0 != v1);
  int vt1 = tria[t1][0];
  if (vt1 == v0 || vt1 == v1) vt1 = tria[t1][1];
  if (vt1 == v0 || vt1 == v1) vt1 = tria[t1][2];
  assert (vt1 != v0 && vt1 != v1);

  
  // before we start, check if we would get special collapsing super triangles
  int super = replaceVertexInEdges(v1,v0,true,vt0,vt1); //simulate only
  if (super > 0) //here special treatment would be needed, e.g. remove everyting inside super tria
  //   once new code is added here, we do not need to treat this further down
  {
     if (verbose > 2)
         cout << " there would be " << super << " collapsing super trias. Skipping edge for now ..." << endl;
     return false;
  }
  
  // possible new edge midpoint
  Vector midpoint = 0.5*(points3d[v0]+points3d[v1]);  
  bool ridge = false;
  if (onridge.size() == points3d.size())
  {
    // if one is on ridge, do not use midpoint, but ridge point
    if (onridge[v0] && ! onridge[v1])
    { midpoint = points3d[v0]; ridge = true;}
    if (onridge[v1] && ! onridge[v0])
    { midpoint = points3d[v1]; ridge = true;}
    if (onridge[v0] && onridge[v1]) 
    {
      ridge = true;
      if (verbose > 2)
         cout << " both endpoints on ridge, skipping..." << endl;
      return false;
    }
  }
  
  
  // also check if any of the new trias will flip their tnormal:
  double eps = cos(M_PI/3.0);
  for (unsigned int t = 0;t<v_to_t[v0].size();t++)
  {
    int tt = v_to_t[v0][t];
    if (tt == t0 || tt == t1) continue;
    int lv = v_to_lv[v0][t];
    assert(tria[tt][lv] == v0);
    int pb = tria[tt][(lv+1)%3];
    int pc = tria[tt][(lv+2)%3];
    //cout << " midpoint " << midpoint << " pb " << points3d[pb] << " pc " << points3d[pc] << endl;
    Vector a = points3d[pb] - midpoint;
    Vector b = points3d[pc] - midpoint;
    //cout << " a: "<< a << "  b: " << b << endl;
    Vector tn = cross(a, b);
    //cout << " cross : " << tn  <<" elementn [" << tt << "] : " << elementn[tt] <<endl;
    double proj = tn * elementn[tt] / tn.norm();
    //cout <<  "  proj: " << proj <<endl;
    if (proj < eps)
    {
       if (verbose > 2) 
         cout<< " tria: " << tt << " at vertex " << v0 << " would flip! Skipping ..." << endl;
       return false;
    
    }
  }
  for (unsigned int t = 0;t<v_to_t[v1].size();t++)
  {
    int tt = v_to_t[v1][t];
    if (tt == t0 || tt == t1) continue;
    int lv = v_to_lv[v1][t];
    assert(tria[tt][lv] == v1);
    int pb = tria[tt][(lv+1)%3];
    int pc = tria[tt][(lv+2)%3];
    Vector tn = cross((points3d[pb] - midpoint), (points3d[pc] - midpoint));
    double proj = tn * elementn[tt] / tn.norm();
    if (proj < eps)
    {
       if (verbose > 2)
         cout<< " tria: " << tt << " at vertex " << v0 << " would flip! Skipping ..." << endl;
       return false;
    
    }
  }

  
  // move point v0 to midpoint:
  points3d[v0] = midpoint;
  if (ridge) onridge[v0] = true;
  // points3d done (v1 will be keept but not used)
  
  // deal with edge structures in each tria
  contractEdgeInTria(eidx,t0);
  contractEdgeInTria(eidx,t1);
  // e_to_t, e_to_v, v_to_e are now correct (one edge was deleted for each tria) !!!
  // and t_to_e is correct in neighbors (no references to t0 and t1) !!!
  
  // replace v1 with v0 in all trias at v1
  for (unsigned int lt = 0; lt < v_to_t[v1].size(); lt++)
  {
    int tnum = v_to_t[v1][lt];
    int lvid = v_to_lv[v1][lt];
    
    //cout << " tnum: " << tnum << " lvid: " << lvid << endl;
    tria[tnum][lvid] = v0;
  }
  // tria is up to date now (except for deleting t0 t1) !!!
  
  // delete t0 and t1 from v_to_t at the four relevant vertices
  rmTriaInVtoT(v0,t0);
  rmTriaInVtoT(v0,t1);
  rmTriaInVtoT(v1,t0);
  rmTriaInVtoT(v1,t1);
  rmTriaInVtoT(vt0,t0);
  rmTriaInVtoT(vt1,t1);
  // merge v_to_t at v0 (with v1)
  v_to_t[v0].insert(v_to_t[v0].end(), v_to_t[v1].begin(), v_to_t[v1].end());
  v_to_lv[v0].insert(v_to_lv[v0].end(),v_to_lv[v1].begin(),v_to_lv[v1].end());
  // clear v_to_t[v1]
  v_to_t[v1].clear();
  v_to_lv[v1].clear();
  // v_to_t and v_to_lv is up to date now !!!
  //  (no more t0, t1) also no info at v1
  
  // rm edge vrom vtoe list
  rmEdgeInVtoE(eidx);
  
  // now vtoe has only valid edges at v1
  // replace v1 with v0 in all edges at v1 (fix e_to_v)
  int special = replaceVertexInEdges(v1,v0);
  
  // rm edge
  e_to_t[eidx].clear();
  e_to_v[eidx].clear();
  
  // rm t0 and t1
  tria[t0].clear();
  tria[t1].clear();
  t_to_e[t0].clear();
  t_to_e[t1].clear();
  
  
//     if (eidx == 364)
//   {
//      Remesher T(*this);
//      T.exportOFF("test-ii-sub-before.off");
//      //getchar();
//   }
// 
  // if something collapsed, fix it:
  int nextv = vt0;
  int sp = 0;
  while ( nextv != -1)
  { 
    nextv = removeCollapsedTrias(nextv,v0); 
    if (nextv != -1) sp ++;
    //cout << nextv << endl;
  }
  nextv = vt1;
  while (nextv != -1)
  {
     nextv = removeCollapsedTrias(nextv,v0);
     if (nextv != -1) sp++;
     // cout << nextv << endl;
  }
  if (special > 0 && verbose > 1)
  {
    cout << " special cases in replaceVertexInEdges: " << special << endl;
    cout << " special cases when removing Collapsed Trias: " << sp << endl;
  }
  assert(sp == special);
  return true;
}


/**  The new vertex will be inserted at w*v1 + (1-w)*v2 
     where pos is a 3d vector: int vertex1, int vertex2, double w.
     This function is virtual to be overloaded by Remesherfunc. */
int Remesher::createVertexOnEdge(const Vector& pos)
{
  int v1 = (int) pos[0];
  int v2 = (int) pos[1];
  double w = pos[2];
  assert(v1 < (int)points3d.size());
  assert(v2 < (int)points3d.size());
  // swap so that v1 < v2
  if (v1 > v2)
  {
    int ti = v1;
    v1 = v2;
    v2 = ti;
    w = 1.0 - w;
  }
  Vector nv = points3d[v1]*w + points3d[v2]*(1.0-w);
  if (std::isnan(nv[0]))
  {
    cerr << " createVertexOnEdge ERROR point nan? " << nv << endl;
    exit(1);
  }
  points3d.push_back(nv);
  return points3d.size()-1;  
}


void Remesher::createElementN(bool force) const
// computes and caches the tria normals
{
  //cout << "Tria3d.createElementN()" << endl;
  if (elementn.size() == tria.size() && ! force) return;
  unsigned int i;
  elementn.resize(tria.size());
  Vector p1,p2,p3,p2mp1,p3mp1, n;
  for (i=0; i<tria.size(); i++)
  {
    p1 = points3d[tria[i][0]];
    p2 = points3d[tria[i][1]];
    p3 = points3d[tria[i][2]];
    p2mp1 = p2 - p1;
    p3mp1 = p3 - p1;
    n = cross(p2mp1,p3mp1);
    if (n.norm() < 0.00000000001) elementn[i] = Vector(0,0,0); // zero tria
    else elementn[i] = (1.0/n.norm()) * n;
  }
  //cout << " Element Normals : " << endl << elementn << endl;
  
}


/** Removes labeled trias:
 ones with one index == -1 or if tria[i].size() == 0.*/
bool Remesher::rmTrias(bool fix)
{
  bool found = false;
  for (unsigned int t = 0;t<tria.size();t++)
  {
    if (tria[t].size() == 0 || tria[t][0] == -1 || tria[t][1] == -1 || tria[t][2] == -1)
    {
        tria[t] = tria.back();
        tria.pop_back();
        t--;
        found = true;
    }
  }
  if (fix || found) rmFreeVertices(fix);
  return found;
}


/**
  Removes hanging edges (with no tria on either side).
  Clears e_to_v and removes edge from v_to_e.
  Note: does not check if it was removed from any t_to_e 
  (datastructe would be invalid if edge does not have trias, but
  still is listed in any t_to_e).
*/
bool Remesher::removeHangingEdge(unsigned int eidx)
{
  assert(eidx < e_to_v.size());
  if (e_to_t[eidx].size() >0) return false; // cannot be deleted
  assert (e_to_v[eidx].size() == 2);
  
  int a = e_to_v[eidx][0];
  int b = e_to_v[eidx][1];
  int e = (int)v_to_e.eraseGetVal(a,b) - 1;
  if (e != (int) eidx)
  {
    cerr << " hanging edge: " << eidx << " e_to_v: " << e_to_v[eidx] << " but e " << e << endl;
  }
  assert (e == (int) eidx);
  e = (int)v_to_e.eraseGetVal(b,a) - 1;
  if (e != (int) eidx)
  {
    cerr << " hanging edge: " << eidx << " e_to_v: " << e_to_v[eidx] << " but e " << e << endl;
  }
  assert (e == (int) eidx);
  
  e_to_v[eidx].clear();
  return true;
}


/**
  Removes triangle tidx. 
  If inside, creates a hole and boundary loop,
  if at the boundary, the boundary edges are deleted
*/
bool Remesher::removeTria(unsigned int tidx)
{
  cout << "Remesher::removeTria( "<<tidx<<" )" << endl;
  assert(tidx < tria.size());
  if (tria[tidx].size() < 3) return false;

  // for each edge in tria
  int found = 0;
  for (unsigned int e = 0;e<3;e++) {
    int eid = t_to_e[tidx][e];
    cout << " e_to_t[ " << eid << " ] before: " << e_to_t[eid] << endl;
    // 1. remove tidx from e_to_t
    assert(e_to_t[eid].size() > 0);
    found = 0;
    for (unsigned int t = 0;t<e_to_t[eid].size(); t++)
      if (e_to_t[eid][t] == (int) tidx)
      {
        e_to_t[eid][t] = e_to_t[eid].back();
        e_to_t[eid].pop_back();
        found++;
        t--;
      }
    assert(found == 1);
    cout << " e_to_t[ " << eid << " ] after: " << e_to_t[eid] << endl;

    // 2. if edge is hanging (no tria), remove it
    if (e_to_t[eid].size() == 0) {
      removeHangingEdge(eid);
      bedgecount--;
    }
    // 3. if edge has become boundary:
    if (e_to_t[eid].size() == 1) bedgecount++;

  }

  // for each vertex in tria:
  for (unsigned int v = 0; v<3; v++)
  {
    int vidx = tria[tidx][v];
    
    //  remove tidx from v_to_t (and v_to_lv)
    assert(v_to_t.size() >1);
    found = 0;
    for (unsigned int t = 0; t<v_to_t[vidx].size(); t++) {
      if (v_to_t[vidx][t] == (int)tidx) {
        v_to_t[vidx][t] = v_to_t[vidx].back();
        v_to_t[vidx].pop_back();
        v_to_lv[vidx][t] = v_to_lv[vidx].back();
        v_to_lv[vidx].pop_back();
        found++;
        t--;
      }
    }
    assert(found == 1);
    
    // hanging vertices are 'deleted' if v_to_t empty,
    //   and if not referenced in tria or e_to_t
    // this should be the case if v_to_t[vidx].size==0
  }


  tria[tidx].clear();
  t_to_e[tidx].clear();
  return true;
}


/**
  Removes reference to tria tidx from 
  v_to_t list (and parallel v_to_lv list).
  Returns true if deleted.
  Returns false if not found.
*/
bool Remesher::rmTriaInVtoT(int vidx, int tidx)
{

   bool found = false;
   for (unsigned int i = 0;i<v_to_t[vidx].size();i++)
     if (v_to_t[vidx][i] == tidx)
     {
       found = true;
       v_to_t[vidx][i] = v_to_t[vidx].back();
       v_to_t[vidx].pop_back();
       v_to_lv[vidx][i] = v_to_lv[vidx].back();
       v_to_lv[vidx].pop_back();
      }
      
   return found;
}


/** 
  Computes the quality of tria #k
**/
double Remesher::computeQuality(int k) const
{
//  cout << "Remesher::computeQuality("<< k << ")" << endl;

  
  Vector p1 = points3d[tria[k][0]];
  Vector p2 = points3d[tria[k][1]];
  Vector p3 = points3d[tria[k][2]];
  
  Vector d1 = p2-p1;
  Vector d2 = p3-p2;
  Vector d3 = p3-p1;
  
  // circumference square:
  double lsumsq = d1.normSquared() + d2.normSquared() + d3.normSquared();
  
  // area
  double a = 0.5* (cross(d1,d3)).norm();
  
    
  double q = 4.0 * sqrt(3.0) * a / lsumsq;
        

  return q;
}


/** Replace vold with vnew in all edges at vold.
 Also fix v_to_e simultaneously.
 Treat special case of collapsing other edges!
 ret value is the number of special cases (collapsing super trias).
 If simulate = true only simulate, do not actually replace anything.
 This is used for counting the special cases that would arise.
 vtip0 and vtip1 must be given when simulating so that we can skip these edges
 as they will be removed later anyway and will otherwise pop up 
 as special cases (but are in the same triangle as vold--vnew). */
int Remesher::replaceVertexInEdges(int vold, int vnew, bool simulate, int vtip0, int vtip1)
{
  // get all edges touching vold from v_to_e:
  std::list < SparseMatrix::entry >  row = v_to_e.getRow(vold);
  list<SparseMatrix::entry>::const_iterator iter;
  int ret = 0;
    for (iter=row.begin(); iter != row.end(); iter++)
    {
      int edge = (int)iter->value -1;
      assert(edge >= 0);
      assert(edge < (int) e_to_v.size());
      if (e_to_v[edge].size() != 2)
      {
        cout << "replaceVertexInEdges ERROR: more than two vertices? e_to_v[" << edge << "] : " << e_to_v[edge] << endl;
        assert(e_to_v[edge].size() == 2);
      }
      int vother = -1;

      //cout << " working edge at v1 ( " << v1 << " ) : " << edge << flush;
      //cout << "  iter->pos: " << iter->pos << flush;
      //if (edge > -1) cout << " etov: " << e_to_v[edge] << endl;
      
      // replace vold with vnew in this edge:
      assert (e_to_v[edge][0] == vold || e_to_v[edge][1] == vold);
      if (e_to_v[edge][0] == vold)
      {
        if (e_to_v[edge][1] == vnew) continue; // skip this is the edge, it  will be deleted later (see also below)
        if (vtip0 >=0 && e_to_v[edge][1] == vtip0) continue; // these are in the same tria (not special cases)
        if (vtip1 >=0 && e_to_v[edge][1] == vtip1) continue; // will be removed anyway later
        if (!simulate) e_to_v[edge][0] = vnew;
        vother = e_to_v[edge][1];
      }
      else
      {
        if (e_to_v[edge][0] == vnew) continue; // skip edge between vold, vnew (else it would collapse onto a point)
        if (vtip0 >=0 && e_to_v[edge][0] == vtip0) continue;
        if (vtip1 >=0 && e_to_v[edge][0] == vtip1) continue;
        if  (!simulate) e_to_v[edge][1] = vnew;
        vother = e_to_v[edge][0];
      }
      //ret++;
      assert(vother != vold);
      assert(vother != vnew);
      assert(vother == iter->pos);
        
      // earse this edge in v_to_e[vold]
      if (! simulate)
      {
        v_to_e.erase(vold,vother);
        v_to_e.erase(vother,vold);
      }
      
      if (!v_to_e.isSet(vother,vnew)) // edge did not exist before
      {
        // add it to v_to_e
        if (!simulate )
        {
          v_to_e.setVal(vother,vnew,edge+1);
          v_to_e.setVal(vnew,vother,edge+1);
        } //                                     x
      }  //                                     |||
      else //                                  |x |
      {  //                                   |/\ |
         //                                  |/  \|
         // special case, edge existed       x----x
         // keep the other one, and change this edge to the other in all trias
         // this makes everything non-manifold
         //cout << " vold: " << vold << " vother: " << vother << " vnew: " << vnew << endl;
         if (verbose > 1)
           cout << " special case, found edge : " << (int)v_to_e.getVal(vother,vnew)-1 << " between: " << vother << " <--> " << vnew << endl;
         if (! simulate) mergeEdges((int)v_to_e.getVal(vother,vnew)-1,edge);    
         ret++;
         //exit(1);  
      }
    }
    return ret;
}


/** Removes reference to edge from v_to_e matrix.
  Returns true if deleted.
  Returns false if not found.*/
bool Remesher::rmEdgeInVtoE(int eidx)
{
   assert(e_to_v[eidx].size() == 2); //should still have the two vertices
   int v0 = e_to_v[eidx][0];
   int v1 = e_to_v[eidx][1];
   
   // try to delete entries:
   int val1 = (int) v_to_e.eraseGetVal(v0,v1) - 1;
   int val2 = (int) v_to_e.eraseGetVal(v1,v0) - 1;
   assert (val1 == val2); // either both were allready empty or they should be the edge
   unused(val2);

   if (val1 == -1) return false; // was allready deleted
   
   assert(val1 == eidx);
   return true;
}


/** Deals with edge structures in this trias neighbors*/
void Remesher::contractEdgeInTria(int eidx, int tidx)
{
  assert(tria[tidx].size() == 3);
  
  // first sort out local indices
  unsigned int lv;
  for (lv =0;lv<3;lv++) 
    if (t_to_e[tidx][lv] == eidx) break;
  assert(lv < 3);
  // lv mark the begining of collapsing edge
  
  // get the two edges at the non-collapsing vertex
  int e0 = t_to_e[tidx][(lv+2)%3];
  int e1 = t_to_e[tidx][(lv+1)%3];
  assert (e0 != e1);
  
  // get the two neighbor triangles
  if (e_to_t[e0].size() != 2)
  {
    cout << "e_to_t[" << e0 <<"] size !=2 : " << e_to_t[e0] << endl;
    assert (e_to_t[e0].size() == 2); // boundary not allowed yet
  }
  int t0 = e_to_t[e0][0];
  int lei = 1; // local idx of tidx
  if ( t0 == tidx ) { t0 = e_to_t[e0][1]; lei = 0;}
  else assert( e_to_t[e0][lei] == tidx );
  assert (e_to_t[e1].size() == 2); // boundary not allowed yet
  int t1 = e_to_t[e1][0];
  if ( t1 == tidx )  t1 = e_to_t[e1][1];
  else assert( e_to_t[e1][1] == tidx );
  assert (t0 != t1);
  
  // replace reference to tidx in e0 with t1
  e_to_t[e0][lei] =  t1;
  // replace e1 in t1 with e0
  unsigned int i;
  for (i = 0;i<3;i++)
     if ( t_to_e[t1][i] == e1 ) { t_to_e[t1][i] = e0; break;}
  assert (i<3);
    
  // delete e1 from v_to_e
  rmEdgeInVtoE(e1);
  // delete e1 from e_to_t and e_to_v 
  e_to_t[e1].clear();
  e_to_v[e1].clear(); 
}


/** Remove triangles with identical vertices
 where two edges are the same and one
 edge is a double edge (same vertices, but two edges).
 Think of it as a triangular jelly bag.
 vidx is the tip, mididx is the vertex from the collapsed edge.
 Returns nextidx the other index in these trias, to be checked next.*/
int Remesher::removeCollapsedTrias(int vidx, int mididx)
{
   if (v_to_t[vidx].size() > 2) return -1;
   if ( v_to_t[vidx].size() == 1) cout << " tria: " << tria[v_to_t[vidx][0]] << endl;
   assert (v_to_t[vidx].size() == 2);
   
   // get two trias
   int t0  = v_to_t[vidx][0];
   int lv0 = v_to_lv[vidx][0];
   int t1  = v_to_t[vidx][1];
   int lv1 = v_to_lv[vidx][1];
   assert (tria[t0][lv0] == tria[t1][lv1]); 
   
   // check if all three vertices agree (note reverse order):
   if (tria[t1][(lv1+2)%3] != tria[t0][(lv0+1)%3]) return -1;
   if (tria[t1][(lv1+1)%3] != tria[t0][(lv0+2)%3]) return -1;
   
   cout << " Remesher::removeCollapsedTrias( "<<vidx<<" , "<<mididx<<" )" << endl;
   cout << " v_to_t[vidx] : " <<  v_to_t[vidx] << endl;
   cout << " Trias are identical : ... removing ... " << endl;
   
   
   // check edges at tip (should only have two trias each)
   assert (e_to_t[t_to_e[t0][lv0]].size() == 2);
   assert (e_to_t[t_to_e[t0][(lv0+2)%3]].size() == 2);
   assert (e_to_t[t_to_e[t1][lv1]].size() == 2);
   assert (e_to_t[t_to_e[t1][(lv1+2)%3]].size() == 2);
   // and should be the same edge
   assert (t_to_e[t0][lv0] == t_to_e[t1][(lv1+2)%3]);
   assert (t_to_e[t1][lv1] == t_to_e[t0][(lv0+2)%3]);


   // so we can delete these two trias
    
   // find edge not containing vidx
   int e0 = t_to_e[t0][(lv0+1)%3];
   int e1 = t_to_e[t1][(lv1+1)%3];
   assert (e_to_v[e0][0] == e_to_v[e1][0] || e_to_v[e0][0] == e_to_v[e1][1] );
   assert (e_to_v[e0][1] == e_to_v[e1][0] || e_to_v[e0][1] == e_to_v[e1][1] );
   assert (e_to_v[e0][0] == mididx || e_to_v[e0][1] == mididx );
   unused(e1);
   
   // find nextidx
   int nextidx = e_to_v[e0][0];
   if (nextidx == mididx) nextidx = e_to_v[e0][1];
   
   // merge edge (making this non manifold here)
   mergeEdges(e0,e1);
   // e1 is gone now
//          int cs = checkStructure();
//       cout << " merge edge check: " << cs << " manifold: " << ismanifold << endl << endl << endl;    
//       assert(cs <= 0 );         


   
   removeTria(t1); // here boundary structres might be destroyed
   removeTria(t0); // now should be OK again
   // can still be non manifold, if nextidx is tip of another collapsing tria
   
//    int cs = checkStructure();
//      cout << "rmTria check: " << cs << " manifold: " << ismanifold << endl << endl << endl;    
//      assert(cs <= 0 && ismanifold );         
   
   cout << " NEXTIDX : " << nextidx << endl;
   return nextidx;
}


/** Merges identical edges by adding triangles at edel to etarg.
  Then deletes edel.
 Both edges need to point at the same vertices
 (onboundary is not correctly adjusted for vertices)
 It can happen (if each edge has more than 1 tria) that 
 this destroys manifold property.*/
bool Remesher::mergeEdges(unsigned int etarg, unsigned int edel)
{
   //cout << "Remesher::mergeEdges( " <<etarg << " , " << edel << " )" << endl;
   //cout << " etarg trias: " << e_to_t[etarg] << endl;
   //cout << " edel  trias: " << e_to_t[edel] << endl;
   if (etarg == edel) return true;
   
   if (! ( (e_to_v[etarg][0] == e_to_v[edel][0] &&
            e_to_v[etarg][1] == e_to_v[edel][1]) ||
           (e_to_v[etarg][0] == e_to_v[edel][1] &&
            e_to_v[etarg][1] == e_to_v[edel][0])))
   {
     cout << "Remesher::mergeEdges( " <<etarg << " , " << edel << " )" << endl;
     cout << " e_to_v[" << etarg << "] : " << e_to_v[etarg] << endl;
     cout << " e_to_v[" << edel  << "] : " << e_to_v[edel] << endl;
     assert( (e_to_v[etarg][0] == e_to_v[edel][0] &&
              e_to_v[etarg][1] == e_to_v[edel][1]) ||
             (e_to_v[etarg][0] == e_to_v[edel][1] &&
              e_to_v[etarg][1] == e_to_v[edel][0]));
    }        
   // adjust bedgecount
   if (e_to_t[edel].size() == 1) bedgecount--;
   if (e_to_t[etarg].size() == 1) bedgecount--;
   
   // in all trias at edel, replace edel with etarg
   for (unsigned int t = 0;t<e_to_t[edel].size();t++)
   {
      bool found = false;
      int tidx = e_to_t[edel][t];
      //cout << " before tria " << tidx << " has these edges " << t_to_e[tidx] << endl;
      for (unsigned int i = 0;i<3;i++)
         if (t_to_e[tidx][i] == (int)edel)
         {
            found = true;
            t_to_e[tidx][i] = etarg;
         }
      //cout << " after tria " << tidx << " has these edges " << t_to_e[tidx] << endl;
      assert (found);
   }
      
   // append trias from edel to etarg
   for (unsigned int td = 0;td < e_to_t[edel].size(); td++)
   {
      bool found = false;
      for (unsigned int tt = 0;tt < e_to_t[etarg].size();tt++)
        if (e_to_t[etarg][tt] == e_to_t[edel][td])
        {
           found = true;
           break;
        }
      if (! found) e_to_t[etarg].push_back(e_to_t[edel][td]);
   
   }
   
   // ensure v_to_e holds the remaining edge (etarg)
   // it can only hold one edge, so it might have been the wrong one
   // not sure if this gives any difficulties before (in collapseedge)??
   v_to_e.setVal(e_to_v[etarg][0],e_to_v[etarg][1],etarg+1);
   v_to_e.setVal(e_to_v[etarg][1],e_to_v[etarg][0],etarg+1);
   
   // free edel
   e_to_t[edel].clear();
   e_to_v[edel].clear(); // v_to_e is allready set correctly above

   return true;
}


void SparseMatrix::redim(int m, int n)
{
  assert(m == n);
  if (row == m && col == n) return;
  entries.resize(n);
  row = col = n;
}


// resets matrix keeping the dimension
void SparseMatrix::fillzero()
{
  for (unsigned int i=0; i< entries.size(); i++) entries[i].clear();
  nonzero = 0;
}


// removes entry i, j
bool SparseMatrix::erase(int i, int j)
{
  assert(i >= 0 && i < row && j >= 0 && j < col);
  bool ret = false;
  std::list<entry>::iterator iter;
  for (iter = entries[i].begin(); iter != entries[i].end(); iter++) {
    if(iter->pos == j) {
      ret = true;
      entries[i].erase(iter);
      nonzero--;
      break;
    }
  }
  return ret;
}


double SparseMatrix::getVal(int i, int j) const
{
  assert(i>=0 && i<row && j>=0 && j<col);
  std::list<entry>::const_iterator iter;
  for (iter= entries[i].begin(); iter != entries[i].end(); iter++) {
    if((*iter).pos == j)  {
      return (*iter).value;
    } else if ((*iter).pos > j) {
      return 0;
    }
  }
  return 0;
}


void SparseMatrix::setVal(int i, int j, double d)
{
  if (i >= row) {
    row = col = i+1;
    entries.resize(i+1);
  }

  std::list<entry>::iterator iter;
  for (iter= entries[i].begin(); iter != entries[i].end() && (*iter).pos < j; iter++) {}

  if (iter != entries[i].end() && (*iter).pos == j) {
   (*iter).value = d;
  } else {
    entry neu;
    neu.pos = j; 
    neu.value = d;
    entries[i].insert(iter, neu);
    nonzero++;
  }    
}


const double& SparseMatrix::getRef(int i, int j) const
{
  assert(i >= 0 && i < row && j >= 0 && j < col);
  std::list<entry>::const_iterator iter;
  for (iter=entries[i].begin(); iter != entries[i].end(); iter++) {
    if((*iter).pos == j) return (*iter).value;
  }

  std::cerr << "const SSMatrix::getRef( " << i << " , " << j << " ) does not exist!" << std::endl;
  exit(1);
}


double& SparseMatrix::getRef(int i, int j)
{
  if(i < 0 || i >= row || j < 0 || j >= col) {
    std::cout << "SparseMatrix::getRef( " << i << " , " << j << " ) with row: " << row << "  col: " << col << std::endl;
    assert(i >= 0);
    assert(i < row);
    assert(j >= 0);
    assert(j < col);
    exit(1);
  }

  std::list<entry>::iterator iter;
  for (iter= entries[i].begin(); iter != entries[i].end() && (*iter).pos < j ; iter++) {};
  if(iter == entries[i].end() || (*iter).pos != j) {
    entry neu;
    neu.pos = j; 
    neu.value = 0;
    entries[i].insert(iter, neu);
    iter--;
    nonzero++;
  }

  return (*iter).value;
}


double SparseMatrix::eraseGetVal(int i, int j)
{
  assert(i>=0 && i<row && j>=0 && j<col);
  double ret = 0;
  std::list<entry>::iterator iter;
  for (iter=entries[i].begin(); iter != entries[i].end(); iter++) {
    if(iter->pos == j) {
      ret = iter->value;
      entries[i].erase(iter);
      nonzero--;
      break;
    }
  }
  return ret;
}


bool SparseMatrix::isSet(int i,int j) const
{
  assert(i>=0 && i<row && j>=0 && j<col);
  std::list<entry>::const_iterator iter;
  for (iter= entries[i].begin(); iter != entries[i].end(); iter++) {
    if((*iter).pos == j) {
      return true;
    } else if ((*iter).pos > j) {
      return false;
    }
  }

  return false;
}


double SparseMatrix::getMax() const
{
  double max = getVal(0,0);
  for (int i=0; i<row; i++) {
    // Iterate through list
    list<entry>::const_iterator iter;
    for (iter=entries[i].begin(); iter != entries[i].end(); iter++) {
      if ( iter->value > max) max = iter->value;
    }
  }
  return max;
}


double SparseMatrix::getMin() const
{
  double min = getVal(0,0);
  for (int i=0; i<row; i++) {
    // Iterate through list
    list<entry>::const_iterator iter;
    for (iter=entries[i].begin(); iter != entries[i].end(); iter++) {
      if ( iter->value < min) min = iter->value;
    }
  }
  return min;
}
