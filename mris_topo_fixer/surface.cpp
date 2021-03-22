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


#include <fstream>

#include "globals.h"
#include "surface.h"
#include "patchdisk.h"
#include "fastloop.h"

Surface::Surface(void) {
  nvertices = maxvertices= 0;
  vertices = NULL;
  nfaces = maxfaces = 0 ;
  faces = NULL;
  type_of_surface = UNKNOWN_TYPE_OF_SURFACE;
  vtrans_to = ftrans_to = vtrans_from = ftrans_from = NULL;
  surface_source = NULL;
  disk=0;
}

Surface::Surface(int nv, int nf) {
  _Allocate(nv,nf);
  type_of_surface = UNKNOWN_TYPE_OF_SURFACE;
  vtrans_to = ftrans_to = vtrans_from = ftrans_from = NULL;
  surface_source = NULL;
  disk=0;
}

Surface::Surface(const string s) {
  nvertices = maxvertices= 0;
  vertices = NULL;
  nfaces = maxfaces = 0 ;
  faces = NULL;
  type_of_surface = UNKNOWN_TYPE_OF_SURFACE;
  vtrans_to = ftrans_to = vtrans_from = ftrans_from = NULL;
  surface_source = NULL;
  disk=0;
  OpenFile(s);
}


Surface *Surface::Clone() const {
  Surface *surf=new Surface(maxvertices,maxfaces);

  surf->nvertices=nvertices;
  surf->nfaces=nfaces;
  surf->euler=euler;
  for (int n = 0 ; n < maxvertices ; n++)
    surf->vertices[n]=vertices[n];
  for (int n = 0 ; n < maxfaces ; n++)
    surf->faces[n]=faces[n];

  surf->disk = disk;

  return surf;
}

Surface::~Surface(void) {
  if (vertices) delete [] vertices;
  if (faces) delete [] faces;
  if (vtrans_to) delete [] vtrans_to;
  if (ftrans_to) delete [] ftrans_to;
  if (vtrans_from) delete [] vtrans_from;
  if (ftrans_from) delete [] ftrans_from;
}

int Surface::_Allocate(int nv, int nf) {
  nvertices = 0 ;
  maxvertices = nv;
  vertices = new Vertex[maxvertices];
  nfaces = 0;
  maxfaces = nf ;
  faces = new Face[maxfaces];

  return 0;
}

void Surface::_OverAlloc(int nextrav,int nextraf) {
  if (nvertices + nextrav > maxvertices) {
    int new_maxvertices = nvertices + nextrav;
    Vertex *new_vertices = new Vertex[new_maxvertices];
    for (int n = 0 ; n < maxvertices ; n++) {
      new_vertices[n] = vertices[n];
      vertices[n].Clear();
    }
    maxvertices = new_maxvertices;
    if (vertices) delete [] vertices;
    vertices = new_vertices;
  };
  if (nfaces + nextraf > maxfaces) {
    int new_maxfaces = nfaces + nextraf;
    Face *new_faces = new Face[new_maxfaces];
    for (int n = 0 ; n < maxfaces ; n++)
      new_faces[n]=faces[n];
    maxfaces = new_maxfaces;
    if (faces) delete [] faces;
    faces = new_faces;
  }
}

void Surface::Expand(int nextrav,int nextraf) {
  _OverAlloc(nextrav,nextraf);
}

int Surface::OpenFile(const string s, int verbose) {
  ifstream file(s.c_str());

  if (s[s.size()-1]=='d') { // CERTIS file format
    file >> nvertices;
    maxvertices = nvertices ;
    if (verbose) cout << nvertices << " vertices - " ;
    vertices = new Vertex[maxvertices];
    for (int n= 0 ; n < nvertices ; n++)
      file >> vertices[n].x >> vertices[n].y >> vertices[n].z ;
    file >> nfaces;
    maxfaces = nfaces ;
    if (verbose) cout << nfaces << " faces " << endl;
    faces = new Face[maxfaces];
    for (int n = 0 ; n < nfaces ; n++)
      file >> faces[n].v[0] >> faces[n].v[1] >>faces[n].v[2] ;
  } else if (s[s.size()-1]=='c') { // ASCII MGH file format
    int tmp;
    file >> nvertices;
    maxvertices = nvertices ;
    if (verbose) cout << nvertices << " vertices - " ;
    vertices = new Vertex[maxvertices];
    for (int n= 0 ; n < nvertices ; n++)
      file >> vertices[n].x >> vertices[n].y >> vertices[n].z >> tmp;

    file >> nfaces;
    maxfaces = nfaces ;
    if (verbose) cout << nfaces << " faces " << endl;
    faces = new Face[maxfaces];
    for (int n = 0 ; n < nfaces ; n++)
      file >> faces[n].v[0] >> faces[n].v[1] >>faces[n].v[2] >> tmp;
  }

  InitSurface();

  return 0;
}

int Surface::WriteFile(const string s, int verbose) const {
  ofstream file(s.c_str());

  file << nvertices << " " << nfaces << endl;
  for (int n= 0 ; n < nvertices ; n++)
    file << vertices[n].x << " " << vertices[n].y << " " << vertices[n].z << " " << 0 << endl;
  for (int n = 0 ; n < nfaces ; n++)
    file << faces[n].v[0] << " " << faces[n].v[1] << " " << faces[n].v[2] << " " << 0 << endl ;
  if (verbose) cout << endl << "surface written in " << s << " " ;
  return 0;
}

int Surface::GetDefectLabels(const string s) {
  OpenCurvatureFile(s);
  for (int n = 0 ; n < nvertices ; n++)
    vertices[n].marked = int(vertices[n].curv);
  return 0 ;
}

int Surface::OpenCurvatureFile(const string s) {

  ifstream file(s.c_str());

  int tmp;
  float a,b,c,curv;
  for (int n = 0 ; n < nvertices ; n++) {
    file >> tmp >> a >> b >> c >> curv;
    vertices[n].curv = curv;
  }
  return 0;
}

int Surface::InitSurface(void) {
  _InitSurfaceConnectivity();
  if (_InitFaceConnectivity())
    type_of_surface = CLOSED_SURFACE;
  else
    type_of_surface = OPEN_SURFACE;
  _InitFaceCoordinates();
  return 0;
}


int Surface::_InitSurfaceConnectivity(void) {
  Vertex *v;

  for (int n = 0 ; n < nvertices ; n++) {
    v=&vertices[n];
    v->fnum=0;
    v->marked=0;
  }

  // counting the number of faces per vertex
  for (int n = 0 ; n < nfaces ; n++)
    for (int i = 0 ; i < 3 ; i++)
      vertices[faces[n].v[i]].fnum++;
  // allocate the list of faces
  for (int n = 0 ; n < nvertices ; n++)
    vertices[n].AllocateFaces(vertices[n].fnum);
  // initialize the list of faces
  for (int n = 0 ; n < nfaces ; n++) {
    for (int i = 0 ; i < 3 ; i++) {
      int vno = faces[n].v[i];
      v=&vertices[vno];
      v->f[v->fnum]=n;
      v->n[v->fnum++]=i;
    }
  }

  // counting the list of vertices
  for (int n = 0 ; n < nvertices ; n++) {
    v=&vertices[n];
    v->vnum=0;
    for (int p = 0 ; p < v->fnum ; p++) {
      Face *face=&faces[v->f[p]];
      for ( int i = 0 ; i < 3 ; i++) {
        int vn=face->v[i];
        if (vn==n) continue;
        if (vertices[vn].marked) continue;
        vertices[vn].marked=1;
        v->vnum++;
      }
    }
    // allocate the list of vertices
    vertices[n].AllocateVertices(v->vnum);
    for (int p = 0 ; p < v->fnum ; p++) {
      Face *face=&faces[v->f[p]];
      for ( int i = 0 ; i < 3 ; i++) {
        int vn=face->v[i];
        if (vn==n) continue;
        if (vertices[vn].marked == 0 ) continue;
        vertices[vn].marked=0;
        v->v[v->vnum++]=vn;
      }
    }
  }

  return 0;
}

// this function assumes that each face is a triangle!
bool Surface::_InitFaceConnectivity(void) {
  bool is_closed = true;
  for (int i = 0 ; i < nfaces ; i++) {
    Face *face=&faces[i];
    for (int j = 0 ; j < 3 ; j++) {
      int n1 = face->v[j] , n2;
      //n2 = face->v[(j+1)%2];
      if (j==2) n2=face->v[0];
      else n2=face->v[j+1];

      int fn = _FindFace(n1,n2,i);
      face->f[j]=fn;
      if (fn == -1) is_closed=false;
    }
  }

  return is_closed;
}

// find the neighboring face of fn, which contains n1 and n2
// return -1 if fn is a border face
int Surface::_FindFace(int n1, int n2, int fn) const {
  Vertex *v1=&vertices[n1],*v2=&vertices[n2];
  for (int k = 0 ; k < v1->fnum ; k++) {
    if (v1->f[k]==fn) continue;
    for (int p = 0 ; p < v2->fnum ; p++)
      if (v2->f[p]==v1->f[k]) return v1->f[k];
  }

  return -1;
}

int Surface::_InitFaceCoordinates(void) {
  for (int k = 0 ; k < nfaces ; k++) {
    double x,y,z;
    x=y=z=0.0f;
    for (int i = 0 ; i < 3 ; i++) {
      x+=vertices[faces[k].v[i]].x;
      y+=vertices[faces[k].v[i]].y;
      z+=vertices[faces[k].v[i]].z;
    }
    faces[k].x=x/3.0;
    faces[k].y=y/3.0;
    faces[k].z=z/3.0;
  }
  return 0 ;
}

void Surface::Center(void) {
  Vertex *v;
  double x,y,z;
  //compute the center of gravity
  x=y=z=0.0;
  for (int n = 0 ; n < nvertices ; n++) {
    v = &vertices[n];
    x += v->x;
    y += v->y;
    z += v->z;
  }
  x /= double(nvertices);
  y /= double(nvertices);
  z /= double(nvertices);

  //center the surface
  for (int n = 0 ; n < nvertices ; n++) {
    v = &vertices[n];
    v->x -= x;
    v->y -= y;
    v->z -= z;
  }

}

void Surface::scale(double scaling_factor) {
  double r=0;
  for (int n = 0 ; n < nvertices ; n++) {
    Vertex *v=&vertices[n];
    r += __norm(v->x,v->y,v->z);
  }
  r /= nvertices;
  r = scaling_factor/r;
  for (int n = 0 ; n < nvertices ; n++) {
    Vertex *v=&vertices[n];
    v->x *= r;
    v->y *= r;
    v->z *= r;
  }
}


int Surface::GetEuler(int mark) {
  int nv,ne,nf;
  return GetEuler(nv,ne,nf,mark);
}

int Surface::GetEuler(int &nv, int &ne, int &nf, int mark) {
  if (mark < 0 ) { //Euler characteristic of the whole surface
    nedges = 0 ;
    for (int n = 0 ; n < nvertices ; n++)
      nedges += vertices[n].vnum ;
    nedges /= 2 ;
    euler = nvertices + nfaces - nedges;
    nv = nvertices;
    ne = nedges;
    nf = nfaces;
    return euler;
  };

  //counting number of vertices
  nv = 0 ;
  for (int n = 0 ; n < nvertices ; n++) {
    Vertex *v=&vertices[n];
    if (v->marked!=mark) continue;
    nv++;
    for (int p = 0 ; p < v->vnum ; p++)
      v->e[p]=0;
  }
  //counting number of faces
  nf = 0 ;
  for (int n = 0 ; n < nfaces ; n++) {
    Face *face = &faces[n];
    bool is_face=true;
    for (int i = 0 ; i < 3 ; i++) {
      if (vertices[face->v[i]].marked!=mark) {
        is_face=false;
        break;
      }
    }
    if (is_face) {
      nf++;
      //add the edges
      for (int i = 0 ; i < 3 ; i++) {
        Vertex *v = &vertices[face->v[i]];
        int v1,v2;
        (i==0) ? v1 = face->v[2] : v1 = face->v[i-1];
        (i==2) ? v2 = face->v[0] : v2 = face->v[i+1];
        for (int p = 0 ; p < v->vnum ; p++) {
          if (v->v[p] == v1 || v->v[p] == v2)
            v->e[p]=1;
        }
      }
    }
  }
  //counting number of edges
  ne = 0 ;
  for (int n = 0 ; n < nvertices ; n++) {
    Vertex *v=&vertices[n];
    if (v->marked!=mark) continue;
    for (int p = 0 ; p < v->vnum; p++) {
      if (v->e[p]) ne++;
    }
  }
  ne = ne/2;

  return (nv-ne+nf);
}

int Surface::GetEuler(const int *list_of_faces, int nfs) {
  int nv,nf,ne;

  //mark and count the vertices
  for (int n = 0 ; n < nvertices ; n++)
    vertices[n].marked = 0;
  nv = 0 ;
  for (int n = 0 ; n < nfs ; n++)
    for (int i = 0 ; i < 3 ; i++) {
      Vertex *v=&vertices[faces[list_of_faces[n]].v[i]];
      if (v->marked==0) {
        nv++;
        v->marked=1;
        for (int p = 0 ; p < v->vnum ; p++)
          v->e[p]=0;
      };
    }
  //mark and count the faces
  for (int n = 0 ; n < nfaces ; n++)
    faces[n].marked=0;
  nf=ne=0;
  for (int n = 0 ; n < nfs ; n++) {
    Face *face= &faces[list_of_faces[n]];
    if (face->marked==0) {
      nf++;
      face->marked=1;
      int vn0,vn1;
      Vertex *v0,*v1;
      for (int i = 0 ; i < 3 ; i++) {
        vn0 = face->v[i];
        if (i==2) vn1 = face->v[0];
        else vn1 = face->v[i+1];
        v0 = &vertices[vn0];
        v1 = &vertices[vn1];
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

  return (nv-ne+nf);
}

void Surface::PrintDefectInfo(int ndefect) {
  if (ndefect<0) {
    int nlabels = 0 ;
    for (int n = 0 ; n < nvertices ; n++)
      if (vertices[n].marked > nlabels)
        nlabels = vertices[n].marked;
    cout << nlabels-1 << " defects in the surface" << endl;

    for (int l = 1 ; l < nlabels ; l++)
      PrintDefectInfo(l);
  } else {
    cout << "   Defect #" << ndefect << ": ";
    int nv,ne,nf,eul;
    eul = GetEuler(nv,ne,nf,ndefect);
    cout << " (" << nv << "," << ne << "," << nf <<") -> X = " << eul << endl;
  }
}

void Surface::SetMarks(int mark) {
  for (int n = 0 ; n < nvertices ; n++)
    vertices[n].marked=mark;
}

void Surface::SetMarks(const int *v, int nv, int mark) {
  for (int n = 0 ; n < nv ; n++)
    vertices[v[n]].marked=mark;
}

void Surface::ExpandMarks(int niters,int mark) {

  while (niters--) {
    for (int n = 0 ; n < nvertices ; n++)
      vertices[n].curv = 0.0;
    for (int n = 0 ; n < nvertices ; n++)
      if (vertices[n].marked==mark) {
        for (int p = 0 ; p < vertices[n].vnum ; p++)
          vertices[vertices[n].v[p]].curv = mark;
      };
    for (int n = 0 ; n < nvertices ; n++)
      if (vertices[n].curv == double(mark))
        vertices[n].marked = mark;
  }

}


void Surface::Smooth(int niters, const int *tab,int nv) {
  double alpha = 0.3;
  while (niters--) {
    for (int n = 0 ; n < nv ; n++) {
      double x=0.,y=0.,z=0.;
      Vertex *v = &vertices[tab[n]];
      for (int p = 0 ; p < v->vnum ; p++) {
        Vertex *vp = &vertices[v->v[p]];
        x += vp->x;
        y += vp->y;
        z += vp->z;
      }
      v->tx = (1 - alpha) * v->x + alpha * x/(double)v->vnum;
      v->ty = (1 - alpha) * v->y + alpha * y/(double)v->vnum;
      v->tz = (1 - alpha) * v->z + alpha * z/(double)v->vnum;
    }
    for (int n = 0 ; n < nv ; n++) {
      Vertex *v = &vertices[tab[n]];
      v->x = v->tx;
      v->y = v->ty;
      v->z = v->tz;
    }
  }
}

void Surface::SmoothMarked(int niters, int mark) {
  double alpha = 0.3;
  while (niters--) {
#if PRINT_MODE
    cout << ".";
#endif
    for (int n = 0 ; n < nvertices ; n++) {
      double x=0.,y=0.,z=0.;
      Vertex *v = &vertices[n];
      if (v->marked!=mark) continue;
      for (int p = 0 ; p < v->vnum ; p++) {
        Vertex *vp = &vertices[v->v[p]];
        x += vp->x;
        y += vp->y;
        z += vp->z;
      }
      v->tx = (1 - alpha) * v->x + alpha * x/(double)v->vnum;
      v->ty = (1 - alpha) * v->y + alpha * y/(double)v->vnum;
      v->tz = (1 - alpha) * v->z + alpha * z/(double)v->vnum;
    }
    for (int n = 0 ; n < nvertices ; n++) {
      Vertex *v = &vertices[n];
      if (v->marked!=mark) continue;
      v->x = v->tx;
      v->y = v->ty;
      v->z = v->tz;
    }
  }
}


void Surface::Smooth(int niters) {
  double alpha = 0.3;
  while (niters--) {
#if PRINT_MODE
    cout << ".";
#endif
    for (int n = 0 ; n < nvertices ; n++) {
      double x=0.,y=0.,z=0.;
      Vertex *v = &vertices[n];
      for (int p = 0 ; p < v->vnum ; p++) {
        Vertex *vp = &vertices[v->v[p]];
        x += vp->x;
        y += vp->y;
        z += vp->z;
      }
      v->tx = (1 - alpha) * v->x + alpha * x/(double)v->vnum;
      v->ty = (1 - alpha) * v->y + alpha * y/(double)v->vnum;
      v->tz = (1 - alpha) * v->z + alpha * z/(double)v->vnum;
    }
    for (int n = 0 ; n < nvertices ; n++) {
      Vertex *v = &vertices[n];
      v->x = v->tx;
      v->y = v->ty;
      v->z = v->tz;
    }
  }
}

Surface *Surface::ExtractPatch(int mark, int nextravertices, int nextrafaces) {
  int nv , nf;

  //counting the number of vertices
  nv = 0 ;
  for (int n = 0 ; n < nvertices ; n++)
    if (vertices[n].marked == mark)
      nv++;
  if (nv==0) return NULL;

  //counting the number of faces
  nf = 0 ;
  for (int n = 0 ; n < nfaces ; n++) {
    bool is_face = true;
    for (int i = 0 ; i < 3 ; i++)
      if (vertices[faces[n].v[i]].marked != mark) {
        is_face = false;
        break;
      };
    if (is_face) nf++;
  }
  /////////////////////////////////////////////////////////////////
  // allocate the new surface
  Surface *surface = new Surface(nv+nextravertices,nf+nextrafaces);
  surface->nvertices=nv;
  surface->nfaces=nf;

  //the corresponding tables
  int *vt_to,*vt_from;
  vt_to = new int[nv];
  vt_from = new int[nvertices];

  nv=0;
  for (int n = 0 ; n < nvertices ; n++)
    if (vertices[n].marked == mark) {
      Vertex *vsrc = &vertices[n];
      Vertex *v = &surface->vertices[nv];
      vt_from[n]=nv;
      vt_to[nv++]=n;
      //copy the strict necessary
      v->x = vsrc->x;
      v->y = vsrc->y;
      v->z = vsrc->z;
    };
  surface->vtrans_to = vt_to;
  surface->vtrans_from = vt_from;

  //the corresponding tables
  int *ft_to,*ft_from;
  ft_to = new int[nf];
  ft_from = new int[nfaces];

  nf=0;
  for (int n = 0 ; n < nfaces ; n++) {
    bool is_face = true;
    for (int i = 0 ; i < 3 ; i++)
      if (vertices[faces[n].v[i]].marked != mark) {
        is_face = false;
        break;
      };
    if (is_face) {
      Face *fsrc = &faces[n];
      Face *f = &surface->faces[nf];
      ft_from[n]=nf;
      ft_to[nf++]=n;
      for (int i = 0 ; i < 3 ; i++)
        f->v[i] = vt_from[fsrc->v[i]];
    }
  }
  surface->ftrans_to = ft_to;
  surface->ftrans_from = ft_from;

  surface->InitSurface();
  surface->type_of_surface = PATCH;

  surface->disk = disk;

  return surface;
}

//testing the validity of the loop
bool Surface::LoopValid(Loop &loop) {

  for (int n = 0 ; n < maxfaces ; n++)
    faces[n].marked=0;

  for (int n = 0 ; n < loop.npoints ; n++)
    faces[loop.points[n]].marked=1;

  for (int n = 0 ; n < loop.npoints ; n++) {
    int nf = 0 ;
    for (int i = 0 ; i < 3 ; i++) {
      if (faces[loop.points[n]].f[i] == -1) continue;
      if (faces[faces[loop.points[n]].f[i]].marked==1) nf++;
    }
    if (nf!=2)  return false;
  }
  return true;
}

void Surface::KnitPatch(Loop &loop , PatchDisk *pdisk) {

  pdisk->Init(); //reset the initial values of the ring

  double x=0.,y=0.,z=0.;
  for (int n = 0 ; n < loop.npoints ; n++) {
    x += vertices[loop.points[n]].x;
    y += vertices[loop.points[n]].y;
    z += vertices[loop.points[n]].z;
  }
  x /= loop.npoints;
  y /= loop.npoints;
  z /= loop.npoints;

  //adding vertices
  for (int n = 0 ; n < pdisk->disk.nvertices ; n++) {
    Vertex *vdst = &vertices[nvertices];
    Vertex *vsrc = &pdisk->disk.vertices[n];
    pdisk->vtrans[n]=nvertices;
    if (vsrc->marked==2) pdisk->ring.Replace(n,nvertices);
    nvertices++;
    vdst->AllocateFaces(vsrc->fnum);
    vdst->AllocateVertices(vsrc->vnum);
    vdst->x = x ;
    vdst->y = y ;
    vdst->z = z;
    vdst->marked=4; //marking vertices with 4
  }
  //adding faces
  for (int n = 0 ; n < pdisk->disk.nfaces ; n++) {
    Face *fdst = &faces[nfaces];
    Face *fsrc = &pdisk->disk.faces[n];
    pdisk->ftrans[n]=nfaces++;
    fdst->v[0] = pdisk->vtrans[fsrc->v[0]];
    fdst->v[1] = pdisk->vtrans[fsrc->v[1]];
    fdst->v[2] = pdisk->vtrans[fsrc->v[2]];
  }
  // adding connectivity stuff
  for (int n = 0 ; n < pdisk->disk.nvertices ; n++) {
    Vertex *vdst = &vertices[pdisk->vtrans[n]];
    Vertex *vsrc = &pdisk->disk.vertices[n];
    vdst->fnum = vsrc->fnum;
    for (int p = 0 ; p < vsrc->fnum ; p++) {
      vdst->f[p] = pdisk->ftrans[vsrc->f[p]];
      vdst->n[p] = vsrc->n[p];
    }
    vdst->vnum = vsrc->vnum;
    for (int p = 0 ; p < vsrc->vnum ; p++) {
      vdst->v[p] = pdisk->vtrans[vsrc->v[p]];
      vdst->e[p] = 0;
    }
  }

  int pos1=0,pos2=0;
  while (pos1<loop.npoints || pos2 < pdisk->ring.npoints) {
    int v0,v1,v2;
    if (pos1 == loop.npoints) { //add the triangle pos1,pos2,pos2+1
      v0=loop[0];
      v1=pdisk->ring[pos2];
      v2=pdisk->ring[(pos2+1)%pdisk->ring.npoints];
      pos2++;
    } else if (pos2 == pdisk->ring.npoints) { //add the triangle pos1,pos1+1,pos2
      v1=loop[pos1];
      v2=pdisk->ring[0];
      v0=loop[(pos1+1)%loop.npoints];
      pos1++;
    } else { //determine which triangle to add
      if (fabs(pdisk->ring.npoints*pos1-loop.npoints*(pos2+1.)) < fabs(pdisk->ring.npoints*(pos1+1.)-loop.npoints*pos2)) {
        v0=loop[pos1];
        v1=pdisk->ring[pos2];
        v2=pdisk->ring[(pos2+1)%pdisk->ring.npoints];
        pos2++;
      } else {
        v1=loop[pos1];
        v2=pdisk->ring[pos2];
        v0=loop[(pos1+1)%loop.npoints];
        pos1++;
      }
    }
    //add the face v0,v1,v2
    Face *fnew = &faces[nfaces];
    fnew->v[0]=v2;
    fnew->v[1]=v1;
    fnew->v[2]=v0;
    vertices[v0].AddFace(nfaces,2);
    vertices[v1].AddFace(nfaces,1);
    vertices[v2].AddFace(nfaces,0);
    nfaces++;
    //update  v0 <--> v2
    vertices[v0].AddEdge(v2);
    vertices[v2].AddEdge(v0);
  }
}

double Surface::_FaceDistance(int fdst, int fsrc) {
  return __norm(faces[fdst].x-faces[fsrc].x,faces[fdst].y-faces[fsrc].y,faces[fdst].z-faces[fsrc].z);
}


double Surface::GetLoopLength(Loop &loop) {
  double length=0;
  for (int n = 0 ; n < loop.npoints ; n++) {
    int f0 = loop.points[n],f1;
    if (n==loop.npoints-1) f1=loop.points[0];
    else f1=loop.points[n+1];
    length += _FaceDistance(f0,f1);
  }
  return length;
}


void Surface::CutLoop(Loop &loop, int very_small_patch) {

  int first_v1 = -1, first_v2 = -1,  first_v3 = -1, first_v4 = -1;

  if (loop.npoints == 0) return;
  ASSERT(LoopValid(loop));

  //which patch should we use ?
  int wd = 0 ;
	if(very_small_patch==1){
		if (loop.npoints < 40 ) wd = 0;//6
		else if (loop.npoints < 80 ) wd = 1;//10
		else if (loop.npoints < 120) wd = 2;//14
		else wd = 3;
	}else{
		if (loop.npoints < 10 ) wd = 0;//6
		else if (loop.npoints < 14 ) wd = 1;//10
		else if (loop.npoints <  20) wd = 2;//14
		else wd = 3;
	}
  PatchDisk *pdisk = &disk[wd];
#if PRINT_ERROR_MODE
  cout << "loop is " << wd << " with " << GetLoopLength(loop) << endl;
#endif

  int *vertex_list = new int[loop.npoints*3];
  int nvertex_list=0;

  for (int n = 0 ; n < maxfaces ; n++)
    faces[n].marked=0;
  for (int n = 0 ; n < loop.npoints ; n++)
    faces[loop.points[n]].marked=1;

  for (int n = 0 ; n < maxvertices ; n++) {
    Vertex *v=&vertices[n];
    v->marked=0; //the new vertices are marked with 2
    for (int p = 0 ; p < v->vnum ; p++)
      v->e[p]=-1;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //verify that we have enough space in the surface
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 2*loop.npoints extra vertices and loop.npoints extra faces for the cut
  // 2*pdisk->disk.nvertices extra vertices for the 2 extra patches
  // 2*pdisk->disk.nfaces extra faces for the 2 patches
  // 2*(loop.npoints + pdisk->ring.npoints) to knit the 2 patches
  _OverAlloc(2*loop.npoints+2*pdisk->disk.nvertices,2*(2*loop.npoints+pdisk->disk.nfaces+pdisk->ring.npoints));

  //add 3 triangles (but reusing one triangle) & 4 points for each face of the defect
  for (int n = 0 ; n < loop.npoints ; n++) {
    int fn = loop.points[n];
    Face * face = &faces[fn];
    //find the vertices
    int vn0,vn1,vn2;
    int k=0;

    for (int i = 0 ; i < 3 ; i++) {
      if (face->f[i]==-1 || faces[face->f[i]].marked==0) {
        k=i;
        break;
      }
    }
    vn0=face->v[(k+2)%3];
    vn1=face->v[(k)%3];
    vn2=face->v[(k+1)%3];
    vertex_list[nvertex_list++]=vn0; //potentially with repeatitions
    vertex_list[nvertex_list++]=vn1;
    vertex_list[nvertex_list++]=vn2;

    //          vn2 <-----------> vn1
    //            \               /
    //             \             /
    //             v4 --------- v3
    //
    //              v2 ------ v1
    //                \      /
    //                 \    /    <- new face = current face
    //                  \  /
    //                  vn0

    //generate four new vertices
    int v1,v2,v3,v4;
    Vertex *vnew;
    //and 3 new faces
    Face *fnew;
    //intermediate variables
    double a;
    int va,vb;

    //v1 = 2/3*vn0+1/3*vn1
    a=2./3.;
    va=vn0;
    vb = vn1;
    v1 = -1;
    for (int p = 0 ; p < vertices[va].vnum ; p++)
      if (vertices[va].v[p] == vb) {
        v1=vertices[va].e[p];
        break;
      };
    if (v1 == -1) {
      v1=nvertices;
      vnew=&vertices[nvertices++];
      vnew->x = a*vertices[va].x + (1.-a)*vertices[vb].x;
      vnew->y = a*vertices[va].y + (1.-a)*vertices[vb].y;
      vnew->z = a*vertices[va].z + (1.-a)*vertices[vb].z;
      vnew->AllocateFaces(6);
      vnew->AllocateVertices(6);
      vnew->v[vnew->vnum]=va;
      vnew->e[vnew->vnum++] = -1;
      vnew->marked=1;
      for (int p = 0 ; p < vertices[va].vnum ; p++)
        if (vertices[va].v[p] == vb) {
          vertices[va].e[p] = v1;
          break;
        };
    }

    //v2 = 2/3*vn0+1/3*vn2
    a=2./3.;
    va=vn0;
    vb = vn2;
    v2 = -1;
    for (int p = 0 ; p < vertices[va].vnum ; p++)
      if (vertices[va].v[p] == vb) {
        v2=vertices[va].e[p];
        break;
      };
    if (v2 == -1) {
      v2=nvertices;
      vnew=&vertices[nvertices++];
      vnew->x = a*vertices[va].x + (1.-a)*vertices[vb].x;
      vnew->y = a*vertices[va].y + (1.-a)*vertices[vb].y;
      vnew->z = a*vertices[va].z + (1.-a)*vertices[vb].z;
      vnew->AllocateFaces(6);
      vnew->AllocateVertices(6);
      //for(int i = 0 ; i < 6 ; i++) vnew->e[i]=0;
      vnew->v[vnew->vnum] = va;
      vnew->e[vnew->vnum++] = -1;
      vnew->marked=1;
      for (int p = 0 ; p < vertices[va].vnum ; p++)
        if (vertices[va].v[p] == vb) {
          vertices[va].e[p] = v2;
          break;
        };
    }
    //update link v1<-->v2
    vertices[v1].v[vertices[v1].vnum]=v2;
    vertices[v1].e[vertices[v1].vnum++]=-1;
    vertices[v2].v[vertices[v2].vnum]=v1;
    vertices[v2].e[vertices[v2].vnum++]=-1;

    //modify face v0,v1,v2 (this face is the current one!!!)
    int novel_face = fn;

    fnew = &faces[fn];
    fnew->v[0] = vn0;
    for (int p = 0 ; p < vertices[vn0].fnum ; p++)
      if (vertices[vn0].f[p] == fn) {
        vertices[vn0].f[p] = novel_face;
        vertices[vn0].n[p] = 0;
      };
    fnew->v[1] = v1;
    vertices[v1].f[vertices[v1].fnum] = novel_face;
    vertices[v1].n[vertices[v1].fnum++] = 1;
    fnew->v[2] = v2;
    vertices[v2].f[vertices[v2].fnum] = novel_face;
    vertices[v2].n[vertices[v2].fnum++] = 2;

    //v3 = 2/3*vn1+1/3*vn0
    a=2./3.;
    va=vn1;
    vb = vn0;
    v3 = -1;
    for (int p = 0 ; p < vertices[va].vnum ; p++)
      if (vertices[va].v[p] == vb) {
        v3=vertices[va].e[p];
        break;
      };
    if (v3 == -1) {
      v3=nvertices;
      vnew=&vertices[nvertices++];
      vnew->x = a*vertices[va].x + (1.-a)*vertices[vb].x;
      vnew->y = a*vertices[va].y + (1.-a)*vertices[vb].y;
      vnew->z = a*vertices[va].z + (1.-a)*vertices[vb].z;
      vnew->AllocateFaces(6);
      vnew->AllocateVertices(6);
      vnew->v[vnew->vnum] = va;
      vnew->e[vnew->vnum++] = -1;
      vnew->marked=1;
      for (int p = 0 ; p < vertices[va].vnum ; p++)
        if (vertices[va].v[p] == vb) {
          vertices[va].e[p] = v3;
          break;
        };
    }

    //v4 = 2/3*vn2+1/3*vn0
    a=2./3.;
    va=vn2;
    vb = vn0;
    v4 = -1;
    for (int p = 0 ; p < vertices[va].vnum ; p++)
      if (vertices[va].v[p] == vb) {
        v4=vertices[va].e[p];
        break;
      };
    if (v4 == -1) {
      v4=nvertices;
      vnew=&vertices[nvertices++];
      vnew->x = a*vertices[va].x + (1.-a)*vertices[vb].x;
      vnew->y = a*vertices[va].y + (1.-a)*vertices[vb].y;
      vnew->z = a*vertices[va].z + (1.-a)*vertices[vb].z;
      vnew->AllocateFaces(6);
      vnew->AllocateVertices(6);
      vnew->v[vnew->vnum] = va;
      vnew->e[vnew->vnum++] = -1;
      vnew->marked=1;
      for (int p = 0 ; p < vertices[va].vnum ; p++)
        if (vertices[va].v[p] == vb) {
          vertices[va].e[p] = v4;
          break;
        };
    }

    //update link v4<-->v3
    vertices[v3].v[vertices[v3].vnum]=v4;
    vertices[v3].e[vertices[v3].vnum++]=-1;
    vertices[v4].v[vertices[v4].vnum]=v3;
    vertices[v4].e[vertices[v4].vnum++]=-1;
    //update link v3 --> vn2
    vertices[v3].v[vertices[v3].vnum]=vn2;
    vertices[v3].e[vertices[v3].vnum++]=-1;

    //create new faces!!!
    //new face vn1,vn2,v3
    novel_face = nfaces;

    fnew = &faces[nfaces++];
    fnew->v[0] = vn1;
    for (int p = 0 ; p < vertices[vn1].fnum ; p++)
      if (vertices[vn1].f[p] == fn) {
        vertices[vn1].f[p] = novel_face;
        vertices[vn1].n[p] = 0;
      };
    fnew->v[1] = vn2;
    for (int p = 0 ; p < vertices[vn2].fnum ; p++)
      if (vertices[vn2].f[p] == fn) {
        vertices[vn2].f[p] = novel_face;
        vertices[vn2].n[p] = 1;
      };
    fnew->v[2] = v3;
    vertices[v3].f[vertices[v3].fnum] = novel_face;
    vertices[v3].n[vertices[v3].fnum++] = 2;
    //new face vn2,v4,v3
    novel_face = nfaces;

    fnew = &faces[nfaces++];
    fnew->v[0] = vn2;
    //add an edge and one face to vn2
    vnew = &vertices[vn2];
    vnew->ExpandFaces(1);
    vnew->ExpandVertices(1);
    vnew->f[vnew->fnum] = novel_face;
    vnew->n[vnew->fnum++] = 0;
    vnew->v[vnew->vnum] = v3;
    vnew->e[vnew->vnum++] = -1;

    fnew->v[1] = v4;
    vertices[v4].f[vertices[v4].fnum] = novel_face;
    vertices[v4].n[vertices[v4].fnum++] = 1;

    fnew->v[2] = v3;
    vertices[v3].f[vertices[v3].fnum] = novel_face;
    vertices[v3].n[vertices[v3].fnum++] = 2;

    if (n==0) {
      first_v1=v2; // adding first v2 for orientation purposes
      first_v2=v1;
      first_v3=v3;
      first_v4=v4;
    }
  }

  ASSERT(nvertices <= maxvertices && nfaces <= maxfaces);

  ASSERT(nvertex_list==loop.npoints*3);
  // now update the vertices
  for (int n = 0 ; n < nvertex_list ; n++) {
    Vertex *v=&vertices[vertex_list[n]];
    for (int p = 0 ; p < v->vnum ; p++) {
      if (v->e[p] >= 0) {
        v->v[p] = v->e[p];
        v->e[p] = -1; //for the next triangle
      }
    }
  }
  if (vertex_list) delete [] vertex_list;

  // finding the two rings of vertices corresponding to the cut
  //first loop
  Loop ring_1(loop.npoints);
  ASSERT(first_v1 != -1 && first_v2 != -1);
  ring_1.AddPoint(first_v1);
  vertices[first_v1].marked=2;
  ring_1.AddPoint(first_v2);
  vertices[first_v2].marked=2;
  int current_vertex = first_v2;
  bool found=true;
  while (found) {
    found=false;
    Vertex *v=&vertices[current_vertex];
    for (int p = 0 ; p < v->vnum ; p++) {
      if (vertices[v->v[p]].marked==1) {
        found = true;
        current_vertex = v->v[p];
        ring_1.AddPoint(v->v[p]);
        vertices[v->v[p]].marked=2;
        break;
      }
    }
  }
  ASSERT(ring_1.npoints == loop.npoints);

  //second loop
  Loop ring_2(loop.npoints);
  ASSERT(first_v3 != -1 && first_v4 != -1);
  ring_2.npoints=0;
  ring_2.AddPoint(first_v3);
  vertices[first_v3].marked=2;
  ring_2.AddPoint(first_v4);
  vertices[first_v4].marked=2;
  current_vertex = first_v4;
  found=true;
  while (found) {
    found=false;
    Vertex *v=&vertices[current_vertex];
    for (int p = 0 ; p < v->vnum ; p++) {
      if (vertices[v->v[p]].marked==1) {
        found = true;
        current_vertex = v->v[p];
        ring_2.AddPoint(v->v[p]);
        vertices[v->v[p]].marked=2;
        break;
      }
    }
  }
  ASSERT(ring_2.npoints == loop.npoints);

  KnitPatch(ring_1,pdisk);
  KnitPatch(ring_2,pdisk);

  ASSERT(nvertices <= maxvertices && nfaces <= maxfaces);

  ExpandMarks(2,4);
  for (int n = 0 ; n < nvertices ; n++)
    if (vertices[n].vnum!=vertices[n].fnum)
      vertices[n].marked=0;//border vertex

  // finding vertices that are marked with 4
  //int * vert = new int[2*pdisk->disk.nvertices],nvert=0;
  int nvert=0;
  for (int n = 0 ; n < nvertices ; n++)
    if (vertices[n].marked==4)
      nvert++;
  int * vert = new int[nvert];
  nvert=0;
  for (int n = 0 ; n < nvertices ; n++)
    if (vertices[n].marked==4)
      vert[nvert++]=n;
  Smooth(150,vert,nvert);
  if (vert) delete [] vert;

  InitSurface();
  GetEuler();

  ASSERT(IsSurfaceValid(0));
}

bool Surface::IsSurfaceValid(int verbose) {
  int nbord=0,npbm=0,nfv1=0,nfv2=0, nadd=0,nmiss=0;


  for (int n = 0 ; n < nvertices ; n++)
    vertices[n].marked=0;
  for (int n = 0 ; n < nfaces ; n++)
    faces[n].marked=0;

  if (verbose) cout << endl;

  for (int n = 0 ; n < nvertices ; n++) {
    Vertex *v=&vertices[n];
    for (int p = 0 ; p < v->vnum ; p++) {
      int vp = v->v[p];
      int fn=_FindFace(n,vp,-1);
      if (fn==-1) {
        npbm++;
        v->marked=2;
        break;
      }
      if (_FindFace(n,vp,fn)==-1) { //border
        faces[fn].marked=1;
        nbord++;
        v->marked=1;
        if (v->vnum!=v->fnum+1)
          nfv2++;
      }
    }
  }
  for (int n = 0 ; n < nvertices ; n++) {
    Vertex *v=&vertices[n];
    if (v->marked) continue;
    if (v->fnum!=v->vnum) nfv1++;
  }
  for (int n = 0 ; n < nfaces ; n++) {
    bool test1 = false;
    if (faces[n].marked) {
      test1=true;
    };
    bool test2=false;
    for (int i = 0 ; i < 3 ; i++)
      if (faces[n].f[i]==-1) test2=true;
    if (test1 && !test2) {
      nadd++;
    }
    if (!test1 && test2) {
      nmiss++;
    }
  }
  if (verbose) {
    cout <<endl << " we have [" << npbm << " , " << nbord/2 << " , " << nadd << " , " << nmiss
    << " , " << nfv1 << " , " << nfv2 << " ]   " << endl ;
    if (nfv2 || nfv1) cout << "pbm";
  }
  if (npbm || nmiss || nadd || nfv1 || nfv2)
    return false;

  return true;
}

void Surface::IncreaseEuler(int nloops,int maxinitface) {
  int nattempts = 10;


  if (maxinitface == -1) maxinitface = nfaces;

  if (euler==1) return;

  FastLoop FL(*this);
  int *defectfaces = new int[nfaces];
  for (int n = 0 ; n < nfaces ; n++)
    defectfaces[n]=n;
  FL.SetDefectList(nfaces,defectfaces);

  Surface *best_patch = 0 ;

  for (int n = 0 ; n < nloops ; n++) {

    int nattemptstofinaloop=10;

    // the seed face is randomly drawn
    Loop *loop;
    do {
      nattemptstofinaloop--;
      if (nattemptstofinaloop==0) {
        if (defectfaces) delete [] defectfaces;
        ErrorExit("could not find a loop in IncreaseEuler");
      }
      int seed_face = Random(maxinitface);
      loop = FL.FindLoop(seed_face);
#if PRINT_MODE
      if (loop==0) cout << "l";
#endif
    } while (loop==0);

    Loop *cutting_loop;
    if (loop[1].npoints) {
      if (Random(2)==1) cutting_loop = &loop[1];
      else cutting_loop = &loop[0];
    } else cutting_loop = &loop[0];

    //if(cutting_loop==&loop[0]) cout << "'0'";

    Surface *cut_patch = Clone();
    cut_patch->CutLoop(*cutting_loop);
    if (loop) delete [] loop;

    //check euler
    cut_patch->GetEuler();
    if (cut_patch->euler != euler + 2) {
      delete cut_patch;
      if (nattempts-- > 0) {
        n--;
#if PRINT_MODE
        cout << "a";
#endif
      }
      continue;
    }

    if (best_patch == 0)
      best_patch = cut_patch;
    else {
      //evaluate if this patch is better than the previous best one
      delete cut_patch;
    }
  }
  if (defectfaces) delete [] defectfaces;
  if (best_patch == 0) return;

  // transfer the faces and vertices
  if (vertices) delete [] vertices;
  vertices = best_patch->vertices;
  nvertices = best_patch->nvertices;
  maxvertices = best_patch->maxvertices;
  best_patch->vertices=0;
  if (faces) delete [] faces;
  faces = best_patch->faces;
  nfaces = best_patch->nfaces;
  maxfaces = best_patch->maxfaces;
  best_patch->faces = 0;
  nedges = best_patch->nedges;
  euler = best_patch->euler;
  if (best_patch) delete best_patch;
}

void Surface::CorrectTopology() {
  int niters = (1-GetEuler())/2;

  int maxinitface = nfaces;
  for (int n = 0 ; n < niters ; n++) {
    IncreaseEuler(1,maxinitface);
    GetEuler();
  }
}

//resultat  = int s:
// >0 if success with initial seed s,
// -1 if failure in finding the loop
// -2 if failure in cutting the patch
// a negative seed specifies the cutting mode
int Surface::CutPatch(int seed, int maxinitface, int nattempts, int very_small_patch) {
  int mode = 0 ;
  if (seed < 0) {
    if (seed == -2) {
      mode=1; // minimal loop mode
    }
    seed=-1;
  }

  // assume that the euler number was computed for this patch
  if (euler==1) return 0;
  int init_euler = euler;

  if (maxinitface < 0 || maxinitface > nfaces) maxinitface = nfaces; //the max value of the seed face
  if (seed < 0 || seed >= maxinitface) seed = -1;
  int seed_face = seed; // the seed face
  if (nattempts < 1) nattempts = 1;

  FastLoop FL(*this);
  int *defectfaces = new int[nfaces];
  for (int n = 0 ; n < nfaces ; n++)
    defectfaces[n]=n;
  FL.SetDefectList(nfaces,defectfaces);

  // finding the loop in at most nattempts attempts
  Loop *loop=0;
  int ntries = 0;
  do {
    ntries++;
    if (ntries > nattempts) {
      if (defectfaces) delete [] defectfaces;
      return -1;
    }
    //find a seed face
    do {
      seed_face = Random(maxinitface);
    } while (seed_face < 0 || seed_face >= maxinitface);
    //check if we use the initial seed face
    if (seed >= 0 && seed < maxinitface) {
      seed_face = seed;
      seed = -2;
    };
    if (mode) {
      if (mode==1) {
        loop = new Loop[2];
        FL.FindMinimalLoop(loop[0], maxinitface , nattempts);
        if (loop->npoints==0) { // did not work ! try again!
          delete [] loop;
          loop=0;
        } else
          nattempts=1; //success -> exit
      }
    } else
      loop = FL.FindLoop(seed_face);
  } while (loop==0);
  if (defectfaces) delete [] defectfaces;

  //picking a cutting loop
  Loop *cutting_loop;
  if (loop[1].npoints) {
    if (Random(2)==1) cutting_loop = &loop[1];
    else cutting_loop = &loop[0];
  } else cutting_loop = &loop[0];

  // cutting the patch
  CutLoop(*cutting_loop,very_small_patch);
  if (loop) delete [] loop;
  //computing the new euler number
  GetEuler();
  //checking the validity of the new cut
  if (IsSurfaceValid() == false || euler != init_euler+2) return -2; //failure

  return (seed_face); //success
}
