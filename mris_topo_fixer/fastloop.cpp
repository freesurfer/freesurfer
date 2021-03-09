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


#include "fastloop.h"
#include "segment.h"

FastLoop::FastLoop(Surface &s):surface(&s) {
  defect_faces=0;
  ndefect_faces=0;

  facedata = new FaceData[surface->nfaces+1]; //one extra face for the border
  vertexdata = new VertexData[surface->nvertices];
  _InitFaceData();

  FM_trial_heap = new FaceHeap(HeapCompare(facedata));
  FCC_trial_heap = new FaceHeap(HeapCompare(facedata));
  Init();
}

FastLoop::~FastLoop(void) {
  delete [] facedata;
  delete [] vertexdata;
  delete FM_trial_heap;
  delete FCC_trial_heap;
}

void FastLoop::_InitFaceData(void) {
  int nfaces = surface->nfaces;

  for (int n = 0 ; n < nfaces ; n++) {
    FaceData *fdst = &facedata[n];
    Face *fsrc = &surface->faces[n];
    fdst->x=fsrc->x;
    fdst->y=fsrc->y;
    fdst->z=fsrc->z;
    for (int i = 0 ; i < 3 ; i++) {
      fdst->v[i]=fsrc->v[i];
      int fn = fsrc->f[i];
      if (fn == -1) fn=nfaces; //the border faces point towards the last face
      fdst->f[i]=fn;
    }
  }
  //update the extra face to avoid abnormal behaviors
  facedata[nfaces].border=1;
  facedata[nfaces].fmState=eForbidden;
  facedata[nfaces].fccState=eForbidden;
  facedata[nfaces].fmFather=-1;
  facedata[nfaces].fccFather=-2;
  facedata[nfaces].nfather=-3;
  facedata[nfaces].fccLabel=-4;
  facedata[nfaces].f[0]=-5;
  facedata[nfaces].f[1]=-5;
  facedata[nfaces].f[2]=-5;
  facedata[nfaces].v[0]=-6;
  facedata[nfaces].v[1]=-6;
  facedata[nfaces].v[2]=-6;
}


void FastLoop::Init(void) {
  while (!FM_trial_heap->empty()) FM_trial_heap->pop();
  while (!FCC_trial_heap->empty()) FCC_trial_heap->pop();
}

void FastLoop::SetSeed(int seed) {
  _AddAliveFace(seed);
  facedata[seed].val=0.0;
  facedata[seed].fmFather=-1;
  for (int i = 0 ; i < 3 ; i++) {
    _UpdateFace(facedata[seed].f[i],seed);
  }
}

int FastLoop::_AddAliveFace(int fno) {
  facedata[fno].fmState = eAlive;
  facedata[fno].fccState = eForbidden;

  for (int i = 0 ; i < 3 ; i++) {
    if (vertexdata[facedata[fno].v[i]].fmState==eFar) {
      vertexdata[facedata[fno].v[i]].fmState=eAlive;
      vertexdata[facedata[fno].v[i]].fmFather=fno;
    }
  }

  return 0;
}

void FastLoop::_UpdateFace(int fdst, int fsrc) {
  const eState st = (eState)facedata[fdst].fmState;

  if (st == eFar) {
    facedata[fdst].fmState = eTrial;
    facedata[fdst].fmFather=fsrc;
    facedata[fdst].nfather=1;
    facedata[fdst].val = facedata[fsrc].val + _Distance(fdst,fsrc);
    FM_trial_heap->push(fdst);
  } else if (st == eTrial) { // already in the list
    facedata[fdst].nfather++;
  }
}

double FastLoop::_Distance(int fdst, int fsrc) {
  return __norm(facedata[fdst].x-facedata[fsrc].x,facedata[fdst].y-facedata[fsrc].y,facedata[fdst].z-facedata[fsrc].z);
}


int FastLoop::_AddTrialFace(int fno) {
  for (int i = 0 ; i < 3 ; i++) {
    if (vertexdata[facedata[fno].v[i]].fmState==eFar)
      return -1;
  }

  if (facedata[fno].nfather>=2) return -1;

  //now evaluate if we have a real problem!!!
  //check if there is a neighbor of fno with nfather>=2
  for (int i = 0 ; i < 3 ; i++) {
    int fn=facedata[fno].f[i];
    if (facedata[fn].fmState==eTrial && facedata[fn].nfather>=2) {
      // add this face to the list of alive faces...
      _AddAliveFace(fn);
      return -1;
    }
  }
  //find problematic face
  int face=-1;
  for (int i = 0 ; i < 3 ; i++) {
    int test=0;
    for (int j = 0 ; !test && j < 3 ; j++) {
      if (facedata[facedata[fno].fmFather].v[j]==facedata[fno].v[i]) {
        test=1;
        break;
      }
    }
    if (!test) {
      face = vertexdata[facedata[fno].v[i]].fmFather;
      return face;
    }
  }
  ASSERT(face!=-1);
  return 0;
}



int FastLoop::_Run(int &stopping_face) {
  while (!FM_trial_heap->empty()) {

    // we remove the face at the top of the stack and add it to the list of alive faces
    int fno = FM_trial_heap->top();


    if (facedata[fno].fmState==eAlive) {//this face was already added!
      FM_trial_heap->pop();
      continue;
    }



    //check if this face (fno) can be added
    int res=_AddTrialFace(fno);
    if (res >= 0) {
      stopping_face = res;
      return fno;
    }

    FM_trial_heap->pop();
    _AddAliveFace(fno);

    // we look at its neighboring faces and we update them
    for (int i = 0 ; i < 3 ; i++)
      _UpdateFace(facedata[fno].f[i],fno);
  }
  return -1;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


void FastLoop::_UpdateSegmentFaces(int which_segment) {
  Segment *segment = &segments[which_segment];

  int segment_add=which_segment,segment_keep=(which_segment+1)%2;

  int euler = surface->GetEuler(segment->GetPointList(),segment->size());
  if (euler < 0 ) {
    segment_add=(which_segment+1)%2;
    segment_keep=which_segment;
  }
#if PRINT_ERROR_MODE
  cout << "segment " << which_segment << " has an euler number of " << euler << endl;

  cout << " we have " << segment_add << ": "
  << surface->GetEuler(segments[segment_add].GetPointList(),segments[segment_add].size())
  << " and " << segment_keep << ": "
  << surface->GetEuler(segments[segment_keep].GetPointList(),segments[segment_keep].size())
  << endl;
#endif

  const int* face_list;
  int list_size;

  //add the faces of segments[segment_add] to the list of alive faces
  if (nsegments[segment_add]==0) {
    face_list = segments[segment_add].GetPointList();
    list_size = segments[segment_add].size();
    for (int n = 0 ; n < list_size ; n++)
      _AddAliveFace(face_list[n]);
    for (int n = 0 ; n < list_size ; n++)
      for (int i = 0 ; i < 3 ; i++)
        _UpdateFace(facedata[face_list[n]].f[i],face_list[n]);
  } else {
#if PRINT_ERROR_MODE
    cout << "m";
#endif
    face_list = segments[segment_add].GetPointList();
    list_size = segments[segment_add].size();
    for (int n = 0 ; n < list_size ; n++) {
      FaceData *face = &facedata[face_list[n]];
      face->fccFather = -1;
      face->val=face->val2;
      face->fccState=eFar;
      face->fccLabel=-1;
    }
  }

  //keep the faces of segments[segment_keep]
  face_list = segments[segment_keep].GetPointList();
  list_size = segments[segment_keep].size();
  for (int n = 0 ; n < list_size ; n++) {
    FaceData *face = &facedata[face_list[n]];
    face->fccFather = -1;
    face->val=face->val2;
    face->fccState=eFar;
    face->fccLabel=-1;
  }
  while (!FCC_trial_heap->empty()) {
    int fn=FCC_trial_heap->top();
    FCC_trial_heap->pop();
    FaceData *face = &facedata[fn];
    face->fccFather = -1;
    face->val=face->val2;
    face->fccState=eFar;
    face->fccLabel=-1;
  }

}

int FastLoop::_FindSeedFaces(int conflicting_face, int& init_fn1, int& init_fn2) {
  int number_of_faces=0;
  init_fn1 = init_fn2=-1;
  FaceData *face = &facedata[conflicting_face];
  for (int i = 0 ; i < 3 ; i++) {
    if (facedata[face->f[i]].fmState==eAlive) continue;
    if (facedata[face->f[i]].fmState==eForbidden) continue; //the border face
    if (number_of_faces==0) {
      init_fn1=face->f[i];
      number_of_faces++;
    } else {
      number_of_faces++;
      init_fn2=face->f[i];
    }
  }
  return number_of_faces; // number of found faces...
}

void FastLoop::_InitSegment(int which_segment, int fn) {

  nsegments[which_segment]=0;
  segments[which_segment].clear();
  segments[which_segment].SetMark(0);
  if (fn==-1) {
    segments[which_segment].SetMark(1);
    return;
  }
  segments[which_segment].AddPoint(fn);
  FaceData *face = &facedata[fn];
  for (int i = 0 ; i < 3 ; i++) {
    if (facedata[face->f[i]].border)
      segments[which_segment].SetMark(1);
  }
  face->fccState=eAlive;
  face->val2=face->val;
  face->val=0.0;
  face->fccFather = -1;
  face->fccLabel=which_segment;
  for (int i = 0 ; i < 3 ; i++)
    _UpdateSegmentFace(face->f[i],fn);

}

void FastLoop::_UpdateSegmentFace(int fdst, int fsrc) {
  int label = facedata[fsrc].fccLabel;
  FaceData *face = &facedata[fdst];

  if (face->fccState == eFar) {
    face->fccState = eTrial;
    face->fccLabel=label;
    face->fccFather = fsrc;
    face->val2 = face->val;
    face->val = facedata[fsrc].val + _Distance(fdst,fsrc);
    FCC_trial_heap->push(fdst);
    nsegments[label]++;
  }

}

int FastLoop::_InitSegmentation(int fn) {

  ASSERT(FCC_trial_heap->empty());

  int init_fn1,init_fn2;

  int number_of_faces = _FindSeedFaces(fn,init_fn1,init_fn2);

#if PRINT_ERROR_MODE
  cout << "init1 = " << init_fn1 << " init2 = " << init_fn2 << endl;
  if (number_of_faces < 2 ) cout << "ONLY "<< number_of_faces << " FACES!!!!!" <<endl;
#endif

  _InitSegment(0,init_fn1);
  _InitSegment(1,init_fn2); //potentially init_fn2 == -1

  return number_of_faces;
}

int FastLoop::_FastSegmentation() {

  while (nsegments[0] && nsegments[1]) {

    int fn = FCC_trial_heap->top();
    FCC_trial_heap->pop();
    int label = facedata[fn].fccLabel;
    nsegments[label]--;

    //test if fn is a neighbor of another label
    int fn2 = _CheckAdjacency(fn);

    if (fn2 >= 0) { // connected component
      //empty the list
      while (FCC_trial_heap->empty()) FCC_trial_heap->pop();
      final_face[label]=fn;
      final_face[(label+1)%2]=fn2;
#if PRINT_ERROR_MODE
      cout << endl << "The face "<< fn << " is Neighbor with " << fn2 << endl;
#endif
      return -1;
    }

    int which_segment=facedata[fn].fccLabel;
    segments[which_segment].AddPoint(fn);
    for (int i = 0 ; i < 3 ; i++) {
      if (facedata[facedata[fn].f[i]].border)
        segments[which_segment].SetMark(1);
    }

    facedata[fn].fccState=eAlive;
    for (int i = 0 ; i < 3 ; i++)
      _UpdateSegmentFace(facedata[fn].f[i],fn);
  }
  if (nsegments[1]) return 0;
  else return 1;
}

int FastLoop::_CheckAdjacency(int fn) {
  FaceData *face=&facedata[fn];
  int label = face->fccLabel;

  //find the smallest face
  double val=1e10;
  int which_face = -1;
  for (int i = 0 ; i < 3 ; i++) {
    FaceData *f = &facedata[face->f[i]];
    if ((f->fccState == eAlive ||f->fccState == eTrial) && f->fccLabel!=label) {
      if (f->val < val) {
        val = f->val;
        which_face=face->f[i];
      }
    }
  }

  return which_face;
}

//extract a loop starting from the two faces init_fn1 and init_fn2
//these two faces are adjacent through a vertex
// return 0 if fails to find path
int FastLoop::_ExtractFirstLoop(Loop &loop, int init_fn1, int init_fn2) {
  int final_fn1=-1,final_fn2=-1;

  //we can use fccLabel to extract the adjacent faces
  //since this variable is not necessary anymore
  for (int n = 0 ; n < ndefect_faces ; n++) {
    facedata[defect_faces[n]].fccLabel=0;
    for (int i = 0 ; i < 3 ; i++)
      vertexdata[facedata[n].v[i]].fmState=0;
  }

  // we first count the number of faces
  int number_of_faces=0,fn;

  fn=init_fn1;
  while (fn!=-1) {
    number_of_faces++;
    //mark this face and its vertices
    facedata[fn].fccLabel = 1;
    for (int i = 0 ; i < 3 ; i++)
      vertexdata[facedata[fn].v[i]].fmState = 1;
    ASSERT(fn!=facedata[fn].fmFather);
    if (fn==facedata[fn].fmFather) {
      //  surface->WriteFile("pp2.3d",1);
      //exit(-1);
	  cout << "ExtractFirstLoop: father same as son..." << endl;	
	  return 0;	
      ErrorExit("ExtractFirstLoop: father same as son...");
    }
    fn=facedata[fn].fmFather;
  }
  fn=init_fn2;
  while (fn!=-1) {
    number_of_faces++;
    //mark this face and its vertices
    facedata[fn].fccLabel=2;
    for (int i = 0 ; i < 3 ; i++)
      vertexdata[facedata[fn].v[i]].fmState=2;
    fn=facedata[fn].fmFather;
  }

#if 1
  int vno = _FindCommonVertex(init_fn1,init_fn2);
  vertexdata[vno].fmState = 0;
#else
  // just making sure we are not doing a stupid mistake
  for (int i = 0 ; i < 3 ; i++)
    vertexdata[facedata[init_fn1].v[i]].fmState=0;
  for (int i = 0 ; i < 3 ; i++)
    vertexdata[facedata[init_fn2].v[i]].fmState=0;
#endif

  loop.Alloc(number_of_faces); // we have at max number_of_faces faces (at least for now!)

  fn=init_fn1;
  while (fn!=-1) {
    //check adjacency
    for (int i = 0 ; i < 3 ; i++) {
      FaceData *face = &facedata[facedata[fn].f[i]];
      if (face->fmState == eAlive && face->fccLabel == 1) {
        int nfn = facedata[fn].f[i];
        //check if we had already added this face
        bool test=false;
        for (int n = 0 ; n < loop.npoints ;n++) {
          if (loop.points[n]==nfn)
            test=true;
        }
        //we had already added this neighboring face!!!
        if (test) {
          while (loop.points[loop.npoints-1]!=nfn)
            loop.npoints--;
          break;
        }
      }
    }
    bool adjacent_vertex = false;
    for (int i = 0 ; i < 3 ; i++) {
      VertexData *v = &vertexdata[facedata[fn].v[i]];
      if (v->fmState==2) {
        adjacent_vertex=true;
        v->fmState=3;
      }
    }
    //we can add this face to the list
    loop.AddPoint(fn);
    final_fn1=fn;
    fn=facedata[fn].fmFather;
    if (adjacent_vertex) break;
  }

  fn=init_fn2;
  while (fn!=-1) {
    for (int i = 0 ; i < 3 ; i++) {
      FaceData *face = &facedata[facedata[fn].f[i]];
      if (face->fmState == eAlive && face->fccLabel == 2) {
        int nfn = facedata[fn].f[i];
        //check if we had already added this face
        bool test=false;
        for (int n = 0 ; n < loop.npoints ;n++) {
          if (loop.points[n]==nfn)
            test=true;
        }
        //we had already added this neighboring face!!!
        if (test) {
          while (loop.points[loop.npoints-1]!=nfn)
            loop.npoints--;
          break;
        }
      }
    }
    bool adjacent_vertex = false;
    for (int i = 0 ; i < 3 ; i++) {
      VertexData *v = &vertexdata[facedata[fn].v[i]];
      if (v->fmState==3) {
        adjacent_vertex=true;
        v->fmState=1;
      }
    }
    //we can add this face
    loop.AddPoint(fn);
    final_fn2=fn;
    fn=facedata[fn].fmFather;
    if (adjacent_vertex) break;
  }

#if PRINT_ERROR_MODE
  cout << "first: " ;
  loop.Print();
#endif

  //finally find path in between init_fn1 and init_fn2
  if(!_FindFacePath(loop,init_fn1,init_fn2)){
	return 0;		
  }
#if PRINT_ERROR_MODE
  cout << "second: ";
  loop.Print();
#endif

  // and final_fn1 and final_fn2
  if(!_FindFacePath(loop,final_fn1,final_fn2)){
	return 0;
  }
#if PRINT_ERROR_MODE
  cout << "third: ";
  loop.Print();
#endif

  return 1;
}

// find the only vertex that is common to init_fn1 and init_fn2
int FastLoop::_FindCommonVertex(int init_fn1,int init_fn2) {
  for (int i = 0 ; i < 3 ; i++) {
    int vn = facedata[init_fn1].v[i];
    for (int j = 0 ; j < 3 ; j ++)
      if (vn == facedata[init_fn2].v[j]) return vn;
  }
  return -1;
}

//find the single next face with vno and that is not marked (through fccLabel)
int FastLoop::_FindNextFace(int next_fn,int vno) {
  ASSERT(next_fn < surface->nfaces);

  if (next_fn == surface->nfaces) return -2;

  for (int i = 0 ; i < 3 ; i++) {
    int fn = facedata[next_fn].f[i];
    //  if(facedata[fn].fccLabel) continue;
    if (facedata[fn].fccLabel==3) return -1;
    if (facedata[fn].fccLabel) continue;
    for (int j = 0 ; j < 3 ; j++) {
      if (facedata[fn].v[j]==vno) return fn;
    }
  }
  return -2; //could not find a face!!
}

// find a path of faces in between two faces
// that have one single common vertex
// return 0 if fails to find path
int FastLoop::_FindFacePath(Loop &loop, int init_fn1,int init_fn2) {
  /////////////////////////
  // Initialization

  //we can use fccLabel to extract the path of adjacent faces
  //since this variable is not necessary anymore
  for (int n = 0 ; n < ndefect_faces ; n++) {
    facedata[defect_faces[n]].fccLabel=0;
  }
  //mark the current loop faces + neighbors!!
  for (int n = 0 ; n < loop.npoints ; n++) {
    facedata[loop.points[n]].fccLabel=1;
#if 0
    continue; //test
    if (loop.points[n]==init_fn1) continue;
    if (loop.points[n]==init_fn2) continue;
    for (int i = 0 ; i < 3 ; i++) {
      if (facedata[loop.points[n]].f[i]>=0)
        facedata[facedata[loop.points[n]].f[i]].fccLabel=1;
    }
#endif
  }
  facedata[init_fn1].fccLabel=2;
  facedata[init_fn2].fccLabel=3;
  //avoid the marked faces
  //find the neighboring faces of init_fn1
  int fn1=-1,fn2=-1;
  for (int i = 0 ; i < 3 ; i++) {
    int n = facedata[init_fn1].f[i];
    if (facedata[n].fccLabel == 0) {
      if (fn1==-1) fn1 = n;
      else fn2 = n;
    }
  }
  //find the common vertex
  int vno = _FindCommonVertex(init_fn1,init_fn2);
  if (vno==-1){
   cout << "_FindFacePath: could not find the common vertex! " << endl;
   return 0;
   ErrorExit("_FindFacePath: could not find the common vertex! ");
  }
  int fnum = surface->vertices[vno].fnum;

  //////////////////////////////////
  // Now we are ready to find and evaluate the two paths!
  bool found;
  int fn , next_fn;
  double path[2];
  int *PathFaces[2],nfaces[2];

  // first path
  path[0] = 0;
  nfaces[0] = 0 ;
  next_fn = init_fn1;
  PathFaces[0] = new int[fnum]; //at max fnum faces !!!
  if (fn2 >= 0) facedata[fn2].fccLabel=1; // we want to avoid this face
  found = false;
  while (!found) {
    fn=next_fn;
    next_fn = _FindNextFace(next_fn,vno);
    if (next_fn==-1)
      found = true;
    else if (next_fn==-2) {
      found=true;
      nfaces[0]=0;
    } else {
      path[0] += _Distance(fn,next_fn);
      facedata[next_fn].fccLabel=4;
      PathFaces[0][nfaces[0]++]=next_fn;
    }
  }

  //second path if ( fn2 >= 0 )
  path[1] = 0;
  nfaces[1] = 0;
  next_fn=init_fn1;
  PathFaces[1] = new int[fnum];
  if (fn2>=0) {
    facedata[fn2].fccLabel=0; // now we want to find this face
    found = false;
    while (!found) {
      fn=next_fn;
      next_fn = _FindNextFace(next_fn,vno);
      if (next_fn==-1)
        found = true;
      else if (next_fn==-2) {
        found=true;
        nfaces[1]=0;
      } else {
        path[1] += _Distance(fn,next_fn);
        facedata[next_fn].fccLabel=4;
        PathFaces[1][nfaces[1]++]=next_fn;
      }
    }
  }

  //////////////////////////////////////////////
  // Finally, pick the shortest path
  int which_path=0;
  if (nfaces[0]*nfaces[1]==0) {
    if (!nfaces[0] && !nfaces[1]) { // should never happen!
      delete [] PathFaces[0];
      delete [] PathFaces[1];
      //surface->WriteFile("./test.3d",1);
	  cout << "_FindFacePath: could not find path!" << endl;
	  return 0;
      ErrorExit("_FindFacePath: could not find path!");
    };
    if (!nfaces[0]) which_path = 1;
  } else {
    if (path[1] < path[0]) which_path = 1;
  }

  //add the faces of which_path to the loop
  for (int n = 0 ; n < nfaces[which_path] ; n++)
    loop.AddPoint(PathFaces[which_path][n]);

  delete [] PathFaces[0];
  delete [] PathFaces[1];
  return 1;
}

//extract a loop starting from the two adjacent faces final_face[i] (i=0 & i=1)
int FastLoop::_ExtractSecondLoop(Loop& loop, int init_fn) {
  int init_fn1=final_face[0];
  int init_fn2=final_face[1];

  // as usual, we start by counting the number of faces
  int number_of_faces=0,fn;

  fn=init_fn1;
  while (fn!=-1) {
    number_of_faces++;
    fn=facedata[fn].fccFather;
  }
  fn=init_fn2;
  while (fn!=-1) {
    number_of_faces++;
    fn=facedata[fn].fccFather;
  }
  // allocate the loop
  loop.Alloc(number_of_faces+1);//add an extra face for the seed face

  fn=init_fn1;
  while (fn!=-1) {
    loop.AddPoint(fn);
    fn=facedata[fn].fccFather;
  }
  fn=init_fn2;
  while (fn!=-1) {
    loop.AddPoint(fn);
    fn=facedata[fn].fccFather;
  }
  loop.AddPoint(init_fn);

  return 0;
}

int FastLoop::_OrderLoop(Loop &loop) {
  for (int n = 0 ; n < ndefect_faces ; n++) {
    facedata[defect_faces[n]].fccState=0;
  }
  for (int n = 0 ; n < loop.npoints ; n++)
    facedata[loop.points[n]].fccState = eTemporary;

  int *newFaces = new int[loop.npoints];
  int npoints=0;
  bool found=true;
  int current=loop.points[0];
  while (found) {
    newFaces[npoints++]=current;
    facedata[current].fccState=0;
    found = false;
    //find the next face
    int next = -1;
    for (int i = 0 ; i < 3 ; i++) {
      if (facedata[facedata[current].f[i]].fccState == eTemporary) {
        next = facedata[current].f[i];
        break;
      }
    }
    if (next >= 0) {
      found=true;
      current = next;
    }
  }
  ASSERT(npoints == loop.npoints);

  delete [] loop.points;
  loop.points = newFaces;
  loop.maxpoints=loop.npoints;

  return loop.npoints;
}

int FastLoop::_SimplifyLoop(Loop &loop) {
  int npoints=0;

  //using fccLabel
  for (int n = 0 ; n < ndefect_faces ; n++) {
    facedata[defect_faces[n]].fccState=0;   //to mark the loop faces (with eTemporary)
    facedata[defect_faces[n]].fmState=0;   //to mark the triple faces (with 1)
    facedata[defect_faces[n]].fmFather=0; //to mark the processed faces
  }

  for (int n = 0 ; n < loop.npoints ; n++)
    facedata[loop.points[n]].fccState = eTemporary;  //to mark the loop faces

  for (int n = 0 ; n < loop.npoints ; n++) {
    int count=0;
    for (int i = 0 ; i < 3 ; i++)
      if (facedata[facedata[loop.points[n]].f[i]].fccState==eTemporary)
        count++;
    if (count==3) {
      npoints++;
      facedata[loop.points[n]].fmState=1;
    }
  }
  if (npoints==0)
    return _OrderLoop(loop);
  ; //no triple face to remove

  //pick the first non-triple face
  int current_face=-1;
  for (int n = 0 ; n < loop.npoints ; n++)
    if (facedata[loop.points[n]].fmState==0) {
      current_face=loop.points[n];
      break;
    };
  ASSERT(current_face!=-1);

  bool found = true;
  while (found) {
    found = false;
    //mark the current face as processed
    facedata[current_face].fmFather = 1;

    //find a neighbor of current_face that is not yet processed
    int next_face = -1;
    for (int i = 0 ; i < 3 ; i++) {
      if (facedata[facedata[current_face].f[i]].fccState != eTemporary) continue; //not a loop face
      if (facedata[facedata[current_face].f[i]].fmFather==0) { //not yet processed
        next_face = facedata[current_face].f[i];
        break;
      }
    }
    if (next_face==-1) break; //no more faces to process!

    // we have found another face
    found = true;

    if (facedata[next_face].fmState == 0 ) { //this face is non-triple -> continue
      current_face = next_face;
      continue;
    }

    // we have found a triple face - we need to carefully analyze the remaining faces

    //we already process this face
    facedata[next_face].fmFather=1;
    current_face=next_face;

    //the current face has two neighboring faces (since it is a triple face)
    int fn1=-1,fn2=-1;
    for (int i = 0 ; i < 3 ; i++) {
      if (facedata[facedata[current_face].f[i]].fmFather == 1 ) continue; // already processed
      if (fn1==-1)
        fn1 = facedata[current_face].f[i];
      else
        fn2 = facedata[current_face].f[i];
    }
    if (facedata[fn1].fmState == 0) { //swap faces if fn1 is non triple
      int tmp=fn1;
      fn1=fn2;
      fn2=tmp;
    }
    // we temporary mark these two faces
    facedata[fn1].fmFather = 2 ;
    facedata[fn2].fmFather = 2 ;

    if (facedata[fn2].fmState == 0) { //only fn1 is a triple face!
      //we need to set all the faces that are neighbors of fn2 to 3
      Loop flist;
      flist.AddPoint(fn2);
      while (flist.npoints) {
        //take the last point
        int fn = flist.End();
        facedata[fn].fmFather = 2;
        flist.Pop();
        for (int i = 0 ; i < 3 ; i++) {
          if (facedata[facedata[fn].f[i]].fccState != eTemporary ) continue; //not a loop face
          if (facedata[facedata[fn].f[i]].fmFather==0) { //not yet processed
            flist.AddPoint(facedata[fn].f[i]);
          }
        }
      }
      current_face = fn1;
      continue;
    }

    // In this worst case situation, we have to decide which face (of fn1 and fn2) is the right one
    // to do so we find the connected components of the fccState==eTemporary && fmFather==0
    // we arbitrarily start with the face fn1
    // we need to find the two neighbors of fn1 (since fn1 is a triple face)
    int fn11=-1,fn12=-1;
    for (int i = 0 ; i < 3 ; i++) {
      if ( facedata[facedata[fn1].f[i]].fmFather ) continue; //already processed
      if (fn11==-1)
        fn11 = facedata[fn1].f[i];
      else
        fn12 = facedata[fn1].f[i];
    }
    //we start arbitrarily with fn11
    Loop flist;
    flist.AddPoint(fn11);
    int ending_face=-1;
    bool found = false;
    //se set fn1 to 3 (temporary value)
    facedata[fn1].fmFather = 3 ;
    while (!found && flist.npoints) {
      int fn = flist.End();
      flist.Pop();
      facedata[fn].fmFather = 3;
      //look at the neighbors
      for (int i = 0 ; i < 3 ; i++) {
        int nfn = facedata[fn].f[i];
        if (facedata[nfn].fccState != eTemporary) continue; //not part of the loop
        // now check what is the status of this point
        if (facedata[nfn].fmFather == 0) { //not yet processed
          flist.AddPoint(nfn);
          continue;
        }
        //we know that this point has already been processed!
        if (facedata[nfn].fmFather == 1) {
          found=true;
          ending_face = nfn;
          break;
        }
      }
    }
    if (ending_face==-1) { // two cases to consider
      if (facedata[fn12].fmFather == 0) { //fn1 and fn12 are the right faces
        facedata[fn1].fmFather = 1;
        for (int n = 0 ; n < loop.npoints ; n++)
          if (facedata[loop.points[n]].fmFather==3)
            facedata[loop.points[n]].fmFather=2;
        current_face = fn12;
        continue;
      } else { //fn2 is the right face
        int fn22 = -1;
        for (int i = 0 ; i < 3 ; i++) {
          if (facedata[facedata[fn2].f[i]].fmFather==0) { // the only remaining face
            fn22=facedata[fn2].f[i];
            break;
          }
        }
        facedata[fn2].fmFather=1;
        for (int n = 0 ; n < loop.npoints ; n++)
          if (facedata[loop.points[n]].fmFather==3)
            facedata[loop.points[n]].fmFather=2;
        current_face = fn22;
        continue;
      }
    } else { //fn1 and fn11 are the right faces
      facedata[fn1].fmFather=1;
      for (int n = 0 ; n < loop.npoints ; n++)
        if (facedata[loop.points[n]].fmFather==3)
          facedata[loop.points[n]].fmFather=0;
      current_face = fn11;
      continue;
    };
  }

  int *newFaces = new int[loop.npoints];
  npoints = 0;
  for (int n = 0 ; n < loop.npoints ; n++) {
    if (facedata[loop.points[n]].fmFather==1)
      newFaces[npoints++]=loop.points[n];
  }
  delete [] loop.points;
  loop.points = newFaces;
  loop.maxpoints=loop.npoints;
  loop.npoints = npoints;

  return _OrderLoop(loop);
}

void FastLoop::_InitDefect() {
  Init();
  for (int n = 0 ; n < ndefect_faces ; n++) {
    FaceData * face = &facedata[defect_faces[n]];
    face->fmState = eFar;
    face->fccState = eFar;
    face->border=0;
    for (int i = 0 ; i < 3 ; i++) {
      vertexdata[face->v[i]].fmState=eFar;
    }
  }
}

void FastLoop::SetDefectList(int nfaces, int *list_of_faces) {
  defect_faces = list_of_faces;
  ndefect_faces = nfaces;
  // setting up the correct flags to all the faces
  for (int n = 0 ; n < surface->nfaces ; n++) {
    FaceData *face = &facedata[n];
    face->border=1;
    face->fccState=eForbidden;
    face->fmState=eForbidden;
    for (int i = 0 ; i < 3 ; i++)
      vertexdata[facedata[n].v[i]].fmState = eForbidden;
  }
}

// return NULL if fails to find loop
Loop* FastLoop::FindLoop(int seed) {
#if PRINT_ERROR_MODE
  cout << endl << "seed point is : " << seed << " ";
#endif
  //cout << " " << seed << ".";
  while (seed >= surface->nfaces || seed < 0) {
    seed = Random(surface->nfaces);
  }

  _InitDefect();

  //setting up the seed
  SetSeed(seed);

  bool loop_not_found=true;
  int stopping_face=-1, conflicting_face=-1;

  while (loop_not_found) {

    ///////////////////////////////////////////////////////////////
    // fast marching on the eTrial + eFar faces
    conflicting_face = _Run(stopping_face);

    if (conflicting_face == -1) {
#if PRINT_ERROR_MODE
      cout << endl << "No loop has been found!!!!" << endl;
#endif
      return NULL;
    };
#if PRINT_ERROR_MODE
    cout << "conflict " << conflicting_face << " stopping " << stopping_face << endl;
#endif
    // at this point, conflicting_face generated a clash with the existing stopping_face
    FM_trial_heap->pop();
    _AddAliveFace(conflicting_face);
    for (int i = 0 ; i < 3 ; i++)
      _UpdateFace(facedata[conflicting_face].f[i],conflicting_face);
    ////////////////////////////////////////////////////////////////
    //now check if the remaining faces (eFar+eTrial) are connected

    //first find the two neighboring faces
#if PRINT_ERROR_MODE
    cout << endl << "Init segmentation " ;
#endif
    _InitSegmentation(conflicting_face);

    //start the segmentation
#if PRINT_ERROR_MODE
    cout << endl << "Starting segmentation " ;
#endif
    int which_segment = _FastSegmentation();

    if (which_segment<0) //we have found a path
      break;

    //check if there is a path outside  //here
    if (segments[which_segment].GetMark()) { //this segment is touching the border
      bool is_loop=false;
      int other_segment=(which_segment+1)%2;
      if (segments[other_segment].GetMark()) //we have a loop
        is_loop=true;
      else {
        //check if remaining faces are neighbors with outside
        while (nsegments[other_segment]) {
          int fn = FCC_trial_heap->top(); //next face in list
          FCC_trial_heap->pop();
          int label = facedata[fn].fccLabel;
          ASSERT(label==other_segment);
          nsegments[label]--;
          segments[label].AddPoint(fn);
          for (int i = 0 ; i < 3 ; i++) {
            if (facedata[facedata[fn].f[i]].border) {
              is_loop=true;
              break;
            }
          }
          if (is_loop) break;
          facedata[fn].fccState=eAlive;
          for (int i = 0 ; i < 3 ; i++)
            _UpdateSegmentFace(facedata[fn].f[i],fn);
        }

      }
      if (is_loop) {
#if PRINT_ERROR_MODE
        cout << "Found only One Single Loop!!!!" << endl;
#endif
        Loop *loops = new Loop[2];
		if(!_ExtractFirstLoop(loops[0],conflicting_face,stopping_face)){
			delete[] loops;
			return NULL;
		}
#if PRINT_ERROR_MODE
        loops[0].Print();
#endif
        _SimplifyLoop(loops[0]);
#if PRINT_ERROR_MODE
        cout << endl << "loop 0 : after simplification: ";
        loops[0].Print();
#endif
        //second loop is empty nothing
        return loops;
      }
    }
    _UpdateSegmentFaces(which_segment);

  }
#if PRINT_ERROR_MODE
  cout << endl << " Path has been found!"<<endl;
#endif
  //extracting path
  Loop *loops = new Loop[2];
  if(!_ExtractFirstLoop(loops[0],conflicting_face,stopping_face)){
	delete [] loops;
	return NULL;
  }
#if PRINT_ERROR_MODE
  loops[0].Print();
#endif
  _ExtractSecondLoop(loops[1],conflicting_face);
#if PRINT_ERROR_MODE
  cout << endl;
  loops[1].Print();
#endif
  _SimplifyLoop(loops[0]);
  _SimplifyLoop(loops[1]);
#if PRINT_ERROR_MODE
  cout << endl << "loop 0 : after simplification: ";
  loops[0].Print();
  cout << endl << "loop 1 : after simplification: ";
  loops[1].Print();
#endif
  return loops;
}


double FastLoop::_GetLoopLength(Loop &loop) {
  double length=0;
  for (int n = 0 ; n < loop.npoints ; n++) {
    int f0 = loop.points[n],f1;
    if (n==loop.npoints-1) f1=loop.points[0];
    else f1=loop.points[n+1];
    length += _Distance(f0,f1);
  }
  return length;
}

void FastLoop::FindMinimalLoop(Loop & minimal_loop, int max_init_face , int nattempts) {

  double minimal_length=_GetLoopLength(minimal_loop);

  if (max_init_face == -1) max_init_face = surface->nfaces;
  if ( nattempts < 1 ) nattempts = 10;

  for (int n = 0 ; n < nattempts ; n++) {
    Loop *loop;

    int seed_face;
    //find a seed face
    do {
      seed_face = Random(max_init_face);
    } while (seed_face < 0 || seed_face >= max_init_face);
    //find a loop
    loop = FindLoop(seed_face);

    if (loop==0) continue; //next attempt

    int which_loop = 0 ;
    double length = _GetLoopLength(loop[0]), olength;
    if (loop[1].npoints) {
      olength = _GetLoopLength(loop[1]);
      if (olength < length) {
        which_loop=1;
        length=olength;
      }
    }

    if (minimal_loop.npoints == 0 || //loop[which_loop].npoints < minimal_loop.npoints){
        length < minimal_length) {
      minimal_length = length;
      minimal_loop = loop[which_loop];
    }
    delete [] loop;
  }
  // cout << "minimal seed is " << minimal_seed << " avec "
  //   << minimal_loop.npoints << " points  et length = " << minimal_length << " ("<< wloop << " -" << nfound << ")" << endl;

}
