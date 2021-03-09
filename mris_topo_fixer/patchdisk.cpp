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


#include "patchdisk.h"

PatchDisk::PatchDisk(void) {
  vtrans=0;
  ftrans=0;
}

PatchDisk::PatchDisk(int which_patch) {
  Create(which_patch);
}

void PatchDisk::_Alloc(int which_patch) {
  Face *face;

  switch (which_patch) {
  case 0: //super small surface
    disk.Expand(3,1);
    disk.nvertices=3;
    disk.nfaces = 1;
    face = &disk.faces[0];
    face->v[0] = 0;
    face->v[1] = 2;
    face->v[2] = 1;
    break;
  case 1:
    disk.Expand(5,4);
    disk.nvertices=5;
    disk.nfaces = 4;
    face = &disk.faces[0];
    face->v[0] = 0;
    face->v[1] = 2;
    face->v[2] = 1;
    face = &disk.faces[1];
    face->v[0] = 0;
    face->v[1] = 3;
    face->v[2] = 2;
    face = &disk.faces[2];
    face->v[0] = 0;
    face->v[1] = 4;
    face->v[2] = 3;
    face = &disk.faces[3];
    face->v[0] = 0;
    face->v[1] = 1;
    face->v[2] = 4;
    break;
  case 2: //small surface
    disk.Expand(21,28);
    disk.nvertices=21;
    disk.nfaces=28;
    face = &disk.faces[0];
    face->v[0] = 2;
    face->v[1] = 3;
    face->v[2] = 4;
    face = &disk.faces[1];
    face->v[0] = 5;
    face->v[1] = 0;
    face->v[2] = 6;
    face = &disk.faces[2];
    face->v[0] = 0;
    face->v[1] = 1;
    face->v[2] = 6;
    face = &disk.faces[3];
    face->v[0] = 4;
    face->v[1] = 3;
    face->v[2] = 7;
    face = &disk.faces[4];
    face->v[0] = 11;
    face->v[1] = 9;
    face->v[2] = 10;
    face = &disk.faces[5];
    face->v[0] = 14;
    face->v[1] = 12;
    face->v[2] = 13;
    face = &disk.faces[6];
    face->v[0] = 8;
    face->v[1] = 12;
    face->v[2] = 14;
    face = &disk.faces[7];
    face->v[0] = 15;
    face->v[1] = 9;
    face->v[2] = 11;
    face = &disk.faces[8];
    face->v[0] = 16;
    face->v[1] = 1;
    face->v[2] = 2;
    face = &disk.faces[9];
    face->v[0] = 4;
    face->v[1] = 16;
    face->v[2] = 2;
    face = &disk.faces[10];
    face->v[0] = 17;
    face->v[1] = 5;
    face->v[2] = 6;
    face = &disk.faces[11];
    face->v[0] = 16;
    face->v[1] = 18;
    face->v[2] = 17;
    face = &disk.faces[12];
    face->v[0] = 16;
    face->v[1] = 17;
    face->v[2] = 6;
    face = &disk.faces[13];
    face->v[0] = 16;
    face->v[1] = 6;
    face->v[2] = 1;
    face = &disk.faces[14];
    face->v[0] = 19;
    face->v[1] = 18;
    face->v[2] = 16;
    face = &disk.faces[15];
    face->v[0] = 19;
    face->v[1] = 16;
    face->v[2] = 4;
    face = &disk.faces[16];
    face->v[0] = 19;
    face->v[1] = 4;
    face->v[2] = 7;
    face = &disk.faces[17];
    face->v[0] = 8;
    face->v[1] = 19;
    face->v[2] = 7;
    face = &disk.faces[18];
    face->v[0] = 10;
    face->v[1] = 5;
    face->v[2] = 17;
    face = &disk.faces[19];
    face->v[0] = 17;
    face->v[1] = 18;
    face->v[2] = 20;
    face = &disk.faces[20];
    face->v[0] = 17;
    face->v[1] = 20;
    face->v[2] = 11;
    face = &disk.faces[21];
    face->v[0] = 17;
    face->v[1] = 11;
    face->v[2] = 10;
    face = &disk.faces[22];
    face->v[0] = 20;
    face->v[1] = 18;
    face->v[2] = 19;
    face = &disk.faces[23];
    face->v[0] = 20;
    face->v[1] = 19;
    face->v[2] = 14;
    face = &disk.faces[24];
    face->v[0] = 20;
    face->v[1] = 14;
    face->v[2] = 13;
    face = &disk.faces[25];
    face->v[0] = 19;
    face->v[1] = 8;
    face->v[2] = 14;
    face = &disk.faces[26];
    face->v[0] = 15;
    face->v[1] = 11;
    face->v[2] = 20;
    face = &disk.faces[27];
    face->v[0] = 20;
    face->v[1] = 13;
    face->v[2] = 15;
    break;
  case 3: //medium surface
  default:
    disk.Expand(33,52);
    disk.nvertices=33;
    disk.nfaces=52;
    face = &disk.faces[0];
    face->v[0] = 12;
    face->v[1] = 0;
    face->v[2] = 1;
    face = &disk.faces[1];
    face->v[0] = 1;
    face->v[1] = 2;
    face->v[2] = 14;
    face = &disk.faces[2];
    face->v[0] = 1;
    face->v[1] = 14;
    face->v[2] = 13;
    face = &disk.faces[3];
    face->v[0] = 1;
    face->v[1] = 13;
    face->v[2] = 12;
    face = &disk.faces[4];
    face->v[0] = 14;
    face->v[1] = 2;
    face->v[2] = 3;
    face = &disk.faces[5];
    face->v[0] = 14;
    face->v[1] = 3;
    face->v[2] = 15;
    face = &disk.faces[6];
    face->v[0] = 14;
    face->v[1] = 15;
    face->v[2] = 16;
    face = &disk.faces[7];
    face->v[0] = 4;
    face->v[1] = 15;
    face->v[2] = 3;
    face = &disk.faces[8];
    face->v[0] = 17;
    face->v[1] = 5;
    face->v[2] = 0;
    face = &disk.faces[9];
    face->v[0] = 17;
    face->v[1] = 0;
    face->v[2] = 12;
    face = &disk.faces[10];
    face->v[0] = 17;
    face->v[1] = 12;
    face->v[2] = 18;
    face = &disk.faces[11];
    face->v[0] = 12;
    face->v[1] = 13;
    face->v[2] = 18;
    face = &disk.faces[12];
    face->v[0] = 16;
    face->v[1] = 15;
    face->v[2] = 19;
    face = &disk.faces[13];
    face->v[0] = 4;
    face->v[1] = 6;
    face->v[2] = 20;
    face = &disk.faces[14];
    face->v[0] = 4;
    face->v[1] = 20;
    face->v[2] = 19;
    face = &disk.faces[15];
    face->v[0] = 4;
    face->v[1] = 19;
    face->v[2] = 15;
    face = &disk.faces[16];
    face->v[0] = 7;
    face->v[1] = 5;
    face->v[2] = 17;
    face = &disk.faces[17];
    face->v[0] = 7;
    face->v[1] = 17;
    face->v[2] = 22;
    face = &disk.faces[18];
    face->v[0] = 7;
    face->v[1] = 22;
    face->v[2] = 21;
    face = &disk.faces[19];
    face->v[0] = 23;
    face->v[1] = 21;
    face->v[2] = 22;
    face = &disk.faces[20];
    face->v[0] = 26;
    face->v[1] = 24;
    face->v[2] = 25;
    face = &disk.faces[21];
    face->v[0] = 20;
    face->v[1] = 6;
    face->v[2] = 8;
    face = &disk.faces[22];
    face->v[0] = 20;
    face->v[1] = 8;
    face->v[2] = 24;
    face = &disk.faces[23];
    face->v[0] = 20;
    face->v[1] = 24;
    face->v[2] = 26;
    face = &disk.faces[24];
    face->v[0] = 9;
    face->v[1] = 7;
    face->v[2] = 21;
    face = &disk.faces[25];
    face->v[0] = 27;
    face->v[1] = 10;
    face->v[2] = 9;
    face = &disk.faces[26];
    face->v[0] = 27;
    face->v[1] = 9;
    face->v[2] = 21;
    face = &disk.faces[27];
    face->v[0] = 27;
    face->v[1] = 21;
    face->v[2] = 23;
    face = &disk.faces[28];
    face->v[0] = 11;
    face->v[1] = 10;
    face->v[2] = 27;
    face = &disk.faces[29];
    face->v[0] = 11;
    face->v[1] = 27;
    face->v[2] = 25;
    face = &disk.faces[30];
    face->v[0] = 11;
    face->v[1] = 25;
    face->v[2] = 24;
    face = &disk.faces[31];
    face->v[0] = 24;
    face->v[1] = 8;
    face->v[2] = 11;
    face = &disk.faces[32];
    face->v[0] = 28;
    face->v[1] = 13;
    face->v[2] = 14;
    face = &disk.faces[33];
    face->v[0] = 16;
    face->v[1] = 28;
    face->v[2] = 14;
    face = &disk.faces[34];
    face->v[0] = 29;
    face->v[1] = 17;
    face->v[2] = 18;
    face = &disk.faces[35];
    face->v[0] = 28;
    face->v[1] = 30;
    face->v[2] = 29;
    face = &disk.faces[36];
    face->v[0] = 28;
    face->v[1] = 29;
    face->v[2] = 18;
    face = &disk.faces[37];
    face->v[0] = 28;
    face->v[1] = 18;
    face->v[2] = 13;
    face = &disk.faces[38];
    face->v[0] = 31;
    face->v[1] = 30;
    face->v[2] = 28;
    face = &disk.faces[39];
    face->v[0] = 31;
    face->v[1] = 28;
    face->v[2] = 16;
    face = &disk.faces[40];
    face->v[0] = 31;
    face->v[1] = 16;
    face->v[2] = 19;
    face = &disk.faces[41];
    face->v[0] = 20;
    face->v[1] = 31;
    face->v[2] = 19;
    face = &disk.faces[42];
    face->v[0] = 22;
    face->v[1] = 17;
    face->v[2] = 29;
    face = &disk.faces[43];
    face->v[0] = 29;
    face->v[1] = 30;
    face->v[2] = 32;
    face = &disk.faces[44];
    face->v[0] = 29;
    face->v[1] = 32;
    face->v[2] = 23;
    face = &disk.faces[45];
    face->v[0] = 29;
    face->v[1] = 23;
    face->v[2] = 22;
    face = &disk.faces[46];
    face->v[0] = 32;
    face->v[1] = 30;
    face->v[2] = 31;
    face = &disk.faces[47];
    face->v[0] = 32;
    face->v[1] = 31;
    face->v[2] = 26;
    face = &disk.faces[48];
    face->v[0] = 32;
    face->v[1] = 26;
    face->v[2] = 25;
    face = &disk.faces[49];
    face->v[0] = 31;
    face->v[1] = 20;
    face->v[2] = 26;
    face = &disk.faces[50];
    face->v[0] = 27;
    face->v[1] = 23;
    face->v[2] = 32;
    face = &disk.faces[51];
    face->v[0] = 32;
    face->v[1] = 25;
    face->v[2] = 27;
    break;
  }
  disk.InitSurface();
  vtrans = new int[disk.nvertices];
  ftrans = new int[disk.nfaces];
}


void PatchDisk::_Init(void) {
  //finding the init_ring of faces in the disk
  for (int n = 0 ; n < disk.nvertices ; n++)
    disk.vertices[n].marked=0;

  int first_v1 = -1 , first_v2 = -1;
  for (int n = 0 ; n < disk.nfaces ; n++)
    for (int i = 0 ; i < 3 ; i++) {
      if (disk.faces[n].f[i]==-1) { //we have a border face
        disk.vertices[disk.faces[n].v[i]].marked=1;
        if (first_v1==-1) first_v1 = disk.faces[n].v[i];
        disk.vertices[disk.faces[n].v[(i+1)%3]].marked=1;
        if (first_v2==-1) first_v2 = disk.faces[n].v[(i+1)%3];
      }
    }

  init_ring.AddPoint(first_v1);
  disk.vertices[first_v1].marked=2;
  init_ring.AddPoint(first_v2);
  disk.vertices[first_v2].marked=2;

  int current_vertex = first_v2;
  bool found=true;
  while (found) {
    found=false;
    Vertex *v=&disk.vertices[current_vertex];
    for (int p = 0 ; p < v->vnum ; p++) {
      if (disk.vertices[v->v[p]].marked==1) {
        found = true;
        current_vertex = v->v[p];
        init_ring.AddPoint(v->v[p]);
        disk.vertices[v->v[p]].marked=2;
        break;
      }
    }
  }
  ring.Alloc(init_ring.npoints);
  ring.npoints = init_ring.npoints;
}

PatchDisk::~PatchDisk(void) {
  if (vtrans) delete [] vtrans;
  if (ftrans) delete [] ftrans;
}

void PatchDisk::Init() {
  for (int n = 0 ; n < ring.npoints ; n++)
    ring.points[n] = init_ring.points[n];
}

void PatchDisk::Create(int which_patch) {
  _Alloc(which_patch);
  _Init();
}
