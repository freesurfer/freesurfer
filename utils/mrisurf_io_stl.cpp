#include "mrisurf.h"
#include "mrisurf_metricProperties.h"

// STLs are simple files that only record a collection of triangles, so when reading,
// we  have to build a list of unique vertices on our own with a simple hash table.

#define hashTableSize (1024 * 1024)
#define someKeysSize 100000

// An xyz hashing function that preserves many of the bits in the float and doesn't care much
// about their range. Cost is not very important compared to all the other work done.
static int pointHash(float x, float y, float z)
{
  int hx = 1000.0f * logf(abs(x) + 1.0f);
  int hy = 1000.0f * logf(abs(y) + 1.0f);
  int hz = 1000.0f * logf(abs(z) + 1.0f);
  int hash = (hx*2017) ^ (hy*1029) ^ (hz*1159);  // No good reason behind these multipliers
  return hash % hashTableSize;
}


// STLs can be either binary or ascii, so the STLMesher class can be used as a universal
// means of registering vertices and faces.
class STLMesher
{
  struct Key  { struct Key*  next; float x, y, z; int vertexNo; };
  struct Keys { struct Keys* next; Key someKeys[someKeysSize]; };
  struct XYZ  { float x,y,z; };

  unsigned int vidx = 0;  // current vertex index
  unsigned int fidx = 0;  // current face index

  int nextVertexNo = 0;
  int faceVertexNo = 3;  // to make sure there are never more than 3 vertices per face

  std::vector<Key*> hashTable;
  std::vector<Key*> vertices;
  std::vector<XYZ>  normals;

  Keys* keys = nullptr;
  int keysInUse = someKeysSize;

public:
  // initialize the vertices and normals vectors with an arbitrary starting size - they
  // can always be extended later on
  STLMesher() : hashTable(hashTableSize, nullptr), vertices(500000), normals(150000) {}

  ~STLMesher() {
    while (keys) {
      Keys* next = keys->next;
      delete keys;
      keys = next;
    }
  }

  // adds a vertex to the mesh if it has not already been added
  void addVertex(float &x, float &y, float &z)
  {
    // hash the vertex position - those already seen will reuse their previously assigned vertex idx
    Key** keyPtrPtr = &hashTable[pointHash(x, y, z)];
    Key* keyPtr;

    // ensure vertex is unique (or not)
    for (; (keyPtr = *keyPtrPtr); keyPtrPtr = &keyPtr->next) {
      if (keyPtr->x == x && keyPtr->y == y && keyPtr->z == z) break;
    }

    if (!keyPtr) {
      if (someKeysSize == keysInUse) {
        Keys* moreKeys = new Keys();
        moreKeys->next = keys;
        keys = moreKeys;
        keysInUse = 0;
      }
      keyPtr = &keys->someKeys[keysInUse++];
      keyPtr->next = nullptr;
      keyPtr->x = x;
      keyPtr->y = y;
      keyPtr->z = z;
      keyPtr->vertexNo = nextVertexNo++;
      *keyPtrPtr = keyPtr;
    }

    // sanity check on the current face geometry
    if (faceVertexNo++ > 2) fs::fatal() << "STL file does not appear to represent a triangular mesh";

    // store the vertex (resize if necessary)
    if (vidx == vertices.size()) vertices.resize(vertices.size() * 2);
    vertices[vidx++] = keyPtr;
  }

  // initiates a new face (represented by its normal) in the mesh
  void addFaceNormal(float &x, float &y, float &z)
  {
    // resize if necessary
    if (fidx == normals.size()) normals.resize(normals.size() * 2);
    normals[fidx].x = x;
    normals[fidx].y = y;
    normals[fidx].z = z;
    fidx++;
    faceVertexNo = 0;
  }

  // build the final surface structure once all data has been read
  MRIS* buildSurface()
  {
    MRIS *mris = MRISalloc(nextVertexNo, fidx);
    mris->type = MRIS_TRIANGULAR_SURFACE;

    // set corresponding face vertices
    int nv = 0;
    for (int faceNo = 0; faceNo < mris->nfaces; faceNo++) {

      FACE* face = &mris->faces   [faceNo];
      XYZ*  xyz  = &normals[faceNo];

      setFaceNorm(mris, faceNo, xyz->x, xyz->y, xyz->z);

      for (int faceVertexNo = 0; faceVertexNo < 3; faceVertexNo++) {
        int undupVertexNo = 3*faceNo + faceVertexNo;
        Key* keyPtr = vertices[undupVertexNo];
        int vertexNo = keyPtr->vertexNo;
        if (vertexNo == nv) {
          MRISsetXYZ(mris,vertexNo, keyPtr->x, keyPtr->y, keyPtr->z);
          nv++;
        }
        face->v[faceVertexNo] = vertexNo;
      }
    }

    // count number of faces that each vertex is part of
    for (int fno = 0; fno < mris->nfaces; fno++) {
      FACE *face = &mris->faces[fno];
      for (int fvno = 0; fvno < VERTICES_PER_FACE; fvno++) {
        VERTEX_TOPOLOGY * const v = &mris->vertices_topology[face->v[fvno]];
        v->num++;
        addVnum(mris,face->v[fvno],2);
      }
    }

    // alloc memory for neighbor list of each vertex
    for (int vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vno];
      v->v = (int *)calloc(v->vnum, sizeof(int));
      if (!v->v) ErrorExit(ERROR_NOMEMORY, "MRISreadSTLfile: could not allocate %dth vertex list.", vno);
      clearVnum(mris, vno);
    }

    // build list of neighbors
    for (int fno = 0; fno < mris->nfaces; fno++) {
      FACE *face = &mris->faces[fno];
      for (int fvno = 0; fvno < 3; fvno++) {
        VERTEX_TOPOLOGY * const v = &mris->vertices_topology[face->v[fvno]];

        // add an edge to other 2 vertices if not already in list
        for (int fvno2 = 0; fvno2 < VERTICES_PER_FACE; fvno2++) {
          
          // don't connect vertex to itself
          if (fvno2 == fvno) continue;

          // check to make sure it's not a duplicate
          int vn = mris->faces[fno].v[fvno2];
          for (int n2 = 0; n2 < v->vnum; n2++) {
            if (v->v[n2] == vn) {
              vn = -1;  // mark it as a duplicate
              break;
            }
          }
          if (vn >= 0) v->v[vnumAdd(mris,face->v[fvno],1)] = vn;
        }
      }
    }

    // allocate face arrays in vertices
    for (int vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
      vt->f = (int *)calloc(vt->num, sizeof(int));
      if (!vt->f) ErrorExit(ERROR_NO_MEMORY, "MRISreadSTLfileICOread: could not allocate %d faces", vt->num);
      vt->n = (unsigned char *)calloc(vt->num, sizeof(unsigned char));
      if (!vt->n) ErrorExit(ERROR_NO_MEMORY, "MRISreadSTLfile: could not allocate %d nbrs", vt->n);
      vt->num = 0;  // for use as counter in next section
      vt->vtotal = vt->vnum;
    }

    // fill in face indices in vertex structures
    for (int fno = 0; fno < mris->nfaces; fno++) {
      FACE* face = &mris->faces[fno];
      for (int fvno = 0; fvno < 3; fvno++) {
        VERTEX_TOPOLOGY * const v = &mris->vertices_topology[face->v[fvno]];
        v->n[v->num] = fvno;
        v->f[v->num++] = fno;
      }
    }

    return mris;
  }
};


// Reads an ascii-formatted STL file into a surface.
static MRIS *mrisReadAsciiSTL(const std::string &filename)
{
  STLMesher mesher;
  int solids = 0;

  FILE *stlfile = fopen(filename.c_str(), "r");
  char line[STRLEN], *cp;
  while ((cp = fgetl(line, STRLEN, stlfile))) {

    // For efficiency purposes, we'll only look at the first character
    // of each line to determine what it is. This should not be problematic
    // for a valid STL file.
    switch (*cp) {

      // vertex
      case 'v': {
        float x, y, z;
        sscanf(cp, "vertex %f %f %f", &x, &y, &z);
        mesher.addVertex(x, y, z);
        continue;
      } break;

      // facet
      case 'f': {
        float x, y, z;
        sscanf(cp, "facet normal %f %f %f", &x, &y, &z);
        mesher.addFaceNormal(x, y, z);
        continue;
      } break;

      // solid
      case 's': {
        if (++solids > 1) {
          fs::error() << "STL file " << filename << " contains multiple solids!";
          return nullptr;
        }
        continue;
      } break;

      case 'e': { continue; } break;  // end specifiers
      case 'o': { continue; } break;  // outer loop
    }

    fs::error() << "ASCII STL file " << filename << " contains unrecognized line: " << line;
    return nullptr;
  }

  return mesher.buildSurface();
}


// Reads a single float from an ifstream.
static float readFloat(std::ifstream &stream)
{
  float value;
  stream.read(reinterpret_cast<char*>(&value), 4);
  return value;
}


// Reads a binary-formatted STL file into a surface.
static MRIS* mrisReadBinarySTL(const std::string &filename) 
{
  STLMesher mesher;

  // load and read 80-byte header (not used)
  std::ifstream stlfile(filename, std::ios::in | std::ios::binary);
  char header[80] = "";
  stlfile.read(header, 80);

  // read total number of triangles
  unsigned int ntriangles;
  stlfile.read(reinterpret_cast<char*>(&ntriangles), 4);

  // for each triangle read the face normal and the associated vertex positions
  for (unsigned int tri = 0 ; tri < ntriangles ; tri++) {

    // face normal
    float x = readFloat(stlfile);
    float y = readFloat(stlfile);
    float z = readFloat(stlfile);
    mesher.addFaceNormal(x, y, z);

    // face vertices
    for (int vnum = 0 ; vnum < 3; vnum++) {
      float x = readFloat(stlfile);
      float y = readFloat(stlfile);
      float z = readFloat(stlfile);
      mesher.addVertex(x, y, z);
    }

    // read 2-byte attribute (not used)
    char attribute[2];
    stlfile.read(attribute, 2);
  }

  return mesher.buildSurface();
}


MRIS* mrisReadSTL(const std::string &filename)
{
  // STL files can be binary or ascii - the easiest way to check is to
  // assume binary and read the number of triangles specified in the header.
  // If the file is binary, it's total size will be equal to (84 + ntriangles * 50) bytes

  // open file
  std::ifstream stlfile(filename, std::ios::in | std::ios::binary);
  if (!stlfile) return nullptr;

  // skip the first 80 bytes
  char header[80] = "";
  stlfile.read(header, 80);

  // read the number of triangles specified and calculate the predicted (binary) filesize
  unsigned int ntriangles;
  stlfile.read(reinterpret_cast<char*>(&ntriangles), 4);
  unsigned int fileSizeIfBinary = 84 + ntriangles * 50;

  // extract the true filesize
  stlfile.seekg(0, std::ios::end);
  unsigned int fileSize = stlfile.tellg();
  stlfile.close();

  // if the true filesize matches the predicted size, it's in binary format
  if (fileSize == fileSizeIfBinary) {
    return mrisReadBinarySTL(filename);
  } else {
    return mrisReadAsciiSTL(filename);
  }
}


// Writes a surface into a binary-formatted STL file.
bool mrisWriteSTL(const MRIS* mris, const std::string &filename)
{
  std::ofstream stlfile(filename.c_str(), std::ios::out | std::ios::binary);
  if (!stlfile) return false;

  // write an arbitrary 80 byte header
  char header[80] = "freesurfer mesh";
  stlfile.write(header, 80);

  // specify the total number of triangles
  unsigned int ntris = mris->nfaces;
  stlfile.write((char*)&ntris, 4);

  for (int fno = 0; fno < mris->nfaces; fno++) {

    // write the face normal
    FaceNormCacheEntry const * const norm = getFaceNorm(mris, fno);
    stlfile.write((char*)&norm->nx, 4);
    stlfile.write((char*)&norm->ny, 4);
    stlfile.write((char*)&norm->nz, 4);

    // write the face vertices
    FACE *face = &mris->faces[fno];
    for (int vno = 0; vno < 3; vno++) {
      VERTEX *vtx = &mris->vertices[face->v[vno]];
      stlfile.write((char*)&vtx->x, 4);
      stlfile.write((char*)&vtx->y, 4);
      stlfile.write((char*)&vtx->z, 4);
    }

    // attribute (not used)
    stlfile.write("0", 2);
  }

  return true;
}
