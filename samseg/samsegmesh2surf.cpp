#include "itkImageFileReader.h"  // needed for itksys::SystemTools::FileExists()

#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshCollection.h"

#include "mrisurf.h"
#include "matrix.h"
#include "mri.h"

int __Mesh2Surf(const char *mesh, const char *nii, const char *surfFile, const char *overlayVol);

// samsegmesh2surf mesh template_nii surf overlayvol
int main(int argc, char** argv)
{
  if (argc < 5)
  {
    printf("Usage: samsegmesh2surf <mesh> <template> <surf> <overlay-vol>\n");
    exit(1);
  }

  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 5; argumentNumber++ ) 
  {
    inputParserStream << argv[ argumentNumber ] << " ";
  }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  std::string  meshCollectionFile;
  std::string  template_nii;
  std::string  surfFile;
  std::string  overlayVol;
  inputStream >> meshCollectionFile >> template_nii >> surfFile >> overlayVol;

  std::cout << "samsegmesh2surf Command line params:" << std::endl;
  std::cout << "  meshCollectionFile:         " << meshCollectionFile << std::endl;
  std::cout << "  template:                   " << template_nii << std::endl;
  std::cout << "  surface file:               " << surfFile << std::endl;
  std::cout << "  overlay volume:             " << overlayVol << std::endl;

  __Mesh2Surf(meshCollectionFile.c_str(), template_nii.c_str(), surfFile.c_str(), overlayVol.c_str());

  return 0;
}


int __Mesh2Surf(const char *meshfile, const char *nii, const char *surfFile, const char *overlayVol)
{
  // read mesh collection
  kvl::AtlasMeshCollection::Pointer  meshStartCollection = 0;
  if ( itksys::SystemTools::FileExists( meshfile, true ) )
  {
    meshStartCollection = kvl::AtlasMeshCollection::New();
    if ( !meshStartCollection->Read( meshfile ) )
    {
      std::cerr << "Couldn't read mesh from file " << meshfile << std::endl;
      return -1;
    } 
    else
    {
      std::cout << "reading mesh " << meshfile << std::endl;
    }
  }

  /* read template, calculate vox2tkras matrix
   * /usr/local/freesurfer/dev/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/template.nii
   * dimensions: 181 x 217 x 181
   *
   * >/usr/local/freesurfer/dev/bin/mri_info template.nii --vox2ras
   *   -1.00000   -0.00000    0.00000   90.00000 
   *   -0.00000    1.00000   -0.00000 -126.00000 
   *    0.00000    0.00000    1.00000  -72.00000 
   *    0.00000    0.00000    0.00000    1.00000 
   * >/usr/local/freesurfer/dev/bin/mri_info template.nii --vox2ras-tkr
   *   -1.00000    0.00000    0.00000   90.50000 
   *    0.00000    0.00000    1.00000  -90.50000 
   *    0.00000   -1.00000    0.00000  108.50000 
   *    0.00000    0.00000    0.00000    1.00000 
   */
  MRI *niiVol = MRIread(nii);
  MATRIX *vox2tkras = niiVol->get_Vox2TkregRAS();

  kvl::AtlasMesh::ConstPointer constMesh = meshStartCollection->GetReferenceMesh();
  const int numClasses = constMesh->GetPointData()->Begin().Value().m_Alphas.Size();
  const int numMeshPoints = constMesh->GetPoints()->Size();

  // count triangle faces
  kvl::AtlasMesh::Pointer mutableMesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );
  kvl::AtlasMesh::CellsContainer::ConstIterator cellsIt    = constMesh->GetCells()->Begin();
  kvl::AtlasMesh::CellsContainer::ConstIterator cellsItEnd = constMesh->GetCells()->End();
  int num_faces = 0;
  for ( ;  cellsIt != cellsItEnd; ++cellsIt )
    if ( cellsIt.Value()->GetType() == kvl::AtlasMesh::CellType::TRIANGLE_CELL )
      num_faces++;

  // the surface codes are borrowed from mrisReadGIFTIdanum()
  // create MRIS structure with given number of vertices and faces
  const int num_vertices = numMeshPoints;  //constMesh->GetPoints()->Size();
  printf("num_vertices = %d, num_faces = %d\n", num_vertices, num_faces);
  MRIS *Surf = MRISalloc(num_vertices, num_faces);  //MRISoverAlloc(nVertices, 0, nVertices, 0);
  Surf->type = MRIS_TRIANGULAR_SURFACE; // defined in mrisurf.h
  Surf->hemisphere = NO_HEMISPHERE;
  Surf->useRealRAS = 0;

  /****** begin setting surface volume geometry ******/
  // Surf->vg & Surf->xctr, Surf->yctr, Surf->zctr
  printf("[INFO] set Surf->vg ...\n");
  Surf->vg.width  = niiVol->width;
  Surf->vg.height = niiVol->height;
  Surf->vg.depth  = niiVol->depth;
  Surf->vg.xsize  = niiVol->xsize; Surf->vg.ysize = niiVol->ysize; Surf->vg.zsize = niiVol->zsize;
  Surf->vg.x_r = niiVol->x_r; Surf->vg.x_a = niiVol->x_a; Surf->vg.x_s = niiVol->x_s;
  Surf->vg.y_r = niiVol->y_r; Surf->vg.y_a = niiVol->y_a; Surf->vg.y_s = niiVol->y_s;
  Surf->vg.z_r = niiVol->z_r; Surf->vg.z_a = niiVol->z_a; Surf->vg.z_s = niiVol->z_s;
  Surf->vg.c_r = niiVol->c_r; Surf->vg.c_a = niiVol->c_a; Surf->vg.c_s = niiVol->c_s;

  Surf->vg.valid = 1;
  /****** end setting surface volume geometry ******/

  // create MRI structure for overlay
  MRI *Vol = new MRI({numMeshPoints, 1, 1, numClasses}, MRI_FLOAT);

  kvl::AtlasMesh::PointsContainer::Iterator pointsIt    = mutableMesh->GetPoints()->Begin(); 
  kvl::AtlasMesh::PointsContainer::Iterator pointsItEnd = mutableMesh->GetPoints()->End();
  for (; pointsIt != pointsItEnd; pointsIt++)
  {
    // (c, r, s)
    kvl::AtlasMesh::PointIdentifier pointid   = pointsIt.Index();
    kvl::AtlasMesh::PointType       thispoint = pointsIt.Value();
    int c = floor(thispoint[0]);
    int r = floor(thispoint[1]);
    int s = floor(thispoint[2]); 

    MATRIX *CRS = MatrixAlloc(4, 1, MATRIX_REAL);
    CRS->rptr[1][1] = c;
    CRS->rptr[2][1] = r;
    CRS->rptr[3][1] = s;
    CRS->rptr[4][1] = 1;

    // Convert the CRS to tkRAS
    MATRIX *tkRAS = MatrixAlloc(4, 1, MATRIX_REAL);
    tkRAS->rptr[4][1] = 1;
    tkRAS = MatrixMultiply(vox2tkras, CRS, tkRAS);

    // set surface vertex xyz
    MRISsetXYZ(Surf, pointid, tkRAS->rptr[1][1], tkRAS->rptr[2][1], tkRAS->rptr[3][1]);
    Surf->vertices[pointid].origarea = -1;

    // assign alphas to overlay volume
    for (int nClass = 0; nClass < numClasses; nClass++)
    {
      float alpha = constMesh->GetPointData()->ElementAt(pointid).m_Alphas[nClass];
      MRIsetVoxVal(Vol, pointid, 0, 0, nClass, alpha);
    }
  }  
  //printf("[INFO] Done setting vertex xyz\n");
  fflush(stdout);
  mrisComputeSurfaceDimensions(Surf);
  //printf("[INFO] Done mrisComputeSurfaceDimensions()\n");
  fflush(stdout);

  //int face_index = 0;
  cellsIt = mutableMesh->GetCells()->Begin();
  for (int face_index = 0; cellsIt != cellsItEnd; ++cellsIt)
  {
    if ( cellsIt.Value()->GetType() == kvl::AtlasMesh::CellType::TRIANGLE_CELL )
    {
      kvl::AtlasMesh::CellType::PointIdIterator  pit = cellsIt.Value()->PointIdsBegin();
      kvl::AtlasMesh::PointIdentifier  pointId0 = *pit;
      pit++;
      kvl::AtlasMesh::PointIdentifier  pointId1 = *pit;
      pit++;
      kvl::AtlasMesh::PointIdentifier  pointId2 = *pit;

      Surf->faces[face_index].v[0] = pointId0;
      Surf->vertices_topology[pointId0].num++;
      Surf->faces[face_index].v[1] = pointId1;
      Surf->vertices_topology[pointId1].num++;
      Surf->faces[face_index].v[2] = pointId2;
      Surf->vertices_topology[pointId2].num++;

      face_index++;
    } // if (kvl::AtlasMesh::CellType::TRIANGLE_CELL)
  }
  //printf("[INFO] Done setting faces\n");
  fflush(stdout);

  // each vertex has a face list (faster than face list in some operations)
  for (int vertex_index = 0; vertex_index < num_vertices; vertex_index++)
  {
    Surf->vertices_topology[vertex_index].f = (int   *)calloc(Surf->vertices_topology[vertex_index].num, sizeof(int));
    Surf->vertices_topology[vertex_index].n = (uchar *)calloc(Surf->vertices_topology[vertex_index].num, sizeof(uchar));
    Surf->vertices_topology[vertex_index].num = 0;  // this gets re-calc'd next...
  }
  //printf("[INFO] Done calloc\n");
  fflush(stdout);

  for (int face_index = 0; face_index < Surf->nfaces; face_index++)
  {
    FACE *face = &Surf->faces[face_index];
    int n;
    for (n = 0; n < VERTICES_PER_FACE; n++)
      Surf->vertices_topology[face->v[n]].f[Surf->vertices_topology[face->v[n]].num++] = face_index;  // note that .num is auto-incremented!
  }
  //printf("[INFO] Done going through faces\n");
  fflush(stdout);

  for (int vertex_index = 0; vertex_index < num_vertices; vertex_index++) {
    int n, m;
    for (n = 0; n < Surf->vertices_topology[vertex_index].num; n++) {
      for (m = 0; m < VERTICES_PER_FACE; m++) {
        if (Surf->faces[Surf->vertices_topology[vertex_index].f[n]].v[m] == vertex_index) {
          Surf->vertices_topology[vertex_index].n[n] = m;
        }
      }
    }
  }
  //printf("[INFO] Done making face list\n");
  fflush(stdout);


  mrisCompleteTopology(Surf);
  //printf("[INFO] Done SurfCompleteTopology()\n");
  fflush(stdout);

  MRIScomputeNormals(Surf);
  //printf("[INFO] Done MRIScomputeNormals()\n");
  fflush(stdout);

  UpdateMRIS(Surf, NULL);

  //printf("[DEBUG] (Surf->dist_nsize=%x (%d)\n", Surf->dist_nsize, Surf->dist_nsize);

  printf("Writing Surface to %s ...\n", surfFile);
  int err = MRISwrite(Surf, surfFile);
  if (err)
  {
    printf("Error outputing surface %s\n", surfFile);
    return 1;
  }

  printf("Writing Overlay Volume to %s ...\n", overlayVol);
  err = MRIwrite(Vol, overlayVol);
  if (err)
  {
    printf("Error outputing volume %s\n", overlayVol);
    return 1;
  }

  return 0;
}

