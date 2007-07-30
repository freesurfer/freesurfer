#include "vtkFSSurfaceSource.h"
#include "vtkFSSurfaceWriter.h"

const char* kfnInput = "lh.vtkFSSurfaceWriterTestFile";
const char* kfnOutput = "lh.vtkFSSurfaceWriterTestFileOut";

char* Progname = "vtkFSSurfaceWriterTester";

int
main ( int argc, char** argv ) {


  // Make a source.
  vtkSmartPointer<vtkFSSurfaceSource> source = 
    vtkSmartPointer<vtkFSSurfaceSource>::New();
  source->MRISRead( kfnInput );
  source->Update();

  // Make a dest.
  vtkSmartPointer<vtkFSSurfaceWriter> writer =
    vtkSmartPointer<vtkFSSurfaceWriter>::New();
  writer->SetFileName( kfnOutput );
  writer->SetInputConnection( source->GetOutputPort() );
  writer->Update();

  MRIS* surf1 = source->GetMRIS();
  assert( surf1 );
  MRIS* surf2 = writer->GetMRIS();
  assert( surf2 );
  
  // Check the number of verts and faces.
  assert( surf1->nvertices == surf2->nvertices );
  assert( surf1->nfaces == surf2->nfaces );

  // Check vertices.
  for( int nVertex = 0; nVertex < surf1->nvertices; nVertex++ ) {

    if( surf1->vertices[nVertex].x != surf2->vertices[nVertex].x ) {
      cerr << "Mismatch x on vertex " << nVertex 
	   << ": in " << surf1->vertices[nVertex].x 
	   << ", out " << surf2->vertices[nVertex].x << endl;
      exit( 1 );
    }
    if( surf1->vertices[nVertex].y != surf2->vertices[nVertex].y ) {
      cerr << "Mismatch x on vertex " << nVertex 
	   << ": in " << surf1->vertices[nVertex].y 
	   << ", out " << surf2->vertices[nVertex].y << endl;
      exit( 1 );
    }
    if( surf1->vertices[nVertex].z != surf2->vertices[nVertex].z ) {
      cerr << "Mismatch x on vertex " << nVertex 
	   << ": in " << surf1->vertices[nVertex].z 
	   << ", out " << surf2->vertices[nVertex].z << endl;
      exit( 1 );
    }
    
  }

  // Check faces.
  for( int nFace = 0; nFace < surf1->nfaces; nFace++ ) {

    for( int nVertex = 0; nVertex < VERTICES_PER_FACE; nVertex++ ) {

      if( surf1->faces[nFace].v[nVertex] != surf2->faces[nFace].v[nVertex] ) {
	cerr << "Mismatch face index on face " << nFace 
	     << ", vertex " << nVertex
	     << ": in " << surf1->faces[nFace].v[nVertex] 
	     << ", out " << surf2->faces[nFace].v[nVertex] << endl;
	exit( 1 );
      }

    }
  }
  
  return 0;
}
