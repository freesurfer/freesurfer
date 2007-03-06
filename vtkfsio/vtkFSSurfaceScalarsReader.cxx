#include <stdexcept>
#include <sstream>

extern "C" {
#include "error.h"
#include "mrisurf.h"
}

#include "vtkFSSurfaceScalarsReader.h"
#include "vtkDataObject.h"
#include "vtkFloatArray.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

using namespace std;

vtkCxxRevisionMacro(vtkFSSurfaceScalarsReader, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkFSSurfaceScalarsReader);

vtkFSSurfaceScalarsReader::vtkFSSurfaceScalarsReader () :
  FileName( "" ),
  NumberOfValues( 0 ) {

  this->SetNumberOfInputPorts( 0 );
}

vtkFSSurfaceScalarsReader::~vtkFSSurfaceScalarsReader () {

}

void
vtkFSSurfaceScalarsReader::SetFileName ( const char* ifn ) {

  FileName = ifn;
}

const char*
vtkFSSurfaceScalarsReader::GetFileName () const {

  return FileName.c_str();
}

int
vtkFSSurfaceScalarsReader::RequestData ( vtkInformation*,
					vtkInformationVector**,
					vtkInformationVector* iOutputVector ){


  // Init a float array.
  vtkSmartPointer<vtkFloatArray> scalars = vtkFloatArray::New();

  // Try to read the scalars.
  float* values = NULL;
  int eRead = MRISreadValuesIntoArray( this->FileName.c_str(),
				       this->NumberOfValues,
				       &values );
  if( ERROR_NONE != eRead ) {
    if( values ) free( values );
    vtkErrorMacro(<< "Could not read scalar file " << this->FileName.c_str() );
    return 0; // 0 is failure
  }
  
  // Allocate our scalars.
  scalars->Allocate( this->NumberOfValues );
  scalars->SetNumberOfComponents( 1 );

  // Copy our array into the scalars.
  for( int nValue = 0; nValue < this->NumberOfValues; nValue ++ )
    scalars->InsertNextValue( values[nValue] );

  // MRISreadValuesIntoArray allocated the array, so we free it.
  free( values );
  
  // Set the scalars in the output.
  vtkPolyData* output = vtkPolyData::GetData( iOutputVector );
  if( !output ) {
    vtkErrorMacro(<< "No output for vtkFSSurfaceScalarsReader" );
    return 0; // 0 is failure
  }
  output->GetPointData()->SetScalars( scalars );

  return 1; // 1 is success
}


void
vtkFSSurfaceScalarsReader::PrintSelf ( ostream& iStream, vtkIndent iIndent ) {
    
  vtkPolyDataAlgorithm::PrintSelf( iStream, iIndent );
}
  
  
