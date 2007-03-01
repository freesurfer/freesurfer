#include "vtkFDTensorGlyph.h"

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>

vtkStandardNewMacro(vtkFDTensorGlyph);

// Construct object with scaling 
vtkFDTensorGlyph::vtkFDTensorGlyph() {
  this->UniformScaling = 0;
  this->FAScaling = 0;
  this->ScaleFactor = 2;
  this->ColorTable = vtkLookupTable::New();
  this->VoxelToMeasurementFrameTransform = vtkTransform::New();
}

vtkFDTensorGlyph::~vtkFDTensorGlyph() {
  this->ColorTable->Delete();
  this->VoxelToMeasurementFrameTransform->Delete();
}

void vtkFDTensorGlyph::Execute() {
  vtkImageData *input = static_cast<vtkImageData *>(this->GetInput());
  vtkPolyData *output = this->GetOutput();

  double inputSpacing[3];
  input->GetSpacing(inputSpacing);
  int *dim = input->GetDimensions();
  int *extent = input->GetExtent();
  int indexArray[3];
  int maxGlyphs = dim[0] * dim[1] * dim[2];
  float *tensorComponents, pt[3];
  double translate[3];
  double eigvec1[3], eigvec2[3], eigvec3[3];
  
  float cubePoints[8][3] = { 
    {0.5, -0.5, 0.5},
    {0.5, 0.5, 0.5},
    {-0.5, 0.5, 0.5},
    {-0.5, -0.5, 0.5},
    {0.5, -0.5, -0.5},
    {0.5, 0.5, -0.5},
    {-0.5, 0.5, -0.5},
    {-0.5, -0.5, -0.5}
  };
  
  int cubePolys[6][4] = {
    {0, 1, 2, 3},
    {0, 1, 5, 4},
    {4, 5, 6, 7},
    {6, 7, 3, 2},
    {2, 1, 5, 6},
    {3, 0, 4, 7}
  };

  float uniformScaleFactors[3] = {1, .25, .25};
  float faval;
  double normalizedFA;
  double normalizedEigenValues[3];

  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  vtkTransform *transform = vtkTransform::New();
  vtkPoints *points = vtkPoints::New();

  vtkFloatArray *pointsArray = vtkFloatArray::New();
  pointsArray->SetNumberOfComponents(3);
  pointsArray->SetNumberOfTuples(maxGlyphs * 8);

  vtkCellArray *polys = vtkCellArray::New();
  polys->Allocate(polys->EstimateSize(maxGlyphs * 6, 4));

  vtkIntArray *colorScalars = vtkIntArray::New();
  colorScalars->SetNumberOfValues(maxGlyphs * 8);

  vtkUnsignedCharArray *colors = vtkUnsignedCharArray::New();
  double wxyz[4];
  static_cast<vtkTransform *>(this->VoxelToMeasurementFrameTransform->GetInverse())->GetOrientationWXYZ(wxyz);
  vtkTransform *measurementToVTKRotation = vtkTransform::New();
  measurementToVTKRotation->RotateWXYZ(wxyz[0], wxyz[1], wxyz[2], wxyz[3]);

  colors->SetNumberOfComponents(4);
  colors->SetNumberOfTuples(maxGlyphs);
  
  double faRange[2];
  double eigenValue1Range[2];

  //would improve performance with the expense of code redundancy if
  //computed both ranges in a single iteration
  this->ComputeScalarRangeGreaterThanZero(FA, faRange);
  this->ComputeScalarRangeGreaterThanZero(EVA, eigenValue1Range);

  int count = 0, colorCount = 0;
  for (int x = 0; x < dim[0]; x++) {
    for (int y = 0; y < dim[1]; y++) {
      for (int z = 0; z < dim[2]; z++) {
	tensorComponents =
	  (float *)input->GetScalarPointer(x + extent[0],
					   y + extent[2],
					   z + extent[4]);
	// Should this be used instead?
	//	tensorComponents = (float *)input->GetScalarPointer(x, y, z);
	faval = tensorComponents[FA];
	if (faval == 0) continue;
	normalizedFA = (faval - faRange[0]) / (faRange[1] - faRange[0]);

	for (int i = 0; i < 3; i++) {
	  normalizedEigenValues[i] = (tensorComponents[EVA + i] > 0) ? tensorComponents[EVA + i] : this->SMALL_EIGENVALUE;
	  normalizedEigenValues[i] = tensorComponents[EVA + i] / eigenValue1Range[1];
	}

	indexArray[0] = x;
	indexArray[1] = y;
	indexArray[2] = z;
	input->GetPoint(input->ComputePointId(indexArray), translate);

	for (int i = 0; i < 3; i++) {
	  eigvec1[i] = tensorComponents[EV1 + i];
	  eigvec2[i] = tensorComponents[EV2 + i];
	  eigvec3[i] = tensorComponents[EV3 + i];
	}

	
	measurementToVTKRotation->TransformPoint(eigvec1, eigvec1);
	measurementToVTKRotation->TransformPoint(eigvec2, eigvec2);
	measurementToVTKRotation->TransformPoint(eigvec3, eigvec3);

// 	vtkMath::Normalize(eigvec1);
// 	vtkMath::Normalize(eigvec2);
// 	vtkMath::Normalize(eigvec3);

// 	indexArray[0] = x + extent[0];
// 	indexArray[1] = y + extent[2];
// 	indexArray[2] = z + extent[4];
// 	this->VoxelToMeasurementFrameTransform->TransformPoint(indexArray, translate);
	//cout << "image x: " << x << ", y: " << y << ", z: " << z << endl;
	//cout << "world x: " << translate[0] << ", y: " << translate[1] << ", z: " << translate[2] << endl;

	//this should no longer work:
	//for incorrect dtreg output, uncomment these and swap y and z in matrix
// 	translate[1] *= -1;
// 	translate[2] *= -1;

	for (int matrixRow = 0; matrixRow < 3; matrixRow++) {
	  matrix->SetElement(matrixRow, 0, eigvec1[matrixRow] * (this->UniformScaling ? uniformScaleFactors[0] : normalizedEigenValues[0]) * (this->FAScaling ? normalizedFA : 1) * inputSpacing[matrixRow] * this->ScaleFactor);
	  matrix->SetElement(matrixRow, 1, eigvec2[matrixRow] * (this->UniformScaling ? uniformScaleFactors[1] : normalizedEigenValues[1]) * (this->FAScaling ? normalizedFA : 1) * inputSpacing[matrixRow] * this->ScaleFactor);
	  matrix->SetElement(matrixRow, 2, eigvec3[matrixRow] * (this->UniformScaling ? uniformScaleFactors[2] : normalizedEigenValues[2]) * (this->FAScaling ? normalizedFA : 1) * inputSpacing[matrixRow] * this->ScaleFactor);
	  matrix->SetElement(matrixRow, 3, translate[matrixRow]);
	}

	transform->SetMatrix(matrix);
	
	for (int i = 0; i < 8; i++) {
	  transform->TransformPoint(cubePoints[i], pt);
	  pointsArray->SetTuple(count + i, pt);
	  colorScalars->SetValue(count + i, colorCount);
	}

	colors->SetTuple4(colorCount, sqrt(faval) * fabs(eigvec1[0]) * 255.0, sqrt(faval) * fabs(eigvec1[1]) * 255.0, sqrt(faval) * fabs(eigvec1[2]) * 255.0 + .5, 255);
	colorCount++;
	
	
	for (int i = 0; i < 6; i++) {
	  polys->InsertNextCell(4);
	  for (int j = 0; j < 4; j++) {
	    polys->InsertCellPoint(cubePolys[i][j] + count);
	  }
	}
	count += 8;
      }
    }
  }

  pointsArray->Resize(count);
  colorScalars->Resize(count);
  polys->Squeeze();

  points->SetData(pointsArray);

  output->SetPoints(points);
  output->SetPolys(polys);
  output->GetPointData()->SetScalars(colorScalars);

  //output->Squeeze();
  if (colorCount > 0) {
    colors->Resize(colorCount);
    this->ColorTable->SetNumberOfTableValues(colorCount);
    this->ColorTable->SetTableRange(0, colorCount);
//     unsigned char *outPtr = this->ColorTable->GetPointer(0);
//     //cout << "out: " << outPtr[1500 * 4] << endl;
//     unsigned char *inPtr = colors->GetPointer(0);
//     //cout << "in: " << inPtr[1500 * 4] << endl;
//     memcpy(outPtr, inPtr, colorCount * 4 - 4);
//     //cout << "new out: " << outPtr[1500 * 4] << endl;
    for (int i = 0; i < colorCount; i++) {
      double *tuple = colors->GetTuple(i);
      this->ColorTable->SetTableValue(i, tuple[0] / 255.0, tuple[1] / 255.0, tuple[2] / 255.0, tuple[3] / 255.0);
    }
//     cout << "colors: " << colors->GetComponent(1500, 0) << endl;
//     cout << "colortable: " << ((int *)(this->ColorTable->MapValue(1500)))[1] << endl;
//     cout << "colortable: " << ((int *)(this->ColorTable->MapValue(1501)))[1] << endl;
//     cout << "colortable: " << ((int *)(this->ColorTable->MapValue(1502)))[1] << endl;
//     cout << "colortable: " << ((int *)(this->ColorTable->MapValue(1503)))[1] << endl;
//     cout << "colortable: " << ((int *)(this->ColorTable->MapValue(1504)))[1] << endl;
//     cout << "table # vals: " << this->ColorTable->GetNumberOfTableValues() << endl;
//     cout << "table range: " << this->ColorTable->GetRange()[1] << endl;
  }
  
  matrix->Delete();
  transform->Delete();
  points->Delete();
  pointsArray->Delete();
  polys->Delete();
  colorScalars->Delete();
  colors->Delete();
}

void vtkFDTensorGlyph::ComputeScalarRangeGreaterThanZero(const int component, double range[2]) {
  vtkImageData *input = this->GetInput();
  float *scalarPointer = static_cast<float *>(input->GetScalarPointer());
  int numberOfComponents = input->GetNumberOfScalarComponents();
  if (component >= numberOfComponents) {
    vtkErrorMacro("The component index requested is out of bounds of the number of components: " << component <<  " >= " << numberOfComponents);
    return;
  }
  int numberOfTuples = input->GetNumberOfPoints();
  int numberOfScalars = numberOfComponents * numberOfTuples;

  range[0] = VTK_DOUBLE_MAX;
  range[1] = VTK_DOUBLE_MIN;

  for (int tupleIndex = 0; tupleIndex < numberOfScalars; tupleIndex += numberOfComponents) {
    double currentValue = scalarPointer[tupleIndex + component];
    if (currentValue > 0) {
      if (currentValue < range[0]) {
	range[0] = currentValue;
      } else if (currentValue > range[1]) {
	range[1] = currentValue;
      }
    }
  } //end iteration of tupleIndex
}


// void vtkFDTensorGlyph::PrintSelf(ostream& os, vtkIndent indent) {
// }
