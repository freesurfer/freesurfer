#ifndef __vtkFDTensorGlyph_h
#define __vtkFDTensorGlyph_h

#include <vtkStructuredPointsToPolyDataFilter.h>
#include <vtkTransform.h>

class vtkLookupTable;

class vtkFDTensorGlyph : public vtkStructuredPointsToPolyDataFilter {

 public:
  vtkTypeMacro(vtkFDTensorGlyph,vtkStructuredPointsToPolyDataFilter);
  //void PrintSelf(ostream& os, vtkIndent indent);

  // Description
  // Construct object with scaling on and 
  static vtkFDTensorGlyph *New();

  // Description:
  // Turn on/off glyph scaling with substitution of 4, 1, 1 factors for eigenvalues.
  vtkSetMacro(UniformScaling,int);
  vtkGetMacro(UniformScaling,int);
  vtkBooleanMacro(UniformScaling,int);

  // Description:
  // Turn on/off glyph scaling factor of sqrt of the FA.
  vtkSetMacro(FAScaling,int);
  vtkGetMacro(FAScaling,int);
  vtkBooleanMacro(FAScaling,int);

  // Description:
  // Specify scale factor to scale object by. (Scale factor always affects
  // output even if scaling is off.)
  vtkSetMacro(ScaleFactor,float);
  vtkGetMacro(ScaleFactor,float);

  // Description:
  // Get/Set vtkLookupTable which holds color values for current output
  //vtkSetObjectMacro(ColorTable,vtkLookupTable);
  vtkGetObjectMacro(ColorTable,vtkLookupTable);
  
  // Description:
  // reorients the tensors to be in the right orientation
  vtkSetObjectMacro(VoxelToMeasurementFrameTransform, vtkTransform);
  vtkGetObjectMacro(VoxelToMeasurementFrameTransform, vtkTransform);
  
  // Description:
  // reorients the tensors to be in the right orientation for a particular slice
  vtkSetObjectMacro(VoxelToSliceTransform, vtkTransform);
  vtkGetObjectMacro(VoxelToSliceTransform, vtkTransform);

protected:
  vtkFDTensorGlyph();
  ~vtkFDTensorGlyph();

  void Execute();

  vtkTransform *VoxelToMeasurementFrameTransform;

  vtkTransform *VoxelToSliceTransform;

  int UniformScaling; // Determine whether eigenvalues are used as scaling
                      // factors in corresponding eigenvector directions
  int FAScaling; // Determine whether sqrt(fa) is used as a scaling factor
  float ScaleFactor; // Scale factor to use to scale geometry

  vtkLookupTable *ColorTable; // color table for current output. indeces match
                              // scalars of output's pointdata

  enum {
    FA = 0,
    EV1 = 1,
    EV2 = 4,
    EV3 = 7,
    EVA = 10,
  };

  static const double SMALL_EIGENVALUE = 1e-8;

private:
  vtkFDTensorGlyph(const vtkFDTensorGlyph&);  // Not implemented.
  void operator=(const vtkFDTensorGlyph&);  // Not implemented.

  void ComputeScalarRangeGreaterThanZero(const int component, double range[2]);
  double GetScalarMean(const int component);
};

#endif
