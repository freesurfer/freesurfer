#ifndef __vtkODFGlyph_h
#define __vtkODFGlyph_h

#include <vtkStructuredPointsToPolyDataFilter.h>

class vtkLookupTable;
class vtkPolyData;

class vtkODFGlyph : public vtkStructuredPointsToPolyDataFilter {

 public:
  vtkTypeMacro(vtkODFGlyph,vtkStructuredPointsToPolyDataFilter);
  //void PrintSelf(ostream& os, vtkIndent indent);

  // Description
  // Construct object
  static vtkODFGlyph *New();

  // Description:
  // Get/Set factor by which to scale each odf
  vtkSetMacro(ScaleFactor,float);
  vtkGetMacro(ScaleFactor,float);

  // Description:
  // Get/Set factor by which to scale each odf
  vtkSetMacro(ScaleByGFA,int);
  vtkGetMacro(ScaleByGFA,int);
  vtkBooleanMacro(ScaleByGFA,int);

  // Description:
  // Get/Set vtkLookupTable which holds color values for current output
  //vtkSetObjectMacro(ColorTable,vtkLookupTable);
  vtkGetObjectMacro(ColorTable,vtkLookupTable);

  // Description:
  // Get/Set range of GFA input scalars
  vtkGetVector2Macro(GFARange,double);
  vtkSetVector2Macro(GFARange,double);

protected:
  vtkODFGlyph();
  ~vtkODFGlyph();

  void Execute();

  double GFARange[2];

  float ScaleFactor; // Factor by which to scale each odf
  int ScaleByGFA; // scale each odf glyph by sqrt(gfa)

  int BrightnessLevels; // # of sets of NUM_SPHERE_POINTS values in ColorTable. Each set at a different brightness gradation.
  vtkLookupTable *ColorTable; // color table for current output. indeces match
                              // scalars of output's pointdata

private:
  vtkODFGlyph(const vtkODFGlyph&);  // Not implemented.
  void operator=(const vtkODFGlyph&);  // Not implemented.

  static const int ODF_SIZE = 752;
  static const int NUM_SPHERE_POINTS = ODF_SIZE;
  static const double SPHERE_POINTS[NUM_SPHERE_POINTS][3];

};

#endif
