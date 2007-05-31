#ifndef vtkKWBltGraph_h
#define vtkKWBltGraph_h

#include <string>
#include <vector>

#include "vtkKWCoreWidget.h"
#include "vtkKWOptions.h"

class vtkKWBltGraph : public vtkKWCoreWidget {

 public:

  static vtkKWBltGraph* New();
  vtkTypeRevisionMacro( vtkKWBltGraph, vtkKWCoreWidget );

  // Description:
  // Plot background color
  virtual void GetPlotBackgroundColor(double *r, double *g, double *b);
  virtual double* GetPlotBackgroundColor();
  virtual void SetPlotBackgroundColor(double r, double g, double b);
  virtual void SetPlotBackgroundColor(double rgb[3])
    { this->SetPlotBackgroundColor(rgb[0], rgb[1], rgb[2]); };

  // Description:
  // Relief
  virtual void SetRelief(int);
  virtual int GetRelief();
  virtual void SetReliefToRaised();
  virtual void SetReliefToSunken();
  virtual void SetReliefToFlat();
  virtual void SetReliefToRidge();
  virtual void SetReliefToSolid();
  virtual void SetReliefToGroove();

  // Description:
  // Border
  virtual void SetBorderWidth(int);
  virtual int GetBorderWidth();

  // Description:
  // X axis title
  virtual void SetXAxisTitle ( const char* );
  vtkGetStringMacro(XAxisTitle);

  // Description:
  // Y axis title
  virtual void SetYAxisTitle ( const char* );
  vtkGetStringMacro(YAxisTitle);

  // Description:
  // Set default element attributes.
  vtkGetStringMacro(DefaultElementSymbol);
  vtkSetStringMacro(DefaultElementSymbol);
  vtkSetVector3Macro(DefaultElementColor,double);

  // Description:
  // Create an element.
  //BTX
  void AddElement ( const char* isLabel, 
		    std::vector<double>& iPoints,
		    const char* isSymbol, 
		    double iRed, double iGreen, double iBlue );
  //ETX
  
  // Description:
  // Delete all elements
  void DeleteAllElements ();

  void Draw ();

 protected:

  vtkKWBltGraph ();
  virtual ~vtkKWBltGraph ();

  // Description:
  void CreateWidget ();
  
  // Description:
  void UpdateXAxisTitle ();
  void UpdateYAxisTitle ();

  // Description:
  char* XAxisTitle;
  char* YAxisTitle;

  // Description:
  // Default element attributes
  char* DefaultElementSymbol;
  double DefaultElementColor[3];
  
  //BTX
  class GraphElement {
  public:
    GraphElement ();
    std::string msLabel;
    std::vector<double> mPoints; // x, y, x, y, etc
    std::string msSymbol;
    double mRed, mGreen, mBlue;
  };
  
  // Description:
  typedef std::vector<vtkKWBltGraph::GraphElement> GraphElementList;
  GraphElementList mElements;
  //ETX
};


#endif
