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
  // The event constants this object will send.
  //BTX
  static unsigned long const MouseoverEnterElementEvent;
  static unsigned long const MouseoverExitElementEvent;
  //ETX

  // A small class returned by the kEventMouseoverElement event,
  // describing the element that was mouseovered.
  //BTX
  class SelectedElementAndPoint {
  public:
    char* msLabel;               // Label of the element
    int mnPointInElement;        // Index of the data point in the element
    double mElementX, mElementY; // The graph x,y of the element
    int mWindowX, mWindowY;      // The window x,y of the event
    double mDistanceToElement;   // Distance from mouse to event in window
  };
  //ETX

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
  // Show or hide the legend.
  void SetLegendVisible ( int ibVisible );
  void SetLegendVisibleToOn ();
  void SetLegendVisibleToOff ();

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
		    const char* isSymbol, int iLineWidth,
		    double iRed, double iGreen, double iBlue );
  //ETX
  
  // Description:
  // Delete all elements
  void DeleteAllElements ();

  void Draw ();

  void MotionCallback ( const char* isElement, int iX, int iY );
  
 protected:

  vtkKWBltGraph ();
  virtual ~vtkKWBltGraph ();

  // Description:
  void CreateWidget ();
  
  void Bind ();
  void UnBind ();

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

  double mMouseoverDistanceToElement;
  bool mbCurrentlyOverElement;
  
  //BTX
  class GraphElement {
  public:
    GraphElement ();
    std::string msLabel;
    std::vector<double> mPoints; // x, y, x, y, etc
    int mLineWidth;
    std::string msSymbol;
    double mRed, mGreen, mBlue;
  };
  
  // Description:
  typedef std::vector<vtkKWBltGraph::GraphElement> GraphElementList;
  GraphElementList mElements;
  //ETX
};


#endif
