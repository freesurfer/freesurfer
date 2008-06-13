#ifndef vtkKWBltGraph_h
#define vtkKWBltGraph_h

#include <string>
#include <vector>

#include "vtkKWCoreWidget.h"
#include "vtkKWOptions.h"

class vtkKWMenu;

class vtkKWBltGraph : public vtkKWCoreWidget {

 public:

  static vtkKWBltGraph* New();
  vtkTypeRevisionMacro( vtkKWBltGraph, vtkKWCoreWidget );

  // Description:
  // The event constants this object will send.
  //BTX
  static unsigned long const MouseoverEnterElementEvent;
  static unsigned long const MouseoverExitElementEvent;
  static unsigned long const ContextualMenuOpening;
  //ETX

  // A small struct returned by the EventMouseoverEnterElement event,
  // describing the element that was mouseovered.
  //BTX
  struct SelectedElementAndPoint {
    char const* msLabel;         // Label of the element
    int mnPointInElement;        // Index of the data point in the element
    double mElementX, mElementY; // The graph x,y of the element
    int mWindowX, mWindowY;      // The window x,y of the event
    double mDistanceToElement;   // Distance from mouse to event in window
  };
  //ETX

  // A small struct returned by the ContextualMenuOpening event,
  // describing the element of which the menu is in the context.
  //BTX
  struct ContextualMenuElement {
    vtkKWMenu* mMenu; 		 // The menu to populate
    char const* msElement;       // The element
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
  // Draws all the elements that have been set.
  void Draw ();

  // Description:
  // Delete <all> element<s>
  void DeleteElement ( const char* isLabel );
  void DeleteAllElements ();

  // Description:
  // These callbacks are for Tk events in the widget.
  void MotionCallback ( const char* isWidget, int iX, int iY );
  void Button3DownCallback ( const char* isWidget, int iX, int iY );

  void SavePostscript ( const char* ifnPS );

 protected:

  vtkKWBltGraph ();
  virtual ~vtkKWBltGraph ();

  // Description:
  // Make our BLT graph widget.
  void CreateWidget ();

  // Create and delete or Tk bindings so we can respond to Tk events.
  void Bind ();
  void UnBind ();

  // Description:
  void UpdateXAxisTitle ();
  void UpdateYAxisTitle ();

  //BTX
  // Description:
  // String titles for our axes.
  char* XAxisTitle;
  char* YAxisTitle;

  // Description:
  // Default element attributes
  char* DefaultElementSymbol;
  double DefaultElementColor[3];

  // Description:
  // Information about the currently moused-over element. This is set
  // when an element is moused over.
  double mMouseoverDistanceToElement;
  bool mbCurrentlyOverElement;
  std::string msCurrentMouseoverElement; // only valid if 
                                         // mbCurrentlyOverElement is true

  // Description:
  // Internal representation of an element.
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
  // Our elements.
  typedef std::vector<vtkKWBltGraph::GraphElement> GraphElementList;
  GraphElementList mElements;
  //ETX
};


#endif
