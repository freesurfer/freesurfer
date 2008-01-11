#include <stdexcept>
#include <sstream>
#include <assert.h>
#include "vtkKWBltGraph.h"
#include "vtkKWMenu.h"
#include "vtkCommand.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"

using namespace std;

vtkStandardNewMacro( vtkKWBltGraph );
vtkCxxRevisionMacro( vtkKWBltGraph, "$Revision: 1.11 $" );

unsigned long const vtkKWBltGraph::MouseoverEnterElementEvent = 
  vtkCommand::UserEvent + 1;
unsigned long const vtkKWBltGraph::MouseoverExitElementEvent = 
  vtkCommand::UserEvent + 2;
unsigned long const vtkKWBltGraph::ContextualMenuOpening = 
  vtkCommand::UserEvent + 3;

vtkKWBltGraph::vtkKWBltGraph() :
  XAxisTitle( NULL ),
  YAxisTitle( NULL ),
  mMouseoverDistanceToElement( 10 ),
  mbCurrentlyOverElement( false ),
  msCurrentMouseoverElement( "" ) {

  DefaultElementSymbol = new char[10];
  strncpy( DefaultElementSymbol, "plus", 9 );

  DefaultElementColor[0] = 1.0;
  DefaultElementColor[1] = 0.0;
  DefaultElementColor[2] = 0.0;
}

vtkKWBltGraph::~vtkKWBltGraph () {

  if( XAxisTitle )
    delete [] XAxisTitle;
  if( YAxisTitle )
    delete [] YAxisTitle;
  if( DefaultElementSymbol )
    delete [] DefaultElementSymbol;
}

void
vtkKWBltGraph::CreateWidget () {

  // Create our blt widget.
  if (!vtkKWCoreWidget::CreateSpecificTkWidget(this, 
	"blt::graph", "-plotbackground white -relief raised -border 2")) {
    vtkErrorMacro("Failed creating widget " << this->GetClassName());
    return;
  }

  // Assign our bindings.
  this->Bind();
}

void
vtkKWBltGraph::Bind () {

  // Assign our bindings.
  this->SetBinding( "<Motion>", this, "MotionCallback %W %x %y" );
  this->SetBinding( "<ButtonPress-3>", this, "Button3DownCallback %W %X %Y" );
}

void
vtkKWBltGraph::UnBind () {

  // Remove bindings.
  this->RemoveBinding( "<Motion>" );
  this->RemoveBinding( "<ButtonPress-3>" );
}

void
vtkKWBltGraph::GetPlotBackgroundColor ( double* iRed, double* iGreen,
					double* iBlue ) {

  this->GetConfigurationOptionAsColor ( "-plotbackground", 
					iRed, iGreen, iBlue );
}

double*
vtkKWBltGraph::GetPlotBackgroundColor () {

  return this->GetConfigurationOptionAsColor ( "-pltbackground" );
}

void
vtkKWBltGraph::SetPlotBackgroundColor ( double iRed, double iGreen, 
					double iBlue ) {

  this->SetConfigurationOptionAsColor ( "-plotbackground", 
					iRed, iGreen, iBlue );
}

void 
vtkKWBltGraph::SetRelief ( int iRelief ) {

  this->SetConfigurationOption( "-relief", 
	     vtkKWOptions::GetReliefAsTkOptionValue(iRelief) );
}

void 
vtkKWBltGraph::SetReliefToRaised () {

  this->SetRelief( vtkKWOptions::ReliefRaised );
}

void 
vtkKWBltGraph::SetReliefToSunken () {

  this->SetRelief( vtkKWOptions::ReliefSunken );
}

void 
vtkKWBltGraph::SetReliefToFlat () {

  this->SetRelief( vtkKWOptions::ReliefFlat );
}

void 
vtkKWBltGraph::SetReliefToRidge () {

  this->SetRelief( vtkKWOptions::ReliefRidge );
}

void 
vtkKWBltGraph::SetReliefToSolid () {

  this->SetRelief( vtkKWOptions::ReliefSolid );
}

void 
vtkKWBltGraph::SetReliefToGroove () {

  this->SetRelief( vtkKWOptions::ReliefGroove );
}

int 
vtkKWBltGraph::GetRelief () {

  return vtkKWOptions::GetReliefFromTkOptionValue
    ( this->GetConfigurationOption("-relief") );
}

void
vtkKWBltGraph::SetBorderWidth ( int iWidth ) {

  this->SetConfigurationOptionAsInt( "-border", iWidth );
}

int
vtkKWBltGraph::GetBorderWidth () {

  return this->GetConfigurationOptionAsInt( "-border" );
}

void
vtkKWBltGraph::SetLegendVisible ( int ibVisible ) {

  if( this->IsCreated() ) {
    this->Script( "%s legend configure -hide %d",
		  this->GetWidgetName(), !ibVisible );
  }
}

void
vtkKWBltGraph::SetLegendVisibleToOn () {

  this->SetLegendVisible( 1 );
}

void
vtkKWBltGraph::SetLegendVisibleToOff () {

  this->SetLegendVisible( 0 );
}

void
vtkKWBltGraph::SetXAxisTitle ( const char* isTitle ) {

  if( NULL == XAxisTitle && NULL == isTitle )
    return;

  if( NULL != XAxisTitle && NULL != isTitle &&
      0 == strcmp( XAxisTitle, isTitle ) )
    return;

  if( NULL != XAxisTitle )
    delete [] XAxisTitle;
  
  if( isTitle ) {
    int len = strlen(isTitle);
    XAxisTitle = new char[len + 1];
    strcpy( XAxisTitle, isTitle );
  } else {
    XAxisTitle = NULL;
  }

  this->Modified();
  
  this->UpdateXAxisTitle();
}

void
vtkKWBltGraph::SetYAxisTitle ( const char* isTitle ) {

  if( NULL == YAxisTitle && NULL == isTitle )
    return;

  if( NULL != YAxisTitle && NULL != isTitle &&
      0 == strcmp( YAxisTitle, isTitle ) )
    return;

  if( NULL != YAxisTitle )
    delete [] YAxisTitle;
  
  if( isTitle ) {
    int len = strlen(isTitle);
    YAxisTitle = new char[len + 1];
    strcpy( YAxisTitle, isTitle );
  } else {
    YAxisTitle = NULL;
  }

  this->Modified();
  
  this->UpdateYAxisTitle();
}

void
vtkKWBltGraph::AddElement ( const char* isLabel, vector<double>& iPoints,
			    const char* isSymbol, int iLineWidth,
			    double iRed, double iGreen, double iBlue ) {

  if( NULL == isLabel )
    throw runtime_error( "isLabel argument must not be provided" );
  
  // See if we already have an element with this label
  GraphElementList::iterator tElement;
  for( tElement = mElements.begin(); tElement != mElements.end(); ++tElement ){
    if( tElement->msLabel == isLabel ) {
      throw runtime_error( "Element already exists!" );
    }
  }

  // Make sure we have the right number of points.
  if( iPoints.size() == 0 || iPoints.size() % 2 != 0 )
    throw runtime_error( "Points must have an even number of elements" );
  
  // Make a graph element structure with this data.
  GraphElement element;
  element.msLabel = isLabel;
  element.mPoints = iPoints;
  if( isSymbol )
    element.msSymbol = isSymbol;
  element.mLineWidth = iLineWidth;
  element.mRed = iRed;
  element.mGreen = iGreen;
  element.mBlue = iBlue;

  // Add it to our list of elements.
  mElements.push_back( element );
		
  this->Modified();
}

void
vtkKWBltGraph::DeleteAllElements () {

  // Delete all our elements.
  GraphElementList::iterator tElement;
  for( tElement = mElements.begin(); tElement != mElements.end(); ++tElement ){
    
    GraphElement& element = *tElement;

    // This will erase it in the graph.    
    this->Script( "%s element delete %s", 
		  this->GetWidgetName(), element.msLabel.c_str() );
  }

  mElements.clear();
}

void
vtkKWBltGraph::UpdateXAxisTitle () {

  if( this->IsCreated() ) {

    // Set the title in the graph.
    if( NULL != XAxisTitle ) {
      this->Script( "%s axis configure x -title %s",
		    this->GetWidgetName(), XAxisTitle );
    } else {
      this->Script( "%s axis configure x -title {}",
		    this->GetWidgetName() );
    }

  }
}

void
vtkKWBltGraph::UpdateYAxisTitle () {

  if( this->IsCreated() ) {

    // Set the title in the graph.
    if( NULL != YAxisTitle ) {
      this->Script( "%s axis configure y -title %s",
		    this->GetWidgetName(), YAxisTitle );
    } else {
      this->Script( "%s axis configure y -title {}",
		    this->GetWidgetName() );
    }

  }
}

void
vtkKWBltGraph::Draw () {

  if( !this->IsCreated() )
    return;

  if( mElements.size() == 0 ) 
    return;

  // For each element...
  GraphElementList::iterator tElement;
  for( tElement = mElements.begin(); tElement != mElements.end(); ++tElement ){
    
    GraphElement& element = *tElement;

    // Make a string with the point data.
    stringstream ssData;
    vector<double>::iterator tPoint;
    for( tPoint = element.mPoints.begin();
	 tPoint != element.mPoints.end(); ++tPoint ) {
      ssData << (*tPoint) << " ";
    }

    // Use a Tcl command to set the element data in the BLT graph widget.
    this->Script( "%s element create %s "
		  "-data {%s} " 
		  "-color #%02x%02x%02x " // format to hex
		  "-linewidth %d "
		  "-symbol %s",
		  this->GetWidgetName(),
		  element.msLabel.c_str(),
		  ssData.str().c_str(),
		  (int)(element.mRed * 255.0),
		  (int)(element.mGreen * 255.0),
		  (int)(element.mBlue * 255.0),
		  element.mLineWidth,
		  element.msSymbol.c_str() );

  }
}

void
vtkKWBltGraph::MotionCallback ( const char* isWidget, int iX, int iY ) {

  // We're getting window coords, try to convert to graph coords.
  const char* sReply = 
    this->Script( "%s invtransform %d %d", this->GetWidgetName(), iX, iY );
  if( sReply ) {

    // Extract the x and y from the reply.
    float x, y;
    stringstream ssReply( sReply );
    ssReply >> x >> y;

    // Go through our elements and try to find the closest point.
    GraphElementList::iterator tElement;
    for( tElement = mElements.begin(); 
	 tElement != mElements.end(); ++tElement ){
    
      GraphElement& element = *tElement;

      double minDist2 = numeric_limits<double>::max();
      string sMinElement = "";
      int nMinPoint = 0;
      double minX = 0, minY = 0;
      
      vector<double>::iterator tPoint;
      bool bX = true;
      double pointX = 0, pointY = 0;
      int nPoint = 0;
      for( tPoint = element.mPoints.begin();
	   tPoint != element.mPoints.end(); ++tPoint ) {
	// Get an x one round, then get a y the next one.
	if( bX ) {
	  pointX = (*tPoint);
	  bX = false;
	} else {
	  pointY = (*tPoint);
	  // Calc the distance and compare.
	  double dist2 = (pointX-x)*(pointX-x) + (pointY-y)*(pointY-y);
	  if( dist2 < minDist2 ) {
	    // If this is a new min, save its info.
	    minDist2 = dist2;
	    sMinElement = element.msLabel;
	    nMinPoint = nPoint;
	    minX = pointX;
	    minY = pointY;
	  }
	  bX = true;
	  nPoint++;
	}
      }

      if( minDist2 < numeric_limits<double>::max() ){
	// Convert our min x and y back to window coords and compare
	// the distance to our mouse over distance. We convert to
	// window coords so we can do the comparison in pixels.
	sReply = this->Script( "%s transform %f %f",
			       this->GetWidgetName(), minX, minY );
	if( sReply ) {
	  float minWindowX, minWindowY;
	  stringstream ssReply( sReply );
	  ssReply >> minWindowX >> minWindowY;
	  double distance = sqrt( (minWindowX-iX)*(minWindowX-iX) +
				  (minWindowY-iY)*(minWindowY-iY) );
	  if( distance < mMouseoverDistanceToElement ) {
	    // We found a good one. Notify our observers.
	    SelectedElementAndPoint foundElement;
	    foundElement.msLabel = sMinElement.c_str();
	    foundElement.mnPointInElement = nMinPoint;
	    foundElement.mElementX = minX;
	    foundElement.mElementY = minY;
	    foundElement.mWindowX = iX;
	    foundElement.mWindowY = iY;
	    foundElement.mDistanceToElement = distance;
	    this->InvokeEvent( MouseoverEnterElementEvent, &foundElement );
	    
	    // Remember information about this moused over element.
	    mbCurrentlyOverElement = true;
	    msCurrentMouseoverElement = sMinElement.c_str();
	    return;
	  }
	}
      }
    }
  }

  // If we got here, we couldn't find an element to mouseover. If we
  // were previously over one, send an event that we no longer are.
  if( mbCurrentlyOverElement ) {
    
    this->InvokeEvent( MouseoverExitElementEvent, NULL );
    mbCurrentlyOverElement = false;
  }
}

void
vtkKWBltGraph::Button3DownCallback ( const char* isWidget, int iX, int iY ) {

  // If we're moused over an element...
  if( mbCurrentlyOverElement ) {
    
    // Make a menu.
    vtkSmartPointer<vtkKWMenu> popup = 
      vtkSmartPointer<vtkKWMenu>::New();
    popup->SetParent( this );
    popup->Create();
    
    // Make a struct with this menu and the currently moused over
    // element, and invoke the event that will let observers populate
    // the menu.
    ContextualMenuElement data;
    data.mMenu = popup;
    data.msElement = msCurrentMouseoverElement.c_str();
    this->InvokeEvent( ContextualMenuOpening, (void*)&data );
    
    // If we have any items in this menu, pop it up.
    if( popup->GetNumberOfItems() > 0 )
      popup->PopUp( iX, iY );
  }

}

vtkKWBltGraph::GraphElement::GraphElement () :
  msLabel( "" ),
  msSymbol( "" ), mRed( 0 ), mGreen( 0 ), mBlue( 0 ) {
}

