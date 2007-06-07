#include <stdexcept>
#include <sstream>

#include "vtkKWBltGraph.h"
#include "vtkCommand.h"
#include "vtkObjectFactory.h"

using namespace std;

vtkStandardNewMacro( vtkKWBltGraph );
vtkCxxRevisionMacro( vtkKWBltGraph, "$Revision: 1.6 $" );

int const vtkKWBltGraph::MouseoverEnterElementEvent = vtkCommand::UserEvent + 1;
int const vtkKWBltGraph::MouseoverExitElementEvent = vtkCommand::UserEvent + 2;

vtkKWBltGraph::vtkKWBltGraph() :
  XAxisTitle( NULL ),
  YAxisTitle( NULL ),
  mMouseoverDistanceToElement( 10 ),
  mbCurrentlyOverElement( false ) {

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

  if (!vtkKWCoreWidget::CreateSpecificTkWidget(this, 
	"blt::graph", "-plotbackground white -relief raised -border 2")) {
    vtkErrorMacro("Failed creating widget " << this->GetClassName());
    return;
  }

  this->Bind();
}

void
vtkKWBltGraph::Bind () {

  this->SetBinding( "<Motion>", this, "MotionCallback %W %x %y" );
}

void
vtkKWBltGraph::UnBind () {
  
  this->RemoveBinding( "<Motion>" );
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
    XAxisTitle = new char[strlen(isTitle) + 1];
    strncpy( XAxisTitle, isTitle, sizeof(XAxisTitle) );
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
    YAxisTitle = new char[strlen(isTitle) + 1];
    strncpy( YAxisTitle, isTitle, sizeof(YAxisTitle) );
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

  assert( isLabel );
  
  // See if this already exists in our vector. 
  GraphElementList::iterator tElement;
  for( tElement = mElements.begin(); tElement != mElements.end(); ++tElement ){
    if( tElement->msLabel == isLabel ) {
      throw runtime_error( "Element already exists!" );
    }
  }

  if( iPoints.size() == 0 || iPoints.size() % 2 != 0 ) {
    throw runtime_error( "Points must have an even number of elements" );
  }
  
  GraphElement element;
  element.msLabel = isLabel;
  element.mPoints = iPoints;
  if( isSymbol )
    element.msSymbol = isSymbol;
  element.mLineWidth = iLineWidth;
  element.mRed = iRed;
  element.mGreen = iGreen;
  element.mBlue = iBlue;

  mElements.push_back( element );
		
  this->Modified();
}

void
vtkKWBltGraph::DeleteAllElements () {

  GraphElementList::iterator tElement;
  for( tElement = mElements.begin(); tElement != mElements.end(); ++tElement ){
    
    GraphElement& element = *tElement;
    
    this->Script( "%s element delete %s", 
		  this->GetWidgetName(), element.msLabel.c_str() );
  }

  mElements.clear();
}

void
vtkKWBltGraph::UpdateXAxisTitle () {

  if( this->IsCreated() ) {

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

  GraphElementList::iterator tElement;
  for( tElement = mElements.begin(); tElement != mElements.end(); ++tElement ){
    
    GraphElement& element = *tElement;
    
    stringstream ssData;
    vector<double>::iterator tPoint;
    for( tPoint = element.mPoints.begin();
	 tPoint != element.mPoints.end(); ++tPoint ) {
      ssData << (*tPoint) << " ";
    }

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
vtkKWBltGraph::MotionCallback ( const char* isElement, int iX, int iY ) {

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
      double pointX, pointY;
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
	    foundElement.msLabel = strdup( sMinElement.c_str() );
	    foundElement.mnPointInElement = nMinPoint;
	    foundElement.mElementX = minX;
	    foundElement.mElementY = minY;
	    foundElement.mWindowX = iX;
	    foundElement.mWindowY = iY;
	    foundElement.mDistanceToElement = distance;
	    this->InvokeEvent( MouseoverEnterElementEvent, &foundElement );
	    free( foundElement.msLabel );

	    mbCurrentlyOverElement = true;
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

vtkKWBltGraph::GraphElement::GraphElement () :
  msLabel( "" ),
  msSymbol( "" ), mRed( 0 ), mGreen( 0 ), mBlue( 0 ) {
}

