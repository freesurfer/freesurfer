#include <stdexcept>
#include <sstream>

#include "vtkObjectFactory.h"
#include "vtkKWBltGraph.h"

using namespace std;

vtkStandardNewMacro( vtkKWBltGraph );
vtkCxxRevisionMacro( vtkKWBltGraph, "$Revision: 1.3 $" );

vtkKWBltGraph::vtkKWBltGraph() :
  XAxisTitle( NULL ),
  YAxisTitle( NULL ) {

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
		  //"-color %f %f %f "
		  "-linewidth %d "
		  "-symbol %s",
		  this->GetWidgetName(),
		  element.msLabel.c_str(),
		  ssData.str().c_str(),
		  // element.mRed, element.mGreen, element.mBlue, 
		  element.mLineWidth,
		  element.msSymbol.c_str() );

  }
}

vtkKWBltGraph::GraphElement::GraphElement () :
  msLabel( "" ),
  msSymbol( "" ), mRed( 0 ), mGreen( 0 ), mBlue( 0 ) {
}
