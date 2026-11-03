#ifndef __kvlAtlasParameterEstimationConsole_h
#define __kvlAtlasParameterEstimationConsole_h

#include "kvlAtlasParameterEstimationConsoleGUI.h"
#include <vector>
#include <string>
#include "kvlAtlasParameterEstimator.h"



namespace kvl
{

class AtlasParameterEstimationConsole : public kvlAtlasParameterEstimationConsoleGUI
{

public: 

  //
  AtlasParameterEstimationConsole();
  
  //
  virtual ~AtlasParameterEstimationConsole();

  //
  void SetLabelImages( const std::vector< std::string >&  fileNames );
    
  // 
  void Show();

  //
  void Estimate();
  
  //
  void HandleEstimatorEvent( itk::Object* object, const itk::EventObject & event );
  
protected:

  //
  void DisplayLabelImage( unsigned int labelImageNumber );
  
  //
  void InitializeMesh();
  
  //
  void Step();

  //
  void Interrupt();

  //
  void Continue();

  //
  void SetPositionEstimationResolution( unsigned int positionEstimationResolution );

private:  
  
  AtlasParameterEstimator::Pointer  m_Estimator;
  
  std::string  m_TotalProgressLabel;
  std::string  m_SubProgressLabel;
  
  bool  m_Interrupted;
  bool  m_Stepping;
  
  
};



} // end namespace kvl

#endif
