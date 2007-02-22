#ifndef COMMANDUPDATE_H_
#define COMMANDUPDATE_H_

#include <itkCommand.h>
#include <itkDiffusionTensor3D.h>

#include "dmri_poistats/datamodel/utils/itkPoistatsFilter.h"

class CommandUpdate : public itk::Command
{

public:
  typedef CommandUpdate Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
  
  typedef itk::DiffusionTensor3D< float > TensorPixelType;
  typedef itk::Image< TensorPixelType, 3 > TensorImageType;
  typedef itk::Image< float, 3 > OutputImageType;
  typedef itk::PoistatsFilter< TensorImageType, OutputImageType > 
    PoistatsFilterType;
  typedef const PoistatsFilterType* PoistatsFilterPointer;
//  typedef const PoistatsFilterType* PoistatsFilterPointer;

  /** This is an itk callback */
  void Execute( itk::Object *caller, const itk::EventObject & event );
  
  void Execute(const itk::Object * object, const itk::EventObject & event);
  
  void PostMessage( const std::string message );

  void PostErrorMessage( const std::string message );
  
  static void WriteMessage( const std::string message, 
    const std::string fileName );

  itkGetMacro( OutputDirectory, std::string );
  itkSetMacro( OutputDirectory, std::string );

  itkGetMacro( LogFileName, std::string );
  itkSetMacro( LogFileName, std::string );
  
protected:
  CommandUpdate() {};

private:
  std::string m_OutputDirectory;
  
  std::string m_LogFileName;
  
};

#endif /*COMMANDUPDATE_H_*/
