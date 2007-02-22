#include <iostream>
#include <fstream>
#include <string>

// I need to include this for ofstream to work for some reason...
#include <itkImageSeriesReader.h>

#include "CommandUpdate.h"

#include "PoistatsEvents.h"

void
CommandUpdate::Execute( itk::Object *caller, const itk::EventObject & event ){
  Execute( (const itk::Object *) caller, event);
}

void 
CommandUpdate::Execute(const itk::Object * object, 
  const itk::EventObject & event){

  if( itk::IterationEvent().CheckEvent( &event ) ) {
    
    PoistatsFilterPointer filter =
      dynamic_cast< PoistatsFilterPointer >( object );
      
    const int currentIteration = filter->GetCurrentIteration();
    const int currentLull = filter->GetCurrentLull();
    const double energyDifference = filter->GetCurrentEnergyDifference();
    const double meanOfEnergies = filter->GetCurrentMeanOfEnergies();
    const double minOfEnergies = filter->GetCurrentMinOfEnergies();
    const double globalMin = filter->GetGlobalMinEnergy();
    const int exchanges = filter->GetExchanges();
    const double elapsedTime = filter->GetElapsedTime();
    
    std::ostringstream output;
    output << currentIteration << "   lull: " << currentLull << "  denergy: " << energyDifference << "  mean: " << meanOfEnergies << "  min: " << minOfEnergies << "  global min: " <<  globalMin << "  exchs: " << exchanges << "  time: " << elapsedTime << "\n";
    
    PostMessage( output.str() );
    
  } else if( itk::PoistatsOdfCalculationEvent().CheckEvent( &event ) ) {

    if( itk::PoistatsOdfCalculationStartEvent().CheckEvent( &event ) ) {
      PostMessage( "constructing odf list\n" );
    } else if( itk::PoistatsOdfCalculationProgressEvent().CheckEvent( &event ) ){
      PostMessage( ". " );
    } else if( itk::PoistatsOdfCalculationEndEvent().CheckEvent( &event ) ) {
      PoistatsFilterPointer filter = 
        dynamic_cast< PoistatsFilterPointer >( object );

      std::ostringstream output;

      output << "\nfinished\n" << "elapsed time: " << filter->GetElapsedTime() << std::endl;
          
      PostMessage( output.str() );
    }

  } else if( itk::SeedsEvent().CheckEvent( &event ) ) {
    
    if( itk::SeedsUsingAllEvent().CheckEvent( &event ) ) {
      
      PostMessage( "Seeds to use not explicitly set, using all seeds in seed volume...\n" );      

    } else if( itk::SeedsFoundInitialEvent().CheckEvent( &event ) ) {

      PoistatsFilterPointer filter = 
        dynamic_cast< PoistatsFilterPointer >( object );
      const int nInitialPoints = filter->GetNumberOfInitialPoints();
      std::ostringstream output;
      output << "Finding " << nInitialPoints << " seeds" << std::endl;
      PostMessage( output.str() );

    }
    
  } else if( itk::GenerateOutputEvent().CheckEvent( &event ) ) {
    
    if( itk::GenerateOptimalPathDensitiesEvent().CheckEvent( &event ) ) {
      
      PostMessage( "calculating optimal path densities\n" );
      
    } else if( itk::GenerateBestReplicaPathDensitiesEvent().
      CheckEvent( &event ) ) {

      if( itk::GenerateBestReplicaPathDensitiesStartEvent().
        CheckEvent( &event ) ) {
          
        PostMessage( "calculating best replica path densities\n" );
        
      } else if( itk::GenerateBestReplicaPathDensitiesProgressEvent().
        CheckEvent( &event ) ) {
          
        PostMessage( ". " );

      } else if( itk::GenerateBestReplicaPathDensitiesEndEvent().
        CheckEvent( &event ) ) {
        PostMessage( "\nfinished\n" );
      }
      
    } else if( itk::GenerateFinalPathProbabilitiesEvent().
      CheckEvent( &event ) ) {
        
      PostMessage( "calculate final path probabilities\n" );
      
    }
    
  }
  
}

void 
CommandUpdate::PostMessage( const std::string message ){
  
  std::cout << message;
  
  std::string fileName = m_OutputDirectory + "/" + m_LogFileName;
  WriteMessage( message, fileName );
  
}

void 
CommandUpdate::PostErrorMessage( const std::string message ){
  
  std::cerr << "ERROR: " + message;
  
  std::string fileName = m_OutputDirectory + "/" + m_LogFileName;
  WriteMessage( message, fileName );
  
}

void 
CommandUpdate::WriteMessage( const std::string message, 
  const std::string fileName ) {

  std::ofstream logFile( fileName.c_str(), std::ios::out | std::ios::app  );
  
  logFile << message;
  
  logFile.close();
    
}
