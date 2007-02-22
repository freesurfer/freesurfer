#ifndef POISTATSEVENTS_H_
#define POISTATSEVENTS_H_

#include <itkEventObject.h>

namespace itk
{

/** Store all poistats events here, so that the CommandUpdate can catch them */

itkEventMacro( PoistatsOdfCalculationEvent, ProgressEvent );
itkEventMacro( PoistatsOdfCalculationStartEvent, PoistatsOdfCalculationEvent );
itkEventMacro( PoistatsOdfCalculationProgressEvent, PoistatsOdfCalculationEvent );
itkEventMacro( PoistatsOdfCalculationEndEvent, PoistatsOdfCalculationEvent );

itkEventMacro( SeedsEvent, ProgressEvent );
itkEventMacro( SeedsUsingAllEvent, SeedsEvent );
itkEventMacro( SeedsFoundInitialEvent, SeedsEvent );

itkEventMacro( GenerateOutputEvent, ProgressEvent );
itkEventMacro( GenerateOptimalPathDensitiesEvent, GenerateOutputEvent );
itkEventMacro( GenerateBestReplicaPathDensitiesEvent, GenerateOutputEvent );
itkEventMacro( GenerateBestReplicaPathDensitiesStartEvent, 
               GenerateBestReplicaPathDensitiesEvent );
itkEventMacro( GenerateBestReplicaPathDensitiesProgressEvent, 
               GenerateBestReplicaPathDensitiesEvent );
itkEventMacro( GenerateBestReplicaPathDensitiesEndEvent, 
               GenerateBestReplicaPathDensitiesEvent );
itkEventMacro( GenerateFinalPathProbabilitiesEvent, GenerateOutputEvent );

}

#endif /*POISTATSEVENTS_H_*/
