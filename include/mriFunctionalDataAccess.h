
#ifndef tkmFunctionalDataAccess_H
#define tkmFunctionalDataAccess_H

#include "xTypes.h"
#include "xVoxel.h"
#include "xList.h"
#include "mriTransform.h"

                                   /* error constants. */
typedef enum {
  FunD_tErr_NoError = 0,
  FunD_tErr_InvalidPtr,
  FunD_tErr_PathNotFound,
  FunD_tErr_CouldntGuessStem,
  FunD_tErr_DataNotFound,
  FunD_tErr_HeaderNotFound,
  FunD_tErr_CouldntAllocateVolume,
  FunD_tErr_UnrecognizedHeaderFormat,
  FunD_tErr_QuestionableHeaderFormat,
  FunD_tErr_CouldntDetermineDataType,
  FunD_tErr_CouldntAllocateStorage,
  FunD_tErr_DataAlreadyRead,
  FunD_tErr_SliceFileNotFound,       
  FunD_tErr_ErrorReadingSliceData,
  FunD_tErr_CouldntAllocateMatrix,
  FunD_tErr_CouldntReadRegisterFile,
  FunD_tErr_CouldntCalculateDeviations,
  FunD_tErr_ErrorAccessingTransform,
  FunD_tErr_InvalidParameter,
  FunD_tErr_InvalidTimeResolution,
  FunD_tErr_InvalidNumPreStimTimePoints,
  FunD_tErr_InvalidFunctionalVoxel,
  FunD_tErr_InvalidConditionIndex,
  FunD_tErr_InvalidTimePoint,
  FunD_tErr_InvalidErrorCode,
  FunD_tErr_knNumErrorCodes
} FunD_tErr;

typedef enum {
  FunD_tDataType_Short,
  FunD_tDataType_Float
} FunD_tDataType;

                                   /* the volume struct. */
typedef struct {

  // location of data sources.
  char mPath[256],                     // the path to the data files
    mStem[256],                        // the stem of all data files
    mHeaderStem[256],                  // stem of header
    mRegistration[256],                // location of registration (path)
    mSubjectName[256];                 // the subjects name

  // the dimensions of the volume.
  FunD_tDataType mDataType;
  int mNumRows, mNumCols, mNumSlices;  // volume dimensions
  int mNumTimePoints;
  int mNumConditions;                  // number of non-null conditions. note
                                       // that 0 is a non-null condition.

  // information about the functional volume.
  float mPixelSize, mSliceThickness, mBrightnessScale;
  float mMaxValue, mMinValue;

  // optional values that may not be used in every volume.
  float mTimeResolution;               // num seconds between time points.
  int mNumPreStimTimePoints;           // num of time points before stimuleus
  int mNumDataValues;                  // num of data points

  // for quickly calcaluating indecies, these are the number of values in 
  // each plane of data.
  int mDataPlaneSizes[5];

  // transformation objects.
  mriTransformRef mIdxToIdxTransform;
  mriTransformRef mRASToIdxTransform;

  // slice condition plane row col

  // data in [condition][time][slice_k][row_j][col_i]
  // data isnt a char, could be floats or shorts.
  char * mData;

  // when we know we have error values present, these are used
  char mIsErrorDataPresent;
  float * mDeviations;  // deviations in [condition][timepoint]
  float * mCovMtx;      // Ch in [condtions*timepts][conditions*timepts]
  float * mSigma;         // used in error calc

  // temp storage in conversions
 xVoxelRef mpAnaRAS;
 xVoxelRef mpFuncRAS;

} mriFunctionalData, *mriFunctionalDataRef;


                                   /* allocates and destroys volumes. new
              volumes require a pathname. allocater
              parses header files and initializes all
              data except for reading in the volume
              data. */

FunD_tErr FunD_New    ( mriFunctionalDataRef* outVolume,
      char*                 inPathName, 
            char*                 inStem,
            char*                 inHeaderStem,
            char*                 inRegistrationPath  );
FunD_tErr FunD_Delete ( mriFunctionalDataRef * ioVolume );

                                   /* looks in directory and guesses a stem. */
FunD_tErr FunD_GuessStem ( char*inPathName, char* outStem );

                                   /* different methods of parsing header
              files. */
FunD_tErr FunD_ParseStemHeader    ( mriFunctionalDataRef this );
FunD_tErr FunD_ParseAnalyseHeader ( mriFunctionalDataRef this );
FunD_tErr FunD_ParseBFileHeader   ( mriFunctionalDataRef this );

                                   /* parsing support for reading keywords
              and values. */
FunD_tErr FunD_ReadKeywordAndValue ( FILE* inFile, 
             char* inExpectedKeyword, 
             char* inValueType, 
             int   inNumValues, 
             int   inValueSize, 
             char* inAssign );

                                   /* reads the register.dat file, allocates
              and initializes matricies */
FunD_tErr FunD_ParseRegistrationAndInitMatricies ( mriFunctionalDataRef this );

                                   /* saves the registration to file, making
              a backup if it already exists */
FunD_tErr FunD_SaveRegistration ( mriFunctionalDataRef this );

                                   /* sets registration to identity matrix */
FunD_tErr FunD_SetRegistrationToIdentity ( mriFunctionalDataRef this );


                                   /* looks for data files to determine the
              data type.*/
FunD_tErr FunD_DetermineDataType ( mriFunctionalDataRef this );

                                   /* reads the data, processing it according
              to the dimensions and values in the
              structure. for example, if our flag
              signalling error values are present is
1              set, it will read every other plane as
              a sigma value. */
FunD_tErr FunD_ParseData ( mriFunctionalDataRef this );

                                   /* calcs the deviations. needs sigma and
              CovMtx, should be done after reading
              the header _and_ data files.*/
FunD_tErr FunD_CalcDeviations ( mriFunctionalDataRef this );

/* applies a transformation to the registration, which can later be
   saved out as a new registration */
FunD_tErr FunD_ApplyTransformToRegistration ( mriFunctionalDataRef this,
                MATRIX*             iTransform );
FunD_tErr FunD_TranslateRegistration        ( mriFunctionalDataRef this,
                float                ifDistance,
                tAxis                iAxis );
FunD_tErr FunD_RotateRegistration           ( mriFunctionalDataRef this,
                float                ifDegrees,
                tAxis                iAxis,
                xVoxelRef   iCenterFuncRAS );
FunD_tErr FunD_ScaleRegistration            ( mriFunctionalDataRef this,
                float                ifFactor,
                tAxis                iAxis );

FunD_tErr FunD_GetDataAtAnaIdx ( mriFunctionalDataRef this,
        xVoxelRef             inVoxel, 
         int                  inPlane, 
         int                  inCondition,
         float*               outData );
FunD_tErr FunD_GetDataAtRAS    ( mriFunctionalDataRef this,
        xVoxelRef             inVoxel, 
         int                  inPlane, 
         int                  inCondition,
         float*               outData );

FunD_tErr FunD_GetDataAtAnaIdxForAllTimePoints ( mriFunctionalDataRef this,
            xVoxelRef             inVoxel, 
             int              inCondition, 
             float*              outData );

                                   /* deviations are calced by multplying
              sigma by the square root of the diagonal
              of the CovMtx, then indexing in by time
              major (time,cond 0,0, 1,0, 2,0...) */
FunD_tErr FunD_GetDeviation ( mriFunctionalDataRef this,
            int                  inCondition, 
            int                  inTimePoint,
            float*               outValue );

FunD_tErr FunD_GetDeviationForAllTimePoints ( mriFunctionalDataRef this, 
                int                 inCondition, 
                float*               outData );

                                   /* converts a time point index to a second
              based on the time resolution and time
              resolution */
FunD_tErr FunD_ConvertTimePointToSecond ( mriFunctionalDataRef this,
            int                  inTimePoint,
            float*               outSecond );
FunD_tErr FunD_ConvertSecondToTimePoint ( mriFunctionalDataRef this,
            float                inSecond,
            int*                 outTimePoint );

                                   /* setting these values changes the way
              time points are converted to seconds. */
FunD_tErr FunD_SetTimeResolution       ( mriFunctionalDataRef this, 
           float                inTimeRes );
FunD_tErr FunD_SetNumPreStimTimePoints ( mriFunctionalDataRef this,
           int                  inNumPoints );

                                   /* accessors for volume information. */
FunD_tErr FunD_GetPath                 ( mriFunctionalDataRef this, 
           char*                out );
FunD_tErr FunD_GetStem                 ( mriFunctionalDataRef this, 
           char*                out );
FunD_tErr FunD_GetSubjectName          ( mriFunctionalDataRef this,
           char*                out );
FunD_tErr FunD_GetNumRows              ( mriFunctionalDataRef this, 
           int*                 out );
FunD_tErr FunD_GetNumCols              ( mriFunctionalDataRef this, 
           int*                 out );
FunD_tErr FunD_GetNumSlices            ( mriFunctionalDataRef this, 
           int*                 out );
FunD_tErr FunD_GetNumTimePoints        ( mriFunctionalDataRef this, 
           int*                 out );
FunD_tErr FunD_GetNumConditions        ( mriFunctionalDataRef this, 
           int*                 out );
FunD_tErr FunD_GetTimeResolution       ( mriFunctionalDataRef this, 
           float*               out );
FunD_tErr FunD_GetNumPreStimTimePoints ( mriFunctionalDataRef this, 
           int*                 out );
FunD_tErr FunD_GetNumDataValues        ( mriFunctionalDataRef this, 
           int*                 out );
FunD_tErr FunD_GetValueRange           ( mriFunctionalDataRef this,
           float*               outMin,
           float*               outMax );
                                   /* gets bounds in anatomical index coords.
              this is actually a bounding box in 
              anatomical space. if the functional
              data is oblique to the anatomical,
              it will have coords that arn't valid
              in functional space, but it is a good
              place to start iterating for the 
              in anatomical space. */
FunD_tErr FunD_GetBoundsInAnatomical ( mriFunctionalDataRef this, 
              xVoxelRef             outBeginCorner,
              xVoxelRef             outEndCorner );

FunD_tErr FunD_DebugPrint ( mriFunctionalDataRef this );

char* FunD_GetErrorString ( FunD_tErr inErr );

                                   /* internal functions */

                                   /* true if we can calc and return
              deviations. */
tBoolean FunD_IsErrorDataPresent ( mriFunctionalDataRef this );

                                   /* convertxs between anatomical and 
              functional voxel coordinates. */
void FunD_ConvertAnaIdxToFuncIdx ( mriFunctionalDataRef this,
          xVoxelRef             inAnatomicalIdx,
          xVoxelRef             outFunctionalIdx );
void FunD_ConvertFuncIdxToAnaIdx ( mriFunctionalDataRef this,
          xVoxelRef             inFunctionalIdx,
          xVoxelRef             outAnatomicalIdx );

                                   /* convert between ras and func idx
              coords */
void FunD_ConvertRASToFuncIdx    ( mriFunctionalDataRef this,
          xVoxelRef             inRAS,
          xVoxelRef             outFunctionalIdx );
void FunD_ConvertFuncIdxToFuncRAS ( mriFunctionalDataRef this,
            xVoxelRef            iFuncIdx,
            xVoxelRef            oFuncRAS );

                                   /* convert between functional voxel 
              coordinates and an index into the
              storage array. */
inline int FunD_CoordsToIndex  ( mriFunctionalDataRef this,
         xVoxelRef            inFunctionalVoxel,
         int                  inConditionIndex, 
         int                  inTimePoint );
inline void FunD_IndexToCoords ( mriFunctionalDataRef this,
         int                  inIndex,
         xVoxelRef            outFunctionVoxel,
         int*                 outConditionIndex, 
         int*                 outTimePoint );
             

                                   /* makes slice file names out of 
              pathnames and stems, and tries to find
              the slice files. stops when it can't
              find the next and notes that as the
              number of slices. */
int FunD_CountSliceFiles ( mriFunctionalDataRef this );

                                   /* builds a slice file name */
void FunD_MakeSliceFileName       ( mriFunctionalDataRef this, 
            int                  inSliceNumber,
            char*                inFileNameToSet );
void FunD_MakeBShortSliceFileName ( char*                inPath, 
            char*                inStem, 
            int                  inSliceNumber, 
            char*                ioFileName );
void FunD_MakeBFloatSliceFileName ( char*                inPath, 
            char*                inStem, 
            int                  inSliceNumber, 
            char*                ioFileName );

                                   /* counts how many average data values are
              in a plane so we can do quick 
              index -> coord conversions with mod and 
              division. */
void FunD_CalcDataPlaneSizes ( mriFunctionalDataRef this );

                                   /* basic access functions, using coordinates
              in functional space. assumes all index
              values are valid. */
inline float FunD_GetValue ( mriFunctionalDataRef this,
          xVoxelRef             inFunctionalVoxel,
           int                  inConditionIndex, 
           int                  inTimePoint );
inline void FunD_SetValue  ( mriFunctionalDataRef this,
          xVoxelRef             inFunctionalVoxel,
           int                  inConditionIndex, 
           int                  inTimePoint, 
           float                inValue);

                                   /* performs basic sanity checks. */
FunD_tErr  FunD_AssertIsValid ( mriFunctionalDataRef this );

                                   /* bounds checking */
tBoolean FunD_IsTimeResolutionValid       ( mriFunctionalDataRef this, 
              float                inTimeRes );
tBoolean FunD_IsNumPreStimTimePointsValid ( mriFunctionalDataRef this, 
              int                inNumPoints );
tBoolean FunD_IsFunctionalVoxelValid      ( mriFunctionalDataRef this,
              xVoxelRef             inVoxel );
tBoolean FunD_IsConditionIndexValid       ( mriFunctionalDataRef this, 
              int            inConditionIndex );
tBoolean FunD_IsTimePointValid            ( mriFunctionalDataRef this, 
              int                  inTimePoint );
tBoolean FunD_IsIndexValid                ( mriFunctionalDataRef this, 
              int                  inIndex );



#endif
