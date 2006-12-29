/**
 * @file  Layer.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.35 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef Layer_h
#define Layer_h

#include "string_fixed.h"
#ifdef HAVE_APPLE_OPENGL_FRAMEWORK
#  include "OpenGL/gl.h"
#else
#  include "GL/gl.h"
#endif
#include <map>
#include "ViewState.h"
#include "IDTracker.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "ScubaWindowToRASTranslator.h"
#include "ScubaToolState.h"
#include "InputState.h"
#include "Listener.h"
#include "Broadcaster.h"
#include "DataCollection.h"

class LayerStaticTclListener : public DebugReporter,
      public TclCommandListener {

public :
  LayerStaticTclListener ();
  ~LayerStaticTclListener ();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );
};


class Layer : public DebugReporter,
      public IDTracker<Layer>,
      public TclCommandListener,
      public Listener,    // dataChanged
      public Broadcaster  // layerChanged, layerInfoSettingsChanged
{

  friend class ScubaViewTester;

public:

  Layer();
  virtual ~Layer();

  static bool kbDefaultReportInfo; // Initial value for whether the
  // layer reports info

  void SetBytesPerPixel ( int icBytes ) {
    mBytesPerPixel = icBytes;
  }

  // Tell the layer to draw its contents into a GL frame buffer.
  virtual void DrawIntoBuffer( GLubyte* iBuffer, int iWidth, int iHeight,
                               ViewState& iViewState,
                               ScubaWindowToRASTranslator& iTranslator );

  // Tell the layer to draw its contents that need openGl commands.
  virtual void DrawIntoGL ( ViewState& iViewState,
                            ScubaWindowToRASTranslator& iTranslator );

  // This is what layers will fill out when GetInfoAtRAS is called.
  class InfoAtRAS {
  public:
    InfoAtRAS ();
    void SetLabel ( std::string is ) {
      msLabel = is;
    }
    void SetValue ( std::string is ) {
      msValue = is;
    }
    void SetTclCallback ( std::string is ) {
      msTclCallback = is;
    }
    void SetInputFilter ( std::string is ) {
      msInputFilter = is;
    }
    void SetShortenHint ( bool ib ) {
      mbShortenHint = ib;
    }
    void Clear();
    std::string GetLabel () {
      return msLabel;
    }
    std::string GetValue () {
      return msValue;
    }
    std::string GetTclCallback () {
      return msTclCallback;
    }
    std::string GetInputFilter () {
      return msInputFilter;
    }
    bool GetShortenHint () {
      return mbShortenHint;
    }
  protected:
    std::string msLabel, msValue;  // label/value to display
    std::string msTclCallback;     // Function to call on input ("" if none)
    // The callback call will be $msTclCallback $mID $input
    std::string msInputFilter;     // Filter to use when calling callback
    bool mbShortenHint;     // Whether client should shorten value
  };

  // Asks the layer to describe a point of data by making InfoAtRAS
  // structs.
  virtual void GetInfoAtRAS ( float iRAS[3],
                              std::list<InfoAtRAS>& ioInfo );

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() {
    return "BaseLayer";
  }

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // Handle broadcast messages.
  virtual void
  DoListenToMessage ( std::string isMessage, void* iData );

  // Get the primary data collection type for this layer. Can be null.
  virtual DataCollection* GetMainDataCollection();

  // Get the selected ROI for this layer. Can be -1. By default, this
  // tries to get the MainDataCollection, and calls GetSelectedROI on
  // that. It returns -1 if there is no MainDataCollection.
  virtual int GetSelectedROI ();


  // Called when the layer's data is changed. Default behavior is to
  // request a redisplay.
  virtual void DataChanged ();

  void SetLabel( std::string isLabel );
  std::string GetLabel();

  virtual void SetOpacity( float iOpacity ) {
    mOpacity = iOpacity;
  }
  float GetOpacity() {
    return mOpacity;
  }

  void SetVisible( bool ibVisible ) {
    mbVisible = ibVisible;
  }
  bool GetVisible() {
    return mbVisible;
  }

  void SetWidth( int iWidth );
  void SetHeight( int iHeight );

  void SetReportInfo( bool ibReport );
  bool GetReportInfo() {
    return mbReportInfoAtRAS;
  }

  virtual void HandleTool ( float iRAS[3], ViewState& iViewState,
                            ScubaWindowToRASTranslator& iTranslator,
                            ScubaToolState& iTool, InputState& iInput );

  // Timer function.
  void Timer ();

  // Redisplay posters.
  void RequestRedisplay() {
    mbPostRedisplay = true;
  }
  bool WantRedisplay() const {
    return mbPostRedisplay;
  }
  void RedisplayPosted() {
    mbPostRedisplay = false;
  }

  // Some drawing tools.
  void DrawPixelIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                             int iWindow[2], int iColor[3], float iOpacity );

  void DrawPixelIntoBuffer ( GLubyte* iAddress,
                             int iColor[3], float iOpacity );

  void DrawAALineIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                              int iFromWindow[2], int iToWindow[2],
                              int iColor[3], int iThickness, float iOpacity );
  void DrawLineIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                            int iFromWindow[2], int iToWindow[2],
                            int iColor[3], int iThickness, float iOpacity );

  // This is the distance to advance 'through' a projection, e.g. for
  // going to the 'next' slice, in RAS coords.
  virtual void GetPreferredThroughPlaneIncrements ( float oIncrements[3] );

  // This is the increment for a slider controlling the brush size.
  virtual float GetPreferredBrushRadiusIncrement ();

  // This is the increment for a slider controlling an unbounded value
  // amount.
  virtual float GetPreferredValueIncrement ();

  // Process a list of options. Just parses them into individual
  // options and hands them to overrideable function.
  // In format option[=value][:option[=value]]...
  void ProcessOptionsList ( std::string isOptionList );

  // Handles an invidual option. Subclass should also call this to
  // handle generic options.
  virtual void ProcessOption ( std::string isOption, std::string isValue );

protected:

  // Overridable timer behavior.
  virtual void DoTimer ();

  int mWidth;
  int mHeight;

  std::string msLabel;

  float mOpacity;

  bool mbVisible;

  static LayerStaticTclListener mStaticListener;

  // Redisplay requested flag.
  bool mbPostRedisplay;

  // Whether to return info in GetInfoAtRAS.
  bool mbReportInfoAtRAS;

  int mBytesPerPixel;
};



#endif
