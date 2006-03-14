#ifndef ScubaWindowToRASTranslator_h
#define ScubaWindowToRASTranslator_h

class ScubaWindowToRASTranslator {
  
 public: 
  virtual ~ScubaWindowToRASTranslator () {};
  virtual void TranslateWindowToRAS( int const iWindow[2], float oRAS[3] );
  virtual void TranslateRASToWindow( float const iRAS[3], int oWindow[2] );
};


#endif
