#ifndef ScubaWindowToRASTranslator_h
#define ScubaWindowToRASTranslator_h

class ScubaWindowToRASTranslator {
  
 public: 
  virtual void TranslateWindowToRAS( int iWindow[2], float oRAS[3] );
  virtual void TranslateRASToWindow( float iRAS[3], int oWindow[2] );
};


#endif
