#ifndef ScubaWindowToRASTranslator_h
#define ScubaWindowToRASTranslator_h

class ScubaWindowToRASTranslator {
  
 public: 
  virtual void TranslateWindowToRAS( int iXWindow, int iYWindow,
				     float& oXRAS, float& oYRAS,
				     float& oZRAS );
  
  virtual void TranslateWindowToRAS( int iXWindow, int iYWindow,
				     float oRAS[3] );
};


#endif
