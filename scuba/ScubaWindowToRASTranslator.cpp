#include "ScubaWindowToRASTranslator.h"

void 
ScubaWindowToRASTranslator::TranslateWindowToRAS( int iXWindow, int iYWindow,
						  float& oXRAS, float& oYRAS,
						  float& oZRAS ) {
  oXRAS = oYRAS = oZRAS = 0;
}

void 
ScubaWindowToRASTranslator::TranslateWindowToRAS( int iXWindow, int iYWindow, 
						  float oRAS[3] ) {
  TranslateWindowToRAS( iXWindow, iYWindow, oRAS[0], oRAS[1], oRAS[2] );
}

void 
ScubaWindowToRASTranslator::TranslateRASToWindow( float iXRAS, float iYRAS,
						  float iZRAS,
					  int& oXWindow, int& oYWindow ) {
  oXWindow = oYWindow = 0;
}
				    
void 
ScubaWindowToRASTranslator::TranslateRASToWindow( float iRAS[3], 
						  int& oXWindow, int& oYWindow ) {
  TranslateRASToWindow( iRAS[0], iRAS[1], iRAS[2], oXWindow, oYWindow );
}
