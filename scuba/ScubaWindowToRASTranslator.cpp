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
