#include <stdlib.h>
#include "xDebug.h"

tBoolean xDbg_gbOutput   = FALSE;
FILE*    xDbg_gStream    = NULL;
int      xDbg_gType      = xDebug_Nothing;
char*    xDbg_gsRequest  = NULL;


void xDbg_Init () {

  /* if env variable is defined, check it and set debug type if something is
     requested. do file if nothing is recognizable. otherwise try to do 
     file. */
  xDbg_gsRequest = getenv( "XDEBUG" );                   
  if( NULL == xDbg_gsRequest ) {                         
    xDbg_gType = xDebug_File;                            
  } else {                                              
    if( strcmp( "file", xDbg_gsRequest ) == 0 )         
      xDbg_gType = xDebug_File;                          
    else if( strcmp( "stderr", xDbg_gsRequest ) == 0 )   
      xDbg_gType = xDebug_Print;                         
    else if( strcmp( "none", xDbg_gsRequest ) == 0 )     
      xDbg_gType = xDebug_Nothing;                       
    else                                                 
      xDbg_gType = xDebug_Print;                         
  }                                                      
  if( xDebug_File == xDbg_gType ) {                      
    xDbg_gStream = fopen( ".tkmedit_xdebug", "w" );      
    if( NULL == xDbg_gStream )                           
      xDbg_gType= xDebug_Nothing;                        
  }                                                      
  if( xDebug_Print == xDbg_gType )                       
    xDbg_gStream = stderr;                              
  if( xDebug_Nothing == xDbg_gType ) {                 
    xDbg_gbOutput = FALSE;                            
  } else {                                             
    xDbg_gbOutput = TRUE;                               
  }                                    
}

void xDbg_ShutDown () {

  /* close file if we opened it */
  if( xDebug_File == xDbg_gType                     
      && NULL != xDbg_gStream )                     
    fclose( xDbg_gStream );
}


void xDbg_PrintStatus () {

  fprintf( stderr, "output = %d\n", (int)xDbg_gbOutput );
  fprintf( stderr, "type = %s\n",
     (xDbg_gType == xDebug_Nothing) ? "nothing" :
     (xDbg_gType == xDebug_Print) ? "print" :
     (xDbg_gType == xDebug_File) ? "file" : "" );
  fprintf( stderr, "env var = %s\n", 
     (xDbg_gsRequest != NULL) ? xDbg_gsRequest : "undefined" );
  if( xDbg_gStream == stderr ) 
    fprintf( stderr, "stream = stderr\n" );
  else if( NULL != xDbg_gStream )
    fprintf( stderr, "stream = probably a file\n" );
  else
    fprintf( stderr, "stream = NULL\n" );
}
