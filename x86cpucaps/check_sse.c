/*
 * SSE support detection routine
 * 
 *  imported from MPlayer-0.90pre7 (cpudetect.c)
 * 
 */

#ifndef _WIN32_

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

extern int x86cpucaps_check_sse_supported(int, int);
int _has_sse;

#if defined(_POSIX_SOURCE) && defined(X86_FXSR_MAGIC)
static void sigill_handler_sse( int signal, struct sigcontext sc)
{
   printf( "SIGILL, " );

   /* Both the "xorps %%xmm0,%%xmm0" and "divps %xmm0,%%xmm1"
    * instructions are 3 bytes long.  We must increment the instruction
    * pointer manually to avoid repeated execution of the offending
    * instruction.
    *
    * If the SIGILL is caused by a divide-by-zero when unmasked
    * exceptions aren't supported, the SIMD FPU status and control
    * word will be restored at the end of the test, so we don't need
    * to worry about doing it here.  Besides, we may not be able to...
    */
   sc.eip += 3;

  _has_sse=0;
}

static void sigfpe_handler_sse( int signal, struct sigcontext sc)
{
   printf( "SIGFPE, " );

   if ( sc.fpstate->magic != 0xffff ) {
      /* Our signal context has the extended FPU state, so reset the
       * divide-by-zero exception mask and clear the divide-by-zero
       * exception bit.
       */
      sc.fpstate->mxcsr |= 0x00000200;
      sc.fpstate->mxcsr &= 0xfffffffb;
   } else {
      /* If we ever get here, we're completely hosed.
       */
      printf( "\n\n" );
      printf( "SSE enabling test failed badly!" );
   }
}
#endif /*_POSIX_SOURCE && X86_FXSR_MAGIC */

int x86cpucaps_check_sse_supported(int i, int verbose)
{
#if defined(_POSIX_SOURCE) && defined(X86_FXSR_MAGIC)
   struct sigaction saved_sigill;
   struct sigaction saved_sigfpe;

   _has_sse = i;
   
   /* Save the original signal handlers.
    */
   sigaction( SIGILL, NULL, &saved_sigill );
   sigaction( SIGFPE, NULL, &saved_sigfpe );

   signal( SIGILL, (void (*)(int))sigill_handler_sse );
   signal( SIGFPE, (void (*)(int))sigfpe_handler_sse );

   /* Emulate test for OSFXSR in CR4.  The OS will set this bit if it
    * supports the extended FPU save and restore required for SSE.  If
    * we execute an SSE instruction on a PIII and get a SIGILL, the OS
    * doesn't support Streaming SIMD Exceptions, even if the processor
    * does.
    */
   if ( _has_sse ) {
      if (verbose) printf( "Testing OS support for SSE... " );

//      __asm __volatile ("xorps %%xmm0, %%xmm0");
      __asm __volatile ("xorps %xmm0, %xmm0");

      if (verbose) {
	 if ( _has_sse ) {
	    printf( "yes.\n" );
	 } else {
	    printf( "no!\n" );
	 }
      }
      
   }

   /* Emulate test for OSXMMEXCPT in CR4.  The OS will set this bit if
    * it supports unmasked SIMD FPU exceptions.  If we unmask the
    * exceptions, do a SIMD divide-by-zero and get a SIGILL, the OS
    * doesn't support unmasked SIMD FPU exceptions.  If we get a SIGFPE
    * as expected, we're okay but we need to clean up after it.
    *
    * Are we being too stringent in our requirement that the OS support
    * unmasked exceptions?  Certain RedHat 2.2 kernels enable SSE by
    * setting CR4.OSFXSR but don't support unmasked exceptions.  Win98
    * doesn't even support them.  We at least know the user-space SSE
    * support is good in kernels that do support unmasked exceptions,
    * and therefore to be safe I'm going to leave this test in here.
    */
    if ( _has_sse ) {
      if (verbose) printf( "Testing OS support for SSE unmasked exceptions... " );

//      test_os_katmai_exception_support();

       if (verbose) {
	  if ( _has_sse ) {
	     printf( "yes.\n" );
	  } else {
	     printf( "no!\n" );
	  }
       }
       
   }

   /* Restore the original signal handlers.
    */
   sigaction( SIGILL, &saved_sigill, NULL );
   sigaction( SIGFPE, &saved_sigfpe, NULL );

   /* If we've gotten to here and the XMM CPUID bit is still set, we're
    * safe to go ahead and hook out the SSE code throughout Mesa.
    */
   if (verbose) {
      if ( _has_sse ) {
	 printf( "Tests of OS support for SSE passed.\n" );
      } else {
	 printf( "Tests of OS support for SSE failed!\n" );
      }
   }
   
#else
   /* We can't use POSIX signal handling to test the availability of
    * SSE, so we disable it by default.
    */
   if (verbose) printf( "Cannot test OS support for SSE, disabling to be safe.\n" );
   _has_sse =0;
#endif /* _POSIX_SOURCE && X86_FXSR_MAGIC */

   return _has_sse;
}

#endif /* _WIN32_ */
