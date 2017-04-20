#if !defined STDAFX_H_
# define STDAFX_H_

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2006 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
/*
 * stdafx.h : include file for standard system include files,
 *  or project specific include files that are used frequently, but
 *      are changed infrequently
 * Used for Windows only
 */
#if defined WIN32
  /*  Windows builds are UNICODE */
  #pragma warning(disable : 4786) /* truncated identifiers in debug symbol table. */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */


  #if defined TECPLOTKERNEL 
    #include "custafx.h"

    #if defined _WIN64
      #if !defined _M_IA64 && !defined _M_AMD64
        #error "UNKNOWN 64-BIT PLATFORM"
      #endif
    #endif

    /* TECPLOT preprocessor defines */
    #if !defined MSWIN
      #define MSWIN
    #endif
    #define THREED
    #if defined _WIN64 
      #define NO_ELM
    #endif
    #if !defined NO_RESTRICT
      #define RESTRICT
      #if !defined NO_ELM
        #define ELM
      #endif
    #endif

    /* SmartHeap must be first so its libraries are scanned first for malloc/free/etc. */
    #ifdef USE_SMARTHEAP
      #define USE_SMARTHEAP_VERSION_8
      #if defined _WIN64
        #error "SmartHeap not available for Win64"
      #else
        #ifdef _DEBUG
          #error "Don't build debug versions with SmartHeap"
        #else
          #ifdef USE_SMARTHEAP_VERSION_8
            #include "\\cedar\workgroups\development\tecplot\Builds\smartheap\8.0.0\include\smrtheap.h"
            #pragma comment(lib, "\\\\cedar\\workgroups\\development\\tecplot\\Builds\\smartheap\\8.0.0\\msvc\\shdw32mt.lib")
          #else
            #include "\\cedar\workgroups\development\tecplot\Builds\smartheap\6.0.3\include\smrtheap.h"
            #pragma comment(lib, "\\\\cedar\\workgroups\\development\\tecplot\\Builds\\smartheap\\6.0.3\\msvc\\shdw32mt.lib")
          #endif
        #endif
      #endif
    #endif /* USE_SMARTHEAP */

    #if defined NDEBUG && !defined CHECKED_BUILD
    #  define NO_ASSERTS
    #endif

    #define VC_EXTRALEAN		/* Exclude rarely-used stuff from Windows headers */

    #include <afxwin.h>         /* MFC core and standard components */
    #include <afxext.h>         /* MFC extensions */
    /* #include <afxdao.h> */
    #include <afxinet.h>
    #include <afxmt.h> /* multi-threaded stuff */
    #include <afxtempl.h> /* templated collection classes */

    #ifndef _AFX_NO_AFXCMN_SUPPORT
    #include <afxcmn.h>			/* MFC support for Windows Common Controls */
    #endif /* _AFX_NO_AFXCMN_SUPPORT */

    #ifndef _AFX
    #error MFC is not defined as far as Tecplot is concerned.
    #endif

    #include <cderr.h> /* comm dialog error codes */
    #include <shlobj.h>
    #include <winver.h>
    #include <mbstring.h>

    #include "MASTER.h"

    /*
     * These are necessary, otherwise you get massive warnings. CAM 03/17/2004
     */
    #pragma warning (disable : 4244) /* conversion: 1856 warnings */
    #pragma warning (disable : 4100) /* unreferenced formal parameter: 331 warnings */

    /*
     * This one doesn't appear to do anything in debug, but supresses
     * countless warnings in release builds. CAM 03/05/2004
     */
    #pragma warning (disable : 4711) /* inline function */

    /*
     * It would like to turn this one back on but it created 218 warnings, more
     * that I had time to deal with. CAM 03/05/2004
     */
    #pragma warning (disable : 4701) /* variable 'may' be used without having been initialized */  

    /* linker settings now in stdafx.cpp */

  #else /* !TECPLOTKERNEL */
    #define AfxIsValidAddress(ptr,bb) ((ptr)!=NULL)
  #endif

  /* 64-bit adjustments */
  #if defined _M_IA64 || defined _M_AMD64
    #define WININT  INT_PTR
    #define WINUINT UINT_PTR
  #else
    #define WININT  int
    #define WINUINT UINT
  #endif
  #define WINCALLBACK CALLBACK
  
  #if defined TECPLOTKERNEL
    /* message tracing */
    #if defined (NDEBUG)
      #define TRACE_WINDOWS_MESSAGE(m,w,l)
    #else
      extern void InternalWindowsTraceMessage(UINT m, WPARAM w, LPARAM l);
      #define TRACE_WINDOWS_MESSAGE(m,w,l) InternalWindowsTraceMessage(m,w,l)
    #endif
  #endif /* TECPLOTKERNEL */
#endif /* WIN32 */

/************************************************************
 * Platform Independant
 ************************************************************/

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

#if !defined TECPLOTKERNEL
  /* For other than Tecplot, TranlsatedString_pa is always a const char * */
  typedef const char *TranslatedString_pa; 
#endif /* !defined TECPLOTKERNEL */


#endif /* STDAFX_H_ */
