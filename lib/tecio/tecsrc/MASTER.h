/*****************************************************************
 *****************************************************************
 *******                                                  ********
 ****** (C) Copyright 1988-2005 by Tecplot, Inc.           *******
 *******                                                  ********
 *****************************************************************
 *****************************************************************/
/* BEGINREMOVEFROMADDON */
/* NOTE: All code contained between comments that look like
 *             BEGINREMOVEFROMADDON
 *             ENDREMOVEFROMADDON
 * are pulled out to create the MASTER.h file used in addons.
 */
/* ENDREMOVEFROMADDON */

#ifndef _MASTER_H_
#define _MASTER_H_


/* BEGINREMOVEFROMADDON */
#include "stdafx.h"
/*
 * The following is a temporary fix for figuring out which product is
 * running.  In the future when Focus and 360 use the same code base,
 * we will have to do this dynamically (either with flags on the compiler
 * or variables within Tecplot).
 */
#define TECPLOT_360
/* ENDREMOVEFROMADDON */

#ifdef NDEBUG
# ifdef _DEBUG
#   error "Both NDEBUG and _DEBUG defined"
# endif
#else
# ifndef _DEBUG
#   define _DEBUG
# endif
#endif

/* Now a requirement */
#define USE_3D_HARDWARE

#ifndef THREED
#  define THREED
#endif

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#if defined QUICKDEMO
#define DEMO
#endif

#if defined ELM5
  #define ELM
#endif

#if defined NETRESTRICT
  #define RESTRICT
#endif

#if defined MicrosoftC
#define DOS
#endif

#if defined CRAYX
#define CRAY
#endif

#if defined IRISX
#define IRIS
#endif

#if defined HPX
#define HPUX
#define HP
#endif

#if defined IBMRS6000X
#define IBMRS6000
#endif

#if defined COMPAQALPHAX
#define COMPAQALPHA
#define COMPAQX
#define COMPAQ
#endif

#if defined DECALPHAX
#define DECALPHA
#define DECX
#endif

#if defined DECX
#define DEC
#endif

#if defined SUNSOLARISX || defined SUNSOLARIS86X
#define SUNX
#endif

#if defined SUNX
#define SUN
#endif

#if defined IRISX || defined CRAYX || defined HPX || defined SUNX || defined CONVEXX
#define UNIXX
#define SYSV
#endif

#if defined DECX || defined LINUX || defined IBMRS6000X || defined COMPAQX || defined DARWIN
#define UNIXX
#endif

/* BEGINREMOVEFROMADDON */
#include <stdarg.h>


/* A bit of OEM stuff */
#define OEM_INVALID_CHECKSUM (LgIndex_t) -1

/* Hide the name of the checksum function */
#if defined NDEBUG
# define DECRYPTTIMEDCODE          FixupPlot
# define CHECKHASHEDCODE           ExpandPlot
# define UPDATECLASSICOEMEHCKSUM   ToggleQuadrants
# define UPDATEOEMCHECKSUM         ComputeAngleFromQuatrant
# define InitOemSettings           InitAngleQuatrantSettings
#endif

#if defined MSWIN
#define USE_TRUETYPEFONTS
#endif
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
#ifdef TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#ifdef __cplusplus // STL
#ifdef MSWIN
#endif //MSWIN
#ifdef MSWIN
#endif //MSWIN
#endif //__cplusplus
#endif /* TECPLOTKERNEL */
/* ENDREMOVEFROMADDON */

#ifdef MSWIN
/* BEGINREMOVEFROMADDON */
#ifdef TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#ifdef _DEBUG
#endif
#endif /* TECPLOTKERNEL */
/* ENDREMOVEFROMADDON */

#ifndef TECPLOTKERNEL
#define Widget long
#endif

#endif /* MSWIN */


#if defined UNIXX && defined ENGINE
typedef void *Widget;
#endif



#if defined SYSV || defined MSWIN
#include <string.h>
#else
#include <strings.h>
#endif

#if defined (MicrosoftC)
#include <stdlib.h>
#define EXECOS
#ifndef FAR
#define FAR
#endif
#define VOID       void
#endif

#include <sys/types.h>
#include <stdlib.h>

#if defined UNIXX
#if !defined ENGINE
#define X11
#define MOTIF
#endif
#define FAR
#define NEAR
#include <unistd.h>
#endif

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
/* OPENGL currently a must have */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#  if defined UNIXX
#  endif
#endif
/* ENDREMOVEFROMADDON */
/*
 * If not building the tecplot kernel then at least
 * include the X Instrinsics.  This will make most
 * development for addons etc work.
 */

/* NOTE: MOTIF not defined if ENGINE is defined */
#if defined MOTIF
#  if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#    if XmVERSION == 1 && XmREVISION == 0
#    endif
#  else
#    include <X11/Intrinsic.h>
#  endif
#endif

#if defined MOTIF
#define CREATE_DIALOG_PARAMS Widget W
typedef Widget ComboBoxWidget_t;
typedef Widget DropDownListWidget_t;
typedef Widget FileDialogWidget_t;
typedef Widget LabelWidget_t;
typedef Widget ListWidget_t;
typedef Widget OptionMenuWidget_t;
typedef Widget PullDownMenuWidget_t;
typedef Widget ScaleWidget_t;
typedef Widget TextFieldWidget_t;
typedef Widget ToggleWidget_t;
typedef Widget ButtonWidget_t;
typedef Widget GridWidget_t;
#endif
#if defined MSWIN
#include <windows.h>
#define CREATE_DIALOG_PARAMS     CWnd *, LaunchDialogMode_e
typedef Widget ComboBoxWidget_t;
typedef Widget DropDownListWidget_t;
typedef Widget FileDialogWidget_t;
typedef Widget LabelWidget_t;
typedef Widget ListWidget_t;
typedef Widget OptionMenuWidget_t;
typedef Widget PullDownMenuWidget_t;
typedef Widget ScaleWidget_t;
typedef Widget TextFieldWidget_t;
typedef Widget ToggleWidget_t;
typedef Widget ButtonWidget_t;
typedef Widget GridWidget_t;
#endif

/* BEGINREMOVEFROMADDON */
#if defined MSWIN && defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined TRACE
#endif
#if defined NDEBUG && !defined (CHECKED_BUILD)
#else
#endif
#endif /* MSWIN */
/* ENDREMOVEFROMADDON */


#if !defined (TRACE) /* Assume that if TRACE is not defined, then none of the TRACE macros are */
#  if !defined (NDEBUG) /* DEBUG */
#    if defined (MSWIN)
        /* If the add-on is running in debug mode but does not
         * use MFC, then no TRACE macro is available. Thus, to make tracing available,
         * map TRACE to the win32 OutpuDebugString() function.
         */
#       define TRACE(str)           do { OutputDebugString(str); } while (0)
#       define TRACE1(str,a1)       do { char s[5000]; sprintf(s,str,a1);       OutputDebugString(s); } while (0)
#       define TRACE2(str,a1,a2)    do { char s[5000]; sprintf(s,str,a1,a2);    OutputDebugString(s); } while (0)
#       define TRACE3(str,a1,a2,a3) do { char s[5000]; sprintf(s,str,a1,a2,a3); OutputDebugString(s); } while (0)
#       define TRACE0(str) TRACE(str)
#    else /* !MSWIN */
#       define TRACE(str)           do {printf(str);          fflush(stdout);} while (0)
#       define TRACE1(str,a1)       do {printf(str,a1);       fflush(stdout);} while (0)
#       define TRACE2(str,a1,a2)    do {printf(str,a1,a2);    fflush(stdout);} while (0)
#       define TRACE3(str,a1,a2,a3) do {printf(str,a1,a2,a3); fflush(stdout);} while (0)
#       define TRACE0(str) TRACE(str)
#    endif
#  else /* !DEBUG */
#    define TRACE(str)           ((void)0)
#    define TRACE0(str)          ((void)0)
#    define TRACE1(str,a1)       ((void)0)
#    define TRACE2(str,a1,a2)    ((void)0)
#    define TRACE3(str,a1,a2,a3) ((void)0)
#  endif
#endif /* !defined (TRACE) */


/*
  Platform independent way for add-ons to know how much space
  to allocate for a filename.
*/
#if !defined MAX_SIZEOFUTF8CHAR

#ifdef MS_UNICODE
  #define MAX_SIZEOFUTF8CHAR 6
#else
  #define MAX_SIZEOFUTF8CHAR 1
#endif

#endif

#if !defined (MaxCharsFilePath)
# if defined (MSWIN)
#   define MaxCharsFilePath (_MAX_PATH*MAX_SIZEOFUTF8CHAR+1) /* Includes traling '\0' */
# else
#   define MaxCharsFilePath 2047 /* ...not really a hard limit for Linux/Unix */
# endif /* MSWIN */
#endif /* !MaxCharsFilePath */

/* BEGINREMOVEFROMADDON */

/* 
 * Under Windows, if we are doing a release build (NDEBUG) that is not a CHECKED_BUILD
 * then NO_ASSERTS should be defined
 */
#if defined MSWIN && defined NDEBUG && !defined NO_ASSERTS && !defined CHECKED_BUILD
   /* intentionally break the compile */
#  error "define NO_ASSERTS for release builds"
#endif

/* 
 * Under Windows, if we are doing a CHECKED_BUILD then it should
 * also be a release build (NDEBUG)
 */
#if defined MSWIN && defined CHECKED_BUILD && !defined NDEBUG
#  error "CHECKED_BUILDS must also be release builds"
#endif


#if defined NO_ASSERTS
#  define USE_MACROS_FOR_FUNCTIONS
#endif
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
/*
 * Under Linux the definition of NULL has a cast that conflicts with our own
 * casting causing warnings that make it tough to find real problems.
 */
#if defined LINUX && defined NULL
# undef NULL
# define NULL 0
#endif
/* ENDREMOVEFROMADDON */
#endif


/* BEGINREMOVEFROMADDON */
 /* 
  * The following define allows symbols in the .h files to be
  * used by the compiler but hidden from doxygen.
  */
/* ENDREMOVEFROMADDON */
#define COMPILE_BUT_HIDE_FROM_DOC
