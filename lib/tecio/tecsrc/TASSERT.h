/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2004 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
/*
 * Provide four levels of assertion control. Assertions provide a mechanism
 * to enforce a contract between a client and service provider. The assertions
 * are listed in order of highest to lowest priority. Assertions can be turned
 * off individually by defining the appropriate name (see preprossessor
 * definitions below), however, lower priority assertions should be turned
 * off prior to higher ones. As confidence in the code increases all assertions
 * can be turned off by defining NO_ASSERTS.
 *
 * The assertions defined below have the following meanings:
 *
 *     INVARIANT - Asserts that a property's state is invariant throughout the
 *                 life of the property's scope. For instance in Tecplot
 *                 the 'CurFrame' variable has global scope and must always
 *                 pass the VALID_REF test both before and after any public
 *                 method call. Stating invariant properties of an application
 *                 provides a deeper understanding of the application's state.
 *                 These statements are usually positioned just ahead of the
 *                 preconditions and just after the postconditions.
 *
 *     REQUIRE   - Asserts that a method's preconditions are within their
 *                 valid domains. Preconditions are conditions placed upon
 *                 any state information relied upon for the call. These
 *                 statements should be as close to the top of the method
 *                 as possible (except for assertions on invariant properties).
 *
 *     ENSURE    - Asserts that a method's postconditions are within their
 *                 valid ranges. Postconditions are conditions placed upon
 *                 any state information modified by the call. These
 *                 statements should be as close to the bottom of the method
 *                 (presumably there is only one exit point) as possible
 *                 (except for assertions on invariant properties).
 *
 *     CHECK     - Any other assertion not covered by the above assertions.
 *                 These are often added within a method body to specify
 *                 something that may not be immediately obvious to the reader
 *                 or to validate your assumptions about a call to a 3rd party
 *                 method that does not use runtime assertions for its
 *                 preconditions or postconditions. Obviously if the 3rd party
 *                 method uses assertions then there is no need for the CHECK.
 *
 * Additionally a convenience macro is available to place in code that is
 * pending implementation.
 *
 *     NOT_IMPLEMENTED - Assertion that always fails during runtime for debug
 *                       builds and always fails at compile time for release
 *                       builds.
 */
#if !defined TASSERT_H
#define TASSERT_H

#if defined (MSWIN)
# include <assert.h>
#endif /* MSWIN */

#if !defined TECPLOTKERNEL && !defined STD_ASSERTS
#define STD_ASSERTS
#endif

#if !defined (MSWIN)
#  include <assert.h>
#  if !defined ASSERT
#    define ASSERT assert
#  endif
#endif

/* BEGINREMOVEFROMADDON */
#define INVALID_REF       ((void *)0x0000FFFF)
/*
 * Chances are low the address 0x11111111 will be used, so we'll risk asserting
 * against it (see unitialized assignment in newmalloc).
 */
#define UNINITIALIZED_REF ((void *)0x11111111)
#define INVALID_FN_REF    ((void *)NULL)
/* ENDREMOVEFROMADDON */

#ifdef UNIXX
/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#  if defined NO_ASSERTS
#  else
#  endif
#endif /* TECPLOTKERNAL */
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
#if !defined TECPLOTKERNEL
/* For add-ons, there is a problem with VALID_REF, so just test for non-NULL */
/* ENDREMOVEFROMADDON */
#  define VALID_REF(p)      ( (p)  != NULL )
#  define VALID_FN_REF(fp)  ( (fp) != NULL )
/* BEGINREMOVEFROMADDON */
#endif /* !defined TECPLOTKERNAL */
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
   /* Widgets are pointers under Motif */
# define VALID_WIDGET(widget)       VALID_REF((widget))
  /* Menu widgets are pointers too */
# define VALID_MENU_WIDGET(widget)   VALID_REF((widget))
/* ENDREMOVEFROMADDON */
#endif /* UNIXX */

#ifdef MSWIN
   /* Under Windows, we can use AfxIsValidAddress to check for */
   /* correct addresses with in the program.  In the future,   */
   /* could use AfxIsValidMemBlock to check for valid heap     */
   /* addresses when we know we are to be looking at the heap. */
/* BEGINREMOVEFROMADDON */
   /* For Windows, we want to be able to make a dll (i.e., tecio.dll) */
   /* from the core source code w/o having to use MFC, so we need     */
   /* a generic WIN32 equivalent if MFC is not defined                */
/* ENDREMOVEFROMADDON */
# if defined(_AFX)
/* VC++ 6.0 puts literal strings (stuff in "") in read-only memory */
/* #define VALID_REF(pointer)         AfxIsValidAddress((void*)(pointer),1,TRUE) */
#define VALID_REF(pointer)         AfxIsValidAddress((void*)(pointer),1,FALSE)
#   define VALID_FN_REF(fn_pointer)   AfxIsValidAddress((fn_pointer),1,FALSE)
# else /* !_AFX */
/* VC++ 6.0 puts literal strings (stuff in "") in read-only memory */
/* #   define VALID_REF(p)     ( p != NULL && !IsBadReadPtr((CONST VOID *)p,1) && !IsBadWritePtr((LPVOID)p,1) ) */
#   define VALID_REF(p)     ( p != NULL && !IsBadReadPtr((CONST VOID *)p,1) )
#   define VALID_FN_REF(pf) ( pf != NULL && !IsBadReadPtr((CONST VOID *)pf,1) )
# endif /* _AFX */
/* BEGINREMOVEFROMADDON */
   /* Widgets are numbers under Windows, so we decode it with GetWindowFromWidget */
# if defined ENGINE
#   define VALID_WIDGET(widget)       ((widget) != NULL)
# else 
#   define VALID_WIDGET(widget)       ((widget) != NULL && GetWindowFromWidget((widget))!=NULL)
# endif // ENGINE


  /* Menu widgets are numbers too, so we just check against zero */
# define VALID_MENU_WIDGET(widget)  ((widget)!=NULL)
/* ENDREMOVEFROMADDON */
#endif /* MSWIN */
/* BEGINREMOVEFROMADDON */
/* handles are not pointers to memory, so the only test we can */
/* perform is to check for 0 */
#define VALID_HANDLE(handle)       ((handle)!=0)
/* ENDREMOVEFROMADDON */
/* other useful validity checks */
#define VALID_BOOLEAN(b)           ((b) == TRUE || (b) == FALSE)
#define VALID_ENUM(value, type)    (0 <= (value) && \
                                         (value) < END_##type)

/* Test a parameter than can be NULL or a valid pointer */
#define VALID_REF_OR_NULL(ptr) IMPLICATION((ptr) != NULL, VALID_REF(ptr))
#define VALID_FN_REF_OR_NULL(ptr) IMPLICATION((ptr) != NULL, VALID_FN_REF(ptr))

/* BEGINREMOVEFROMADDON */
/**
 * These macros are a little complicated but it allows one to
 * write a simple assertion regardless of the zone type or
 * selected plane:
 *
 *   REQUIRE(VALID_CELL_INDEX(CZData, CellIndex, Plane)));
 *
 * Prior to using the macros a call to SetupXxx,
 * or at a minimum SetupCZData, must be called to setup
 * the globals defining the dataset structure.
 */
#define VALID_FE_CELL_INDEX(CZData, CellIndex) \
            (/* CellIndex range test */ \
             0 <= (CellIndex) && \
                  (CellIndex) < (CZData)->NumElements)

#define VALID_IPLANE_CELL_INDEX(CZData,CellIndex) \
            (/* CellIndex range test */ \
             (CellIndex) >= 0 && \
             IINDEX((CZData),CellIndex) <= MAX((CZData)->NumIPtsM1,1) && \
             JINDEX((CZData),CellIndex) <  MAX((CZData)->NumJPtsM1,1) && \
             KINDEX((CZData),CellIndex) <  MAX((CZData)->NumKPtsM1,1))

#define VALID_JPLANE_CELL_INDEX(CZData,CellIndex) \
            (/* CellIndex range test */ \
             (CellIndex) >= 0 && \
             IINDEX((CZData),CellIndex) <  MAX((CZData)->NumIPtsM1,1) && \
             JINDEX((CZData),CellIndex) <= MAX((CZData)->NumJPtsM1,1) && \
             KINDEX((CZData),CellIndex) <  MAX((CZData)->NumKPtsM1,1))

#define VALID_KPLANE_CELL_INDEX(CZData,CellIndex) \
            (/* CellIndex range test */ \
             (CellIndex) >= 0 && \
             IINDEX((CZData),CellIndex) <  MAX((CZData)->NumIPtsM1,1) && \
             JINDEX((CZData),CellIndex) <  MAX((CZData)->NumJPtsM1,1) && \
             KINDEX((CZData),CellIndex) <= MAX((CZData)->NumKPtsM1,1))

#define VALID_ORDERED_CELL_INDEX(CZData, CellIndex, Plane) \
            (/* macro preconditions */ \
             ((IJKPlanes_e)(Plane) == Planes_I || \
              (IJKPlanes_e)(Plane) == Planes_J || \
              (IJKPlanes_e)(Plane) == Planes_K || \
              (IJKPlanes_e)(Plane) == Planes_Volume) && \
\
             /* CellIndex range test */ \
             (IMPLICATION(((IJKPlanes_e)(Plane) == Planes_I || \
                           (IJKPlanes_e)(Plane) == Planes_Volume), \
                          VALID_IPLANE_CELL_INDEX((CZData),CellIndex)) && \
              IMPLICATION(((IJKPlanes_e)(Plane) == Planes_J || \
                           (IJKPlanes_e)(Plane) == Planes_Volume), \
                          VALID_JPLANE_CELL_INDEX((CZData),CellIndex)) && \
              IMPLICATION(((IJKPlanes_e)(Plane) == Planes_K || \
                           (IJKPlanes_e)(Plane) == Planes_Volume), \
                          VALID_KPLANE_CELL_INDEX((CZData),CellIndex))))

#define VALID_CELL_INDEX(CZData, CellIndex, Plane) \
           ((CZData)->NM != NULL ? \
              VALID_FE_CELL_INDEX((CZData), (CellIndex)) : \
              VALID_ORDERED_CELL_INDEX((CZData), (CellIndex), (Plane)))

#define VALID_DATASET(D) (((D) != NULL) && ((D)->NumZones > 0))



#if defined MSWIN
  #if defined CHECKED_BUILD

    extern void TWinCheckedFailedLine(const char *expr,
                                           const char *filename,
                                           int LineNum);
  #endif /* CHECKED_BUILD */

   /* Here is a more specific check in Windows for a valid
      pointer to an MFC Window object. 
      Note that GetSafeHwnd() works even if pWnd is NULL, because
      it checks the 'this' pointer first */
  #define VALID_WND(pWnd) (::IsWindow((pWnd)->GetSafeHwnd())) 

#else /* !MSWIN */
  #define VALID_WND(pWnd) /* Should not be used in Motif */ 
#endif /* MSWIN */
/* ENDREMOVEFROMADDON */

/* Check for a non-zero length string */
#define VALID_NON_ZERO_LEN_STR(str) (VALID_REF(str) && !ISEMPTYSTRING(str) )

/* Check for valid stdio file handle */
#define VALID_FILE_HANDLE(stream) ((stream) != NULL)

/* To check colors and pen numbers */
/* BEGINREMOVEFROMADDON */
#define VALID_BASIC_COLOR(BColor) \
          (FirstBasicColor<=(BColor) && (BColor)<=LastBasicColor)
#define VALID_CONTOUR_COLOR(Color) \
          (ContourColorOffset<=(Color) && \
           (Color)<ContourColorOffset+GeneralBase.Limits.MaxNumContourLevels+1)
#define VALID_PLOTTING_COLOR(Color) \
          (VALID_BASIC_COLOR(Color) || VALID_CONTOUR_COLOR(Color))
#define VALID_INTERFACE_SPECIFIC_COLOR(BColor) \
          (FirstInterfaceColor<=(BColor) && (BColor)<=LastInterfaceColor)
#define VALID_INTERFACE_COLOR(Color) \
          (VALID_PLOTTING_COLOR(Color) || VALID_INTERFACE_SPECIFIC_COLOR(Color))
#define VALID_MULTICOLOR_COLOR(Color) \
          (((Color) == MultiColor_C) || ((Color) == MultiColor2_C) || \
           ((Color) == MultiColor3_C) || ((Color) == MultiColor4_C))
#define VALID_RGB_COLOR(Color) \
          ((Color) == RGBColor_C)
#define VALID_ASSIGNABLE_COLOR(C) \
        (VALID_BASIC_COLOR(C)      || \
         VALID_MULTICOLOR_COLOR(C) || \
         VALID_RGB_COLOR(C))
#define VALID_PEN_OFFSET(PenOffset) \
          (Black_C<=(PenOffset) && (PenOffset)<=NumPlotterPens)
#define VALID_PEN_OFFSET_FOR_OBJECT(PenOffset) \
          (FirstObjectPen<=(PenOffset) && (PenOffset)<=LastObjectPen)


/* to check FE cells */
#define VALID_ELEMENT_TYPE(element_type) \
          ((element_type) == ZoneType_FETriangle || \
           (element_type) == ZoneType_FEQuad     || \
           (element_type) == ZoneType_FETetra    || \
           (element_type) == ZoneType_FEBrick    || \
           (element_type) == ZoneType_FELineSeg)



/*
 * Test validity of zone and variable names. A valid name is one that has a
 * valid reference, is not padded with spaces and is within the maximum
 * specified length.
 */
#define VALID_NAME(Name, MaxLength) \
          (VALID_REF(Name) && \
           (ISEMPTYSTRING(Name) || \
            (!isspace((Name)[0]) && !isspace((Name)[strlen(Name)-1]))) && \
           strlen(Name) <= (MaxLength))
#define VALID_ZONE_NAME(Name) VALID_NAME((Name), MaxChrsZnTitle)
#define VALID_VAR_NAME(Name)  VALID_NAME((Name), MaxChrsVarName)


/* Special test for lighting effect (don't allow "none" in some cases) */
#define VALID_LIGHTINGEFFECT(L) \
          (((L) == LightingEffect_Paneled) || ((L) == LightingEffect_Gouraud))


/* type definition for assert failure notification function */
typedef void (*TAssertFailureNotifyFunc)(
    const char *expression, /* text representation of the assertion */
    const char *file_name,  /* name of the file containing the assertion */
    int        line);       /* line number in the file of the assertion */

#if !defined STD_ASSERTS
/* external function prototypes */
extern void TAssert(
    const char *expression, /* text representation of the assertion */
    const char *file_name,  /* name of the file containing the assertion */
    int        line);       /* line number in the file of the assertion */

extern TAssertFailureNotifyFunc InstallTAssertFailureNotify(
    TAssertFailureNotifyFunc new_function); /* new notification function */
#endif /* !STD_ASSERTS */
/* ENDREMOVEFROMADDON */

#if defined NO_ASSERTS
/* BEGINREMOVEFROMADDON */
#   define TASSERT(EXPR)
/* ENDREMOVEFROMADDON */
#   define INVARIANT(EXPR)
#   define REQUIRE(EXPR)
#   define ENSURE(EXPR)
#   define CHECK(EXPR)
#   ifdef VERIFY
#     undef VERIFY
#   endif
#   define VERIFY(EXPR)    ((void)(EXPR))
    /* 
     * Only define IGNORENOTIMPLEMENTED if building a "test" release build
     * that you are fully aware may contain unimplemented features.
     */
#   if defined IGNORENOTIMPLEMENTED
#     define NOT_IMPLEMENTED() CHECK(FALSE)
#   else
#     if defined MSWIN
      /* 
       * NOT_IMPLEMENTED is defined using a parameter, but should be called with none,
       * this will then throw a warning and not break the compile. Unix doesn't pick
       * up this warning, so break the compile under Unix
       */
#       define NOT_IMPLEMENTED(x)  ASSERT("Not Implemented")
#     endif
#     if defined UNIXX
#       define NOT_IMPLEMENTED()  not implemented /* intentionally break the compile */
#     endif
#   endif
#elif defined STD_ASSERTS
/* BEGINREMOVEFROMADDON */
#   define TASSERT(EXPR)         assert(EXPR)
/* ENDREMOVEFROMADDON */
#   define INVARIANT(EXPR)       assert(EXPR)
#   define REQUIRE(EXPR)         assert(EXPR)
#   define ENSURE(EXPR)          assert(EXPR)
#   define CHECK(EXPR)           assert(EXPR) 
#   ifndef VERIFY
#     define VERIFY(EXPR)        assert(EXPR)
#   endif /* VERIFY */
#   define NOT_IMPLEMENTED()     assert(!("Not Implemented"))
#else
/* BEGINREMOVEFROMADDON */
#   if defined MSWIN
#     if defined(NDEBUG)
#       if defined (CHECKED_BUILD)
#         define TASSERT(EXPR) ((void)((EXPR) || (TWinCheckedFailedLine(#EXPR, THIS_FILE, __LINE__),1)))
#       else
#         define TASSERT(EXPR) (void(0))
#       endif /* !CHECKED_BUILD */
#     else
#       define TASSERT(EXPR) ASSERT(EXPR)
#      endif

#   else /* !MSWIN */
#     define TASSERT(EXPR) (void)((EXPR) ||\
                                  (TAssert(#EXPR, __FILE__, __LINE__), 0))
#   endif
#   if defined NO_INVARIANTS
#     define INVARIANT(EXPR)
#   else
#     define INVARIANT(EXPR) TASSERT(EXPR)
#   endif

#   if defined NO_PRECONDITIONS
#     define REQUIRE(EXPR)
#   else
#     define REQUIRE(EXPR) TASSERT(EXPR)
#   endif

#   if defined NO_POSTCONDITIONS
#     define ENSURE(EXPR)
#   else
#     define ENSURE(EXPR) TASSERT(EXPR)
#   endif

#   if defined VERIFY
#     undef VERIFY
#   endif

#   if defined NO_CHECKS
#     define CHECK(EXPR)
#     define VERIFY(EXPR)  ((void)(EXPR))
#   else
#     define CHECK(EXPR)  TASSERT(EXPR)
#     define VERIFY(EXPR) TASSERT(EXPR)
#   endif

#   if defined NICE_NOT_IMPLEMENTED
#     define NOT_IMPLEMENTED() NiceNotImplemented()
#   else
#     define NOT_IMPLEMENTED() TASSERT(!("Not Implemented"))
#   endif
/* ENDREMOVEFROMADDON */
#endif
/* BEGINREMOVEFROMADDON */
#if !defined STD_ASSERTS
extern void TecplotMopupOnAssert(void);
#endif /* !STD_ASSERTS */

#if defined NICE_NOT_IMPLEMENTED
extern void NiceNotImplemented(void);
#endif
/* ENDREMOVEFROMADDON */

/* convenience macros for implication, P -> Q, and equivalence, P <-> Q. */
#define IMPLICATION(P,Q) (!(P) || (Q))
#define EQUIVALENCE(P,Q) ((P) == (Q))

#endif /* TASSERT_H */
