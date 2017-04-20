#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 1988-2002  by AMTEC ENGINEERING INC. *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#define TASSERTMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined (MSWIN)
#endif
#endif

#define MAX_ERRMSG_LENGTH 2096

/* the mopup from assert and the writing out of crash.lay are */
/* used by TUASSERT and thus are needed even if NO_ASSERTS */
/* is set */
#if !defined NO_TU_ASSERTS || !defined NO_ASSERTS

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#ifdef ELM
#endif
#endif /* TECPLOTKERNEL */

#endif /* Mopup function needed ... */



#if !defined STD_ASSERTS
/*
 * Define the final assertion notification function.
 */
#if defined UNIXX && !defined NO_ASSERTS
/*******************************************************************
 *                                                                 *
 *                          UNIX                                   *
 *                                                                 *
 *******************************************************************/


#  if defined NDEBUG
  /*
   * if NDEBUG is defined __assert is NOT defined so we must supply 
   * our own assertion notification function.....
   */
# define ASSERT assert
static void UnixAssert(const char *expression,
                       const char *file_name,
                       int        line)
{
  char buffer[MAX_ERRMSG_LENGTH + 1];
  fprintf(stderr,"Assertion: %s\n"
                 "Tecplot version: %s\n"
                 "File Name: %s\n"
                 "Line Number: %d\n",
          expression,TecVersionId, file_name,line);
  exit(ExitCode_AssertionFailure);
}
static TAssertFailureNotifyFunc assert_failure_notify = UnixAssert;
#  else
  /*
   * NDEBUG is not defined so __assert is available....
   */
#    if defined LINUX
#      define LOWLEVELASSERTFUNCTION __linuxassertproxy
       /*
        * In linux, __assert does not exist but rather
        * __assert_fail which has a differnt API.  Thus
        * a proxy is provided
        */
      static void __linuxassertproxy(const char *__assertion,
                                     const char *__file,
                                     int         __line)
      {
        __assert_fail(__assertion,__file,__line,__ASSERT_FUNCTION);
      }
#    elif defined DARWIN
#      define LOWLEVELASSERTFUNCTION __darwinassertproxy
       /*
        * In Darwin (Mac OS X), __assert is #defined to a call to __eprintf,
        * which also has a different API. Another proxy...
        */
      static void __darwinassertproxy(const char *__assertion,
                                      const char *__file,
                                      int         __line)
      {
        __eprintf ("Assertion: %s\n"
                   "Tecplot version: %s\n"
                   "File Name: %s\n"
                   "Line Number: %d\n",
                   __assertion, TecVersionId, __file, (unsigned)__line);
      }
#   else
#     define LOWLEVELASSERTFUNCTION __assert
#   endif

static TAssertFailureNotifyFunc assert_failure_notify = (TAssertFailureNotifyFunc) LOWLEVELASSERTFUNCTION;

#  endif
#endif /* UNIXX */

#if defined UNIXX && !defined NO_ASSERTS
/*
 * Replace the current assert failure notification function and
 * return the current one.
 *
 * Assumptions:
 *     new function points to valid function (not null) that
 *     conforms to the specified prototype
 *
 * Guarantees:
 *     result is a pointer to the previously installed
 *     function (not null)
 */
TAssertFailureNotifyFunc InstallTAssertFailureNotify(
    TAssertFailureNotifyFunc new_function) /* new notification function */
{
  TAssertFailureNotifyFunc result = 0; /* old function address */

  ASSERT(new_function != 0);

  result = assert_failure_notify;
  assert_failure_notify = new_function;

  ASSERT(result != 0);
  return result;
}



static void AddToAssertTable(const char *Message)
{
  /* No asserts, we don't want recursive calls */
  if ( AssertTable )
    {
      Boolean_t AssertionFound = FALSE;
      int       ii;
      for ( ii = 0; ii < NumAssertTableEntries; ii++ )
        {
          if ( strcmp(Message, AssertTable[ii].assertion) == 0 )
            {
              AssertionFound = TRUE;
              AssertTable[ii].count++;
              break;
            }
        }
      if ( !AssertionFound )
        {
          if ( NumAssertTableEntries < MaxAssertTableEntries )
            {
              AssertTable[NumAssertTableEntries].assertion = NULL;
              AssertTable[NumAssertTableEntries].assertion = (char *)ALLOC_ARRAY(strlen(Message)+1, 
                                                                                 char, 
                                                                                 "AssertTable.assertion");
              if ( AssertTable[NumAssertTableEntries].assertion != NULL )
                {              
                  strcpy(AssertTable[NumAssertTableEntries].assertion, Message);
                  AssertTable[NumAssertTableEntries].count = 1;
                  NumAssertTableEntries++;
                }
            }
          else
            {
              static Boolean_t ShowOnlyOnce = TRUE;
              if ( ShowOnlyOnce )
                {
                  const char *MsgString = "The allowed amount of Assertion "
                                          "Entries has been exceeded. No more "
                                          "assertions will be recorded";
#                 if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#                 else
                    fprintf(stderr, MsgString);
#                 endif
                }
              ShowOnlyOnce = FALSE;
            }
        }
    }
}



/*
 * Perform the installed assert failure notification action.
 */
void TAssert(const char *expression, /* text representation of the assertion */
             const char *file_name,  /* name of the file containing the assertion */
             int        line)        /* line number in the file of the assertion */
{
  static Boolean_t InTAssert = FALSE;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
  char Message[MAX_ERRMSG_LENGTH+1];
  ASSERT(expression != 0 && strlen(expression) != 0);
  ASSERT(file_name != 0 && strlen(file_name) != 0);
  ASSERT(line >= 1);

  /* check for recursion */
  if (InTAssert)
    {
      fprintf(stderr, "Already in assert!\n");
      fprintf(stderr,"Assertion: %s\n"
                     "Tecplot version: %s\n"
                     "File Name: %s\n"
                     "Line Number: %d\n",
                     expression, TecVersionId, file_name, line);
      PrintCurBacktrace(stderr, 100);
      ASSERT(FALSE); /*... really exit */
    }

  InTAssert = TRUE;

  sprintf(Message, "Assertion: %s\n"
                   "Tecplot version: %s\n"
                   "File Name: %s\n"
                   "Line Number: %d\n",
                    expression, TecVersionId, file_name, line);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined (MSWIN)
  TWinCopyAssertToClipboard(expression,file_name,line);
#endif

  if ( AssertTable )    
    AddToAssertTable(Message);
  else
    {
#     if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#     else
        fprintf(stderr, "%s", Message);  
#     endif
    }

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

      (*assert_failure_notify)(expression, file_name, line);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

  InTAssert = FALSE; /* just in case assert_failure_notify has an ignore */
}
#endif /* defined MSWIN || (defined UNIXX && !defined NO_ASSERTS) */
#endif /* STD_ASSERTS */


#if defined MSWIN && defined CHECKED_BUILD

void TWinCheckedFailedLine(const char *Expr,
                           const char *FileName,
                           int LineNum)
{
  static bool bInTWinCheckedFailedLine = false;

  if ( !bInTWinCheckedFailedLine ) /* Avoid recursion */
    {
      bInTWinCheckedFailedLine = true;

      CString strAssert;
      strAssert.Format(_T("Checked build assertion failure:\r\n%s\r\n")
                       _T("Tecplot version: %s\r\n")
                       _T("File Name: %s\r\n")
                       _T("Line Number: %d\r\n"),
                       Utf8Api::Utf8ToWideChar(Expr),
                       Utf8Api::Utf8ToWideChar(TecVersionId),
                       Utf8Api::Utf8ToWideChar(FileName),
                       LineNum);

  
      if ( ::OpenClipboard( ::IsWindow(TPMainWnd->GetSafeHwnd()) ? TPMainWnd->GetSafeHwnd() : NULL ) )
        {
            ::EmptyClipboard();

            CSharedFile sf(GMEM_MOVEABLE|GMEM_DDESHARE|GMEM_ZEROINIT);
            sf.Write(strAssert,strAssert.GetLength()*sizeof(TCHAR));
            HGLOBAL hMem = sf.Detach();

            if ( hMem )
              {
                SetClipboardData(CF_UNICODETEXT,
                                 hMem);

                /* Don't release hMem since the clipbord now owns it */

              }

            ::CloseClipboard();
        }


      ::MessageBeep(0);
      if ( ::MessageBox(NULL,
                   strAssert,
                   _T("Checked build assertion FAILURE!"),
                   MB_OKCANCEL|MB_ICONSTOP|MB_TASKMODAL) != IDOK)
        {
          exit(-1); /* We're outta here! */
        }
    } /* !bInTWinCheckedFailedLine */
}
#endif /* MSWIN && CHECKED_BUILD */


#if defined NICE_NOT_IMPLEMENTED
static Boolean_t NotImplementedCalled = FALSE;
void NiceNotImplemented(void)
{
  if (!NotImplementedCalled)
    {
      Warning("Not Implemented!");
      NotImplementedCalled = TRUE;
    }
}
#endif
