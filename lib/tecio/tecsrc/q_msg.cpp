#include "stdafx.h"
#include "MASTER.h"

#define TECPLOTENGINEMODULE

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2005 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#define Q_MSGMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"
#include "ARRLIST.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#if defined MOTIF
#endif
#if defined MSWIN
#endif
#endif
#endif

#include "STRUTIL.h"


#define MAXCHARACTERSPERLINE 60
/*
 * wrap a string so it contains at most CharactersPerLine
 * characters.  embedded newlines are left alone.  Spaces
 * following newlines are also left alone.
 */
Boolean_t WrapString(const char  *OldString,
                     char       **NewString)
{
  size_t L;
  if (OldString == NULL)
    return (FALSE);

  /*
   *  Assume Old string has ample spaces.  Will only be
   *  replacing some spaces with newlines and removing
   *  other spaces.  New string can be allocated to be
   *  same length as old string.
   */

  L = strlen (OldString);
  *NewString = ALLOC_ARRAY(L+1,char,"new error message string");
  if (*NewString == NULL)
    return (FALSE);

  strcpy(*NewString,OldString);

  if (L > MAXCHARACTERSPERLINE)
    {
      char *LineStart = *NewString;
      char *LastWord  = *NewString;
      char *WPtr      = *NewString;
      while (WPtr && (*WPtr != '\0'))
        {
          size_t CurLen;
          /*
           * Find next hard newline.  If there is one befor the
           * line should be chopped then reset the Last Word to
           * be at the first word after the newline.
           */
          WPtr = strchr(LineStart,'\n');
          if (WPtr && ((WPtr-LineStart) < MAXCHARACTERSPERLINE))
            {
              WPtr++;
              while (*WPtr == '\n')
                WPtr++;
              LineStart = WPtr;
              /*
               * Skip over trailing spaces.  Special handling to
               * allow indent after hard newline.
               */
              while (*WPtr == ' ')
                WPtr++;
              LastWord  = WPtr;
              continue;
            }
          /*
           *  Find next "word"
           */
          WPtr = strchr(LastWord,' ');
          if (WPtr != NULL)
            {
              while (*WPtr == ' ')
                WPtr++;
            }

          if (WPtr == NULL)
            {
              CurLen = strlen(LineStart);
            }
          else
            {
              CurLen = WPtr-LineStart;
            }

          if (CurLen > MAXCHARACTERSPERLINE)
            {
              /*
               * Line is too long.  Back up to previous
               * word and replace preceeding space with
               * a newline.
               */
              if (LastWord == LineStart)
                {
                  /*
                   * Bad news, line has very long word.
                   */
                  if (WPtr && (*WPtr != '\0'))
                    {
                      *(WPtr-1) = '\n';
                      LastWord = WPtr;
                    }
                }
              else
                {
                  *(LastWord-1) = '\n';
                }
              LineStart = LastWord;
            }
          else
            LastWord = WPtr;
        }
    }
  return (TRUE);
}


static void SendWarningToFile(FILE *F,
                              const char *S)
{
  char *S2;
  REQUIRE(VALID_REF(F));
  REQUIRE(VALID_REF(S));
  if (WrapString(S,&S2))
    {
      fprintf(F,"Warning: %s\n",S2);
      FREE_ARRAY(S2,"temp warning string");
    }
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 * Show the warning message. Note that the Format string can be the only
 * argument, in which case it is essentially the warning message itself.
 *
 * param Format
 *     C Format string or a simple message.
 * param ...
 *     Zero or more variable arguments. The number of arguments must correspond
 *     to the placeholders in the format string.
 */
void Warning(TranslatedString_pa Format,
             ...) /* zero or more arguments */
{
  REQUIRE(VALID_REF(Format));

  static Boolean_t InWarning = FALSE; /* ...used to prevent recursive deadlock */
  if (!InWarning)
    {
      #if defined TECPLOTKERNEL
        if ( Thread_ThreadingIsInitialized() )
          Thread_LockMutex(ErrMsgWarningMutex);
      #endif

      InWarning = TRUE;
        {
          char    *MessageString;
          va_list  Arguments;
        
          va_start(Arguments, Format);
        
          /*
           * Attempt to format the string. Failing that simply use the original format
           * string argument which, if we ran out of memory while formatting, is
           * probably just an warning message saying that we ran out of memory in some
           * previous operation anyway.
           */
          MessageString = vFormatString(_ts2char(Format), Arguments);
          if (MessageString == NULL)
            MessageString = (char *)_ts2char(Format);
        
          va_end(Arguments);
        
          #if defined TECPLOTKERNEL
            {
              #ifdef MSWIN
                {
                  /* Always use const char * for TRACE("%s") */
                  
                  if ( strlen(MessageString) < 512 )
                    TRACE("Warning: %s\n",_CSTR(MessageString));
                }
              #endif
        
              if (!IsProcessingMacroCommand && GeneralBase.Interface.EnableWarnings)
                {
                  if (InBatchMode)
                    { 
                      if (WriteBatchLog)
                        SendWarningToFile(BatchFile,MessageString);
                    }
                  else if (!InMiddleOfInteractiveViewChange())
                    {
                      #if defined UNIXX
                        {
                          if (GraphicsAreInitialized)
                            {
                              MainlineInvoke(ShowWarningMsgDialogJob, MessageString);
                              if (GeneralBase.Interface.PrintDebug)
                                SendWarningToFile(StdError,MessageString);
                            }
                          else
                            {
                              SendWarningToFile(stderr,MessageString);
                            }
                        }
                      #endif
                      #if defined MSWIN
                        {
                          MainlineInvoke(ShowWarningMsgDialogJob, (void*)_nts(MessageString));
                          if (GeneralBase.Interface.PrintDebug)
                            SendWarningToFile(StdError,MessageString);
                        }
                      #endif
                    }
                  else
                    SendWarningToFile(StdError,MessageString);
                }
            }
          #else /* !TECPLOTKERNEL */
            {
              SendWarningToFile(stderr,MessageString);
            }
          #endif
        
          /* cleanup if we allocated the string */
          if (MessageString != _ts2char(Format))
            FREE_ARRAY(MessageString, "MessageString");
        
        }
      InWarning = FALSE;

      #if defined TECPLOTKERNEL
        if ( Thread_ThreadingIsInitialized() )
          Thread_UnlockMutex(ErrMsgWarningMutex);
      #endif
    }
}


static void SendErrToFile(FILE       *File,
                          const char *Msg)
{
  char *FormattedMsg;
  REQUIRE(VALID_REF(File));
  REQUIRE(VALID_REF(Msg));
  if ( WrapString(Msg,&FormattedMsg) )
    {
      fprintf(File, "Err: %s\n", FormattedMsg);
      FREE_ARRAY(FormattedMsg, "temp error string");
    }
  else
    fprintf(File, "Err: %s\n", Msg);
}


/* Fall-back ErrMsg procedure when nothing else works */
static void DefaultErrMsg(char *Msg)
{
  REQUIRE(VALID_REF(Msg));
  #ifdef MSWIN
    #ifdef TECPLOTKERNEL
      MessageBox(NULL, _WSTR(Msg), _ts2tchar(_ts("Tecplot Error")), MB_OK|MB_ICONERROR);
    #else
      MessageBox(NULL, Msg, "Error", MB_OK|MB_ICONERROR);
    #endif
  #else
    SendErrToFile(stderr, Msg);
  #endif
}

/**
 * Show the error message. Note that the Format string can be the only
 * argument, in which case it is essentially the error message itself.
 *
 * @param Format
 *   C Format string or a simple message.
 * @param ...
 *   Zero or more variable arguments. The number of arguments must correspond
 *   to the placeholders in the format string.
 *
 * NOTES:
 *   This function is thread safe in that it may be safely called by multiple
 *   threads however when running interactively only the first error message is
 *   queued for display on idle. In batch mode all messages are sent to the
 *   batch log file.
 */
void ErrMsg(TranslatedString_pa Format,
            ...) /* zero or more arguments */
{
  REQUIRE(VALID_REF(Format));

  static Boolean_t InErrMsg = FALSE; /* ...used to prevent recursive deadlock */
  if (!InErrMsg)
    {
      #if defined TECPLOTKERNEL
        if ( Thread_ThreadingIsInitialized() )
          Thread_LockMutex(ErrMsgWarningMutex);
      #endif
    
      InErrMsg = TRUE;
        {
          char    *MessageString;
          va_list  Arguments;
    
          va_start(Arguments, Format);
    
          /*
           * Attempt to format the string. Failing that simply use the original format
           * string argument which, if we ran out of memory while formatting, is
           * probably just an error message saying that we ran out of memory in some
           * previous operation anyway.
           */
          MessageString = vFormatString(_ts2char(Format), Arguments);
          if (MessageString == NULL)
            MessageString = (char *)_ts2char(Format);
    
          va_end(Arguments);
    
          #if defined TECPLOTKERNEL
            {
              #ifdef MSWIN
                {
                  /* TRACE will assert if the message is too big (in this case 512 bytes)...*/
                  if ( strlen(MessageString) < 512 )
                    TRACE("Error: %s\n",MessageString);
                }
              #endif
    
              if (InBatchMode)
                {
                  /* There may be an error before BatchFile is opened or after */
                  /* BatchFile is closed */
                  if ( BatchFile )
                    {
                      if (WriteBatchLog)
                        SendErrToFile(BatchFile,MessageString);
                    }
                  else
                    DefaultErrMsg(MessageString);
                }
              else if (!InMiddleOfInteractiveViewChange())
                {
                  #if defined UNIXX
                    {
                      #if !defined ENGINE
                      /* TODO(BDP)-M 01/10/2006 ENGINE:  May need to make this a global
                       *                                 state that can be queried by the
                       *                                 application via TecUtil call.
                       */
                      PendingErrors = TRUE;
                      #endif
                      if (GraphicsAreInitialized)
                        {
                          MainlineInvoke(ShowErrMsgDialogJob, MessageString);
                          if (GeneralBase.Interface.PrintDebug)
                            SendErrToFile(StdError,MessageString);
                        }
                      else
                        SendErrToFile(stderr,MessageString);
                    }
                  #elif defined MSWIN
                    {
                      MainlineInvoke(ShowErrMsgDialogJob, (void*) _nts(MessageString));
                      if (GeneralBase.Interface.PrintDebug)
                        SendErrToFile(StdError,MessageString);
                    }
                  #endif
                }
              else
                DefaultErrMsg(MessageString);
            }
          #else /* !TECPLOTKERNEL */
            {
              DefaultErrMsg(MessageString);
            }
          #endif
    
          /* cleanup if we allocated the string */
          if (MessageString != _ts2char(Format))
            FREE_ARRAY(MessageString, "MessageString");
        }
      InErrMsg = FALSE;
    
      #if defined TECPLOTKERNEL
        if ( Thread_ThreadingIsInitialized() )
          Thread_UnlockMutex(ErrMsgWarningMutex);
      #endif
    }
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#endif
#if !defined ENGINE
#endif
#if !defined ENGINE
#endif
#if !defined ENGINE
#if defined MOTIF
#endif /* MOTIF */
#endif
#if !defined ENGINE
#endif
#endif /* TECPLOTKERNEL */

