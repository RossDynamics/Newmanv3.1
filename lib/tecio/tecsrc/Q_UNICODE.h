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
  How to use the translate string functions:

  +---------------+---------------------+---------------------------+---------------------------------------------+
  | Function Name | Converts from:      |To:                        | Notes                                       |
  |===============+=====================+===========================+=============================================+
  | _ts(s)        | const char*         | TranslatedString_pa       | Use for any string needing to be translated |
  +-------------- +---------------------+---------------------------+---------------------------------------------+
  | _nts(s)       | const char*         | TranslatedString_pa       | Use for any string NOT to be translated     |
  +---------------+---------------------+---------------------------+---------------------------------------------+
  | _ts_static(s) | static const char * | static TranslatedString_pa| Same as _ts, but use for static strings only|
  +---------------+---------------------+---------------------------+---------------------------------------------+
  | _ts2char(s)   | TranslatedString_pa | const char *              |                                             |
  +---------------+---------------------+---------------------------+---------------------------------------------+
  | _WCHAR(s)     | TranslatedString_pa | const wchar_t *           | Used only for Windows API                   |
  +---------------+---------------------+---------------------------+---------------------------------------------+

Notes:

    1. Functions may be nested. For example:
      _ts2char(_ts("Message to be tranlsated"));

    char Buffer;
    sprintf(Buffer, "The '%s' is translated as '%s'",
            "Hello"
            _ts2char(_ts("Hello"));

    Buffer will be filled with something like "The word 'Hello' is translated as 'Hola'."


    2. The _ts and _nts functions have an optional second parameter which is any notes you wish
       to give the human translator. For example:

      _ts("Error in Frame 1","'Frame' is a Tecplot frame");


    3. _ts_static is used in rare cases where you need to decare a string message which must be
       valid for the duration of the Tecplot session. Typically these are declared at the top of the
       file where the message is displayed multiple times:

       static TranslatedString_pa ErrMsg = _ts_static("Out of memory");

      
*/

#if !defined Q_UNICODE_H_
# define Q_UNICODE_H_

#if defined EXTERN
  #undef EXTERN
#endif
#if defined Q_UNICODEMODULE
  #define EXTERN
#else
  #define EXTERN extern
#endif

#if defined MS_UNICODE
  #include <string>
  #include <vector>
  #include <map>
  using namespace std;
#endif /* MS_UNICODE */


/* Translation */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else /* !TECPLOTKERNEL */
inline TranslatedString_pa  _ts(const char *str,
                                const char *Notes = (const char*)NULL,
                                Widget      W     = (Widget)NULL)
  {
    return (TranslatedString_pa)str;
  }


inline const char*          _ts2char( const TranslatedString_pa str,
                                      const char *Notes = (const char*)NULL,
                                      Widget W = NULL)
  {
    return (const char *)str;
  }


inline TranslatedString_pa _nts( const char *str,
                                 const char *Notes = (const char*)NULL,
                                 Widget W = (Widget)NULL)
  {
    return (TranslatedString_pa)str;
  }
#endif /* !TECPLOTKERNEL */


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined MS_UNICODE

EXTERN LPCTSTR          _ts2tchar( const TranslatedString_pa str,
                                   const char *Notes = (const char*)NULL,
                                   Widget W = NULL);

EXTERN TranslatedString_pa _nts( LPCTSTR      str,
                                 const char   *Notes = (const char*)NULL,
                                 Widget       W = (Widget)NULL);

#else
  #define _ts2tchar _ts2char
#endif /* !MS_UNICODE */



#if defined MS_UNICODE
  #define VALID_UTF8_SBCHAR(c) Utf8Api::IsValidSingleByteUtf8Char(c)
  #define VALID_UTF8_MBCHAR(c) Utf8Api::IsValidMultiByteUtf8Char(c)
#else
  #define VALID_UTF8_SBCHAR(ch) ((ch) >= ' ' && (ch) <= (char)0x7F) 
  #define VALID_UTF8_MBCHAR(ch) FALSE
#endif

#define VALID_UTF8_CHAR(c) (VALID_UTF8_SBCHAR(c) || VALID_UTF8_MBCHAR(c))


#if defined UNIXX
  #define _T(str) (str)             /* NOP */
#endif

#if defined MS_UNICODE
  #define _WSTR_FROM_ANSI(str) Utf8Api::MultiByteToWideChar(str,CP_ACP)
  #define _WSTR(str) _ts2tchar(_nts(str))
  #define _CSTR(str) _ts2char(_nts(str))
#else
  #define _WSTR_FROM_ANSI(str)  (str)
  #define _WSTR(str)            (str)
  #define _CSTR(str)            (str)
#endif


#if defined MS_UNICODE
# include "W_UNICODE.h"
#endif

#endif /* Q_UNICODE_H_ */
