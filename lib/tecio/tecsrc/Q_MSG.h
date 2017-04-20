/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2004 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#if defined EXTERN
#undef EXTERN
#endif
#if defined Q_MSGMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

#define MAX_STATUS_LINE_MSG_LEN 255

EXTERN Boolean_t WrapString(const char  *OldString,
                            char       **NewString);
EXTERN void Warning(TranslatedString_pa Format,
                    ...); /* zero or more arguments */
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
EXTERN void ErrMsg(TranslatedString_pa Format,
                   ...); /* zero or more arguments */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#endif
#if !defined ENGINE
#if defined MOTIF
#endif
#endif
#if !defined ENGINE
#endif
#if defined Q_MSGMODULE
#else
#endif
#endif
