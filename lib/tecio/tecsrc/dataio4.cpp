#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2004 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#define DATAIO4MODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"

#include "AUXDATA.h"
#include "DATASET.h"
#include "FILESTREAM.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#include "GEOM2.h"
#include "GEOM.h"
#include "INPUT.h"
#include "SET.h"
#include "TEXT.h"
#include "DATAIO4.h"
#include "DATASET0.h"

#include "STRUTIL.h"
#include "ARRLIST.h"
#include "STRLIST.h"
#include "Q_MSG.h"

#if defined IRIS
  #include <ieeefp.h>
#endif

#if !defined(TECPLOTKERNEL) && defined(MSWIN)
# pragma warning(disable : 4244)
#endif

/*END HEADER*/

/*
 * This module contains mostly low level i/o functions.
 */

#if defined DECALPHA || defined COMPAQALPHA
  #define _IEEE_FP_INEXACT
  #define _IEEE_FP
#endif

#if defined SUN41
  #define SEEK_SET 0
  #define SEEK_CUR 1
  #define SEEK_END 2
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**********************************************************************
 **********************************************************************
 **********************        INPUT         **************************
 **********************************************************************
 **********************************************************************/

static char FilterFloatChar(float X)
{
  char C;
  if (((X >= 32.0) &&  (X <= 127.0)) || 
      ((X >= 160.0) && (X <= 255.0)) ||
      (X == 0.0))
    C = (char)X;
  else
    C = '?';
  return (C);
}


double GetNextValue(FileStream_s   *FileStream,
                    FieldDataType_e FieldDataType,
                    double          VMin,
                    double          VMax,
                    Boolean_t      *IsOk)
{
  double X = 0.0;
  
  REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));
  REQUIRE(!(*IsOk) || VALID_FIELD_DATA_TYPE(FieldDataType));
  REQUIRE(!(*IsOk) || VALID_REF(FileStream));
  
  if ( *IsOk )
    {
      switch (FieldDataType)
        {
          case FieldDataType_Float :
            {
              float XX;

              *IsOk = (FREAD(&XX,4,1,FileStream->File) == 1);

              if (!FileStream->IsByteOrderNative)
                REVERSE_4_BYTES(&XX);

              /* Convert from IEEE to Cray format */
              #if defined CRAY
                {
                  float X2;
                  int   i      = 1;
                  int   type   = 2;
                  int   bitoff = 0;
                  long  I;
                  char  *buf =  (char *)&XX;
                  IEG2CRAY(&type,&i,buf,&bitoff,&X2);
                  XX = X2;
                }
              #endif

              if (*IsOk)
                X = XX;
              else
                X = 0.0;
            } break;
          case FieldDataType_Double :
            {
              double XX;

              *IsOk = (FREAD(&XX,sizeof(double),1,FileStream->File) == 1);
              if (!FileStream->IsByteOrderNative)
                REVERSE_8_BYTES(&XX);

              /* Convert from IEEE to Cray format */
              #if defined CRAY
                {
                  int    i      = 1;
                  int    type   = 8;  
                  int    bitoff = 0;
                  double X2;
                  long   I;
                  char   *buf =  (char *)&XX;
                  IEG2CRAY(&type,&i,buf,&bitoff,&X2);
                  XX = X2;
                }
              #endif

              if (*IsOk)
                X = XX;
              else
                X = 0.0;
            } break;
          case FieldDataType_Int32  :
            {
              Int32_t L;
              *IsOk = (FREAD(&L,sizeof(Int32_t),1,FileStream->File) == 1);
              if (!FileStream->IsByteOrderNative)
                REVERSE_4_BYTES(&L);
              if (*IsOk)
                X = (double)L;
            } break;
          case FieldDataType_Int16  :
            {
              Int16_t S;
              *IsOk = (FREAD(&S,sizeof(Int16_t),1,FileStream->File) == 1);
              if (!FileStream->IsByteOrderNative)
                REVERSE_2_BYTES(&S);
              if (*IsOk)
                X = (double)S;
            } break;
          case FieldDataType_Byte  :
            {
              Byte_t B;
              *IsOk = (FREAD(&B,sizeof(Byte_t),1,FileStream->File) == 1);
              if (*IsOk)
                X = (double)B;
            } break;
          case FieldDataType_Bit :
            {
              /*
               * Important note:
               *   Reading bit data a value at a time is only valid for a
               *   single bit value. If the file contains a block of more than
               *   one bit value and you attempt to read it a bit at a time it
               *   will not work as Tecplot does not buffer the read. In order
               *   to read a block of bits you must perform a block read.
               */
              Byte_t B;
              *IsOk = (FREAD(&B,sizeof(Byte_t),1,FileStream->File) == 1);
              if (*IsOk)
                X = (double)(B & (Byte_t)01);
            } break;
          default: CHECK(FALSE); break;
        }

      if (*IsOk)
        {
          if ((X < VMin) || (X > VMax))
            {
              *IsOk = FALSE;
            }
        }
    }

  return X;
}


LgIndex_t GetNextI(FileStream_s *FileStream,
                   Boolean_t    *IsOk)
{
  LgIndex_t I = 0;

  REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));
  REQUIRE(!(*IsOk) || (VALID_REF(FileStream) && VALID_REF(FileStream->File)));

  if (*IsOk)
    {
      Int32_t Int32Val;
      *IsOk = (FREAD((void *)&Int32Val, 4, 1, FileStream->File) == 1);
      if (!FileStream->IsByteOrderNative)
        REVERSE_4_BYTES(&Int32Val);

      #if defined CRAY
        {
          int i      = 1;
          int type   = 1;
          int bitoff = 0;
          /* Note, buf may need to be 8 bytes??? */
          IEG2CRAY(&type,&i,(char *)&Int32Val,&bitoff,&I);
        }
      #else
        I = Int32Val;
      #endif
    }
  return I;
}


LgIndex_t GetIoFileInt(FileStream_s *FileStream,
                       short         Version,
                       LgIndex_t     IMin,
                       LgIndex_t     IMax,
                       Boolean_t    *IsOk)
{
  LgIndex_t I = 0;

  REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));
  REQUIRE(!(*IsOk) || (0 < Version && Version <= TecplotBinaryFileVersion));
  REQUIRE(!(*IsOk) || (VALID_REF(FileStream) && VALID_REF(FileStream->File)));
  REQUIRE(!(*IsOk) || IMin <= IMax);

  if (!(*IsOk))
    return (0);

  if (Version <= 63)
    {
      float  X;
      if (*IsOk)
        {
          X = (float)GetNextValue(FileStream,FieldDataType_Float,
                                  (double)IMin-1.0e-10,
                                  (double)IMax+1.0e-10,IsOk);
          if (*IsOk)
            {
              if (ABS(X) < (float)MAXINDEX)
                I = (LgIndex_t)X;
              else
                *IsOk = FALSE;
            }
          else
            *IsOk = FALSE;
        }
    }
  else
    {
      I = GetNextI(FileStream,IsOk);
    }

  if ((I < IMin) || (I > IMax))
    *IsOk = FALSE;

  return (I);
}

/**
 * Basically this is "realloc" but apparently "realloc" doesn't behave reliably
 * on all platforms.
 */
static Boolean_t ReallocString(char      **String,
                               LgIndex_t NewLength)
{
  Boolean_t IsOk;
  char *NewString;

  REQUIRE(VALID_REF(String));
  REQUIRE(*String == NULL || VALID_REF(*String));
  REQUIRE((*String != NULL && NewLength >= (LgIndex_t)strlen(*String)) ||
          (*String == NULL && NewLength >= 0));

  NewString = ALLOC_ARRAY(NewLength+1, char, "reallocated string");
  IsOk = (NewString != NULL);
  if (IsOk)
    {
      if (*String == NULL)
        {
          NewString[0] = '\0';
        }
      else
        {
          strcpy(NewString, *String);
          FREE_ARRAY(*String, "old string");
        }
      *String = NewString;
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  ENSURE(IMPLICATION(IsOk, VALID_REF(*String)));
  return IsOk;
}


/**
 * Reads a string from all versions of a binary plt file.
 *
 * param FileStream
 *     Open file stream positioned at the string to read.
 * param IVersion
 *     Binary file version number.
 * param MaxCharacters
 *     IVersion < 63
 *       This value is exactly the number of characters (actually floats, each
 *       one representing a character's ordinal value) to read from the file.
 *     IVersion >= 63 and ProcessData == TRUE
 *       If non-zero, this value represents the largest string to be returned.
 *       In other words, larger strings are read from the file but an allocated
 *       string of up to MaxCharacters is returned. A zero value indicates that
 *       the string size is unlimited and determined only by the actual length
 *       of the string in the file.
 * param TargetStr
 *     Pointer to hold the allocated string if ProcessData == TRUE.
 * param ProcessData
 *     Indicates if the read string should be retrieved.
 *
 * return
 *     TRUE if the read was successful with an allocated string *TargetStr
 *     containing the string read (no larger than requested or necessary),
 *     FALSE otherwise.
 */
Boolean_t ReadInString(FileStream_s  *FileStream,
                       short          IVersion,
                       int            MaxCharacters,
                       char         **TargetStr,
                       Boolean_t      ProcessData)
{
  Boolean_t IsOk = TRUE;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(0 < IVersion && IVersion <= TecplotBinaryFileVersion);
  REQUIRE(IMPLICATION(IVersion < 63 || ProcessData, MaxCharacters >= 0));
  REQUIRE(IMPLICATION(ProcessData, VALID_REF(TargetStr)));
  REQUIRE(VALID_BOOLEAN(ProcessData));

  if (IVersion < 63)
    {
      /*
       *  One word per character.  Read Exactly "MaxCharacters" number of words.
       */
      float X;

      if (ProcessData)
        {
          *TargetStr = ALLOC_ARRAY(MaxCharacters+1,char,"target string");
          IsOk = (*TargetStr != NULL);
        }

      if (IsOk)
        {
          LgIndex_t I;
          for (I = 0; IsOk && I < MaxCharacters; I++)
            {
              X = (float)GetNextValue(FileStream,FieldDataType_Float,0.0,127.0,&IsOk);
              if (!IsOk)
                break;
              if (ProcessData)
                (*TargetStr)[I] = FilterFloatChar(X);
            }
          if (ProcessData)
            (*TargetStr)[I] = '\0';
        }
      else
        {
          ErrMsg(_ts("Cannot allocate memory for string during read",
            "'string' meaning the computer science data type"));
        }
    }
  else
    {
      #define      MAX_STRBUFFER_LEN 4095
      static char StrBuffer[MAX_STRBUFFER_LEN+1];
      LgIndex_t   StrBufferLen = 0;
      LgIndex_t   TargetStrLen = 0;
      LgIndex_t   I = 0;
      LgIndex_t   CharValue = 0;

      if (ProcessData)
        *TargetStr = NULL;

      do
        {
          CharValue = GetIoFileInt(FileStream, IVersion, 0, 255, &IsOk);
          if (IsOk && ProcessData)
            {
              /* massage the character if necessary */
              if ((CharValue < 32 && CharValue != '\0' && CharValue != '\n') ||
                  (CharValue >= 128 && CharValue < 160))
                CharValue = ' ';

              /*
               * if the limit is not exceded, stuff the
               * character into the buffer
               */
              if (CharValue != '\0' &&
                  (I < MaxCharacters || MaxCharacters == 0))
                {
                  StrBuffer[StrBufferLen] = (char)CharValue;
                  StrBufferLen++;
                }

              if (CharValue == '\0' ||
                  StrBufferLen == MAX_STRBUFFER_LEN)
                {
                  if (StrBufferLen != 0 || *TargetStr == NULL)
                    {
                      StrBuffer[StrBufferLen] = '\0';
                      TargetStrLen += StrBufferLen;
                      IsOk = ReallocString(TargetStr, TargetStrLen);
                      if (IsOk)
                        strcat(*TargetStr, StrBuffer);
                      StrBufferLen = 0; /* reset the string buffer */
                    }
                }
            }

          I++;
        } while (IsOk && (char)CharValue != '\0');

      /* if we failed cleanup if necessary */
      if (!IsOk       &&
          ProcessData &&
          *TargetStr != NULL)
        {
          FREE_ARRAY(*TargetStr, "failed read string");
          *TargetStr = NULL;
        }
    }

  ENSURE(IMPLICATION(ProcessData,
                     (VALID_REF(*TargetStr) || *TargetStr == NULL)));
  ENSURE(VALID_BOOLEAN(IsOk));
  return (IsOk);
}

/**
 */
static void ReadDoubleBlock(FileStream_s *FileStream,
                            Boolean_t     DoRead,
                            double       *Buffer,
                            LgIndex_t     StartIndex,
                            LgIndex_t     NumValues,
                            Boolean_t    *IsOk)
{
  if (DoRead)
    {
      double *DPtr = Buffer + StartIndex;
      *IsOk = (FREAD(DPtr,sizeof(double),NumValues,FileStream->File) == (size_t)NumValues);
      if (!FileStream->IsByteOrderNative && *IsOk)
        {
          LgIndex_t N;
          for (N = 0; N < NumValues; N++)
            {
              REVERSE_8_BYTES(&DPtr[N]);
              #if defined CRAY
                {
                  int    i      = 1;
                  int    type   = 8;  
                  int    bitoff = 0;
                  double X2;
                  long   I;
                  char   *buf =  (char *)&DPtr[N];
                  IEG2CRAY(&type,&i,buf,&bitoff,&X2);
                  DPtr[N] = X2;
                }
              #endif
            }
        }
     }
   else
     *IsOk = (FSEEK(FileStream->File,NumValues*sizeof(double),SEEK_CUR) == 0);
}

/**
 */
static void ReadFloatBlock(FileStream_s *FileStream,
                           Boolean_t     DoRead,
                           float        *Buffer,
                           LgIndex_t     StartIndex,
                           LgIndex_t     NumValues,
                           Boolean_t    *IsOk)
{
  if (DoRead)
    {
      float *FPtr = Buffer + StartIndex;
      *IsOk = (FREAD(FPtr,sizeof(float),NumValues,FileStream->File) == (size_t)NumValues);
      if (!FileStream->IsByteOrderNative && *IsOk)
        {
          LgIndex_t N;
          for (N = 0; N < NumValues; N++)
            {
              REVERSE_4_BYTES(&FPtr[N]);
              #if defined CRAY
                {
                  float X2;
                  int   i      = 1;
                  int   type   = 2;
                  int   bitoff = 0;
                  long  I;
                  char  *buf =  (char *)&FPtr[N];
                  IEG2CRAY(&type,&i,buf,&bitoff,&X2);
                  FPtr[N] = X2;
                }
              #endif
            }
        }
    }
  else
    *IsOk = (FSEEK(FileStream->File,NumValues*sizeof(float),SEEK_CUR) == 0);
}


/**
 */
static void ReadBitBlock(FileStream_s *FileStream,
                         Boolean_t     DoRead,
                         Byte_t       *Buffer,
                         LgIndex_t     NumValues,
                         Boolean_t    *IsOk)
{
  /*
   * Do not allow reading of bit values if startindex is not 0.
   * (This means geometries cannot use bit data.
   */
  LgIndex_t NumBytes = (NumValues+7)/8;
  if (DoRead)
    {
      *IsOk = (FREAD(Buffer,
                     sizeof(Byte_t),
                     NumBytes,
                     FileStream->File) == (size_t)NumBytes);
    }
  else
    *IsOk = (FSEEK(FileStream->File,NumBytes*sizeof(Byte_t),SEEK_CUR) == 0);
}

/**
 */
static void ReadByteBlock(FileStream_s *FileStream,
                          Boolean_t     DoRead,
                          Byte_t       *Buffer,
                          LgIndex_t     StartIndex,
                          LgIndex_t     NumValues,
                          Boolean_t    *IsOk)
{
  if (DoRead)
    {
      *IsOk = (FREAD(Buffer+StartIndex,
                     sizeof(Byte_t),
                     NumValues,
                     FileStream->File) == (size_t)NumValues);
    }
  else
    *IsOk = (FSEEK(FileStream->File,NumValues*sizeof(Byte_t),SEEK_CUR) == 0);
}


/**
 */
static void ReadInt16Block(FileStream_s *FileStream,
                           Boolean_t     DoRead,
                           Int16_t      *Buffer,
                           LgIndex_t     StartIndex,
                           LgIndex_t     NumValues,
                           Boolean_t    *IsOk)
{
  if (DoRead)
    {
      Int16_t *IntPtr = Buffer+StartIndex;
      *IsOk = (FREAD(IntPtr,
                     sizeof(SmInteger_t),
                     NumValues,
                     FileStream->File) == (size_t)NumValues);

      if (!FileStream->IsByteOrderNative && *IsOk)
        {
          LgIndex_t N;
          for (N = 0; N < NumValues; N++)
            REVERSE_2_BYTES(&IntPtr[N]);
        }
    }
  else
    *IsOk = (FSEEK(FileStream->File,NumValues*sizeof(SmInteger_t),SEEK_CUR) == 0);
}

/**
 */
static void ReadInt32Block(FileStream_s *FileStream,
                           Boolean_t     DoRead,
                           Int32_t      *Buffer,
                           LgIndex_t     StartIndex,
                           LgIndex_t     NumValues,
                           Boolean_t    *IsOk)
{
  if (DoRead)
    {
      Int32_t *IntPtr = Buffer+StartIndex;
      *IsOk = (FREAD(IntPtr,
                     sizeof(LgIndex_t),
                     NumValues,
                     FileStream->File) == (size_t)NumValues);

      if (!FileStream->IsByteOrderNative && *IsOk)
        {
          LgIndex_t N;
          for (N = 0; N < NumValues; N++)
            REVERSE_4_BYTES(&IntPtr[N]);
        }
    }
  else
    *IsOk = (FSEEK(FileStream->File,NumValues*sizeof(LgIndex_t),SEEK_CUR) == 0);
}

/**
 */
void ReadPureBlock(FileStream_s   *FileStream,
                   Boolean_t       DoRead,
                   void           *Buffer,
                   FieldDataType_e FieldDataType,
                   LgIndex_t       StartIndex,
                   LgIndex_t       NumValues,
                   Boolean_t      *IsOk)
{
  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_BOOLEAN(DoRead));
  REQUIRE(!DoRead || VALID_REF(Buffer));
  REQUIRE(VALID_FIELD_DATA_TYPE(FieldDataType));
  REQUIRE(StartIndex >= 0);
  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));

  switch (FieldDataType)
    {
      case FieldDataType_Float :
        {
          ReadFloatBlock(FileStream,
                         DoRead,
                         (float *)Buffer,
                         StartIndex,
                         NumValues,
                         IsOk);
        } break;
      case FieldDataType_Double :
        {
          ReadDoubleBlock(FileStream,
                          DoRead,
                          (double *)Buffer,
                          StartIndex,
                          NumValues,
                          IsOk);
        } break;
      case FieldDataType_Bit :
        {
          if (StartIndex != 0)
            {
              ErrMsg(_ts("Internal Error: Attempt to read bit data at non-zero offset",
                "see Tecplot User's manual for a definition of 'bit' data"));
              *IsOk = FALSE;
            }
          else
            ReadBitBlock(FileStream,
                         DoRead,
                         (Byte_t *)Buffer,
                         NumValues,
                         IsOk);
        } break;
      case FieldDataType_Byte :
        {
          ReadByteBlock(FileStream,
                        DoRead,
                        (Byte_t *)Buffer,
                        StartIndex,
                        NumValues,
                        IsOk);
        } break;
      case FieldDataType_Int16 :
        {
          ReadInt16Block(FileStream,
                         DoRead,
                         (Int16_t *)Buffer,
                         StartIndex,
                         NumValues,
                         IsOk);
        } break;
      case FieldDataType_Int32 :
        {
          ReadInt32Block(FileStream,
                         DoRead,
                         (Int32_t *)Buffer,
                         StartIndex,
                         NumValues,
                         IsOk);
        } break;
      case FieldDataType_IJKFunction : /* Not used yet */
      case FieldDataType_Int64 :       /* Not used yet */
      default: CHECK(FALSE); break;
    }
  ENSURE(VALID_BOOLEAN(*IsOk));
}

/**
 */
void ReadBlock(FileStream_s   *FileStream,
               FieldData_pa    FieldData,
               Boolean_t       DoRead,
               FieldDataType_e FieldDataTypeInFile,
               LgIndex_t       StartIndex,
               LgIndex_t       EndIndex,
               Boolean_t      *IsOk)
{
  REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));
  REQUIRE(IMPLICATION(IsOk, VALID_REF(FileStream)));
  REQUIRE(IMPLICATION(IsOk, VALID_FIELD_DATA_TYPE(FieldDataTypeInFile)));
  REQUIRE(VALID_BOOLEAN(DoRead));
  REQUIRE(IMPLICATION(DoRead, VALID_REF(FieldData)));

  /*
   * Bit data is packed into bytes. Since Tecplot doesn't buffer reads it can
   * not perform bit by bit value reads and therefore must only perform block
   * reads of bit data.
   */
  Boolean_t ReadByBlock = IMPLICATION(DoRead, GetFieldDataType(FieldData) == FieldDataTypeInFile);
  REQUIRE(ReadByBlock || (FieldDataTypeInFile != FieldDataType_Bit));
  
  if ( *IsOk )
    {
      LgIndex_t NumValues = EndIndex-StartIndex+1;
      if ( ReadByBlock )
        {
          void *data_array;
          if ( DoRead )
            data_array = GetFieldDataVoidPtr(FieldData);
          else
            data_array = NULL;
          ReadPureBlock(FileStream,
                        DoRead,
                        data_array,
                        FieldDataTypeInFile,
                        StartIndex,
                        NumValues,
                        IsOk);
          /*TODO(CAM): Potentially remove this test as input data should not be reversed */
          if ( DoRead && !FieldData->IsNativeOrdering )
            {
              if ( FieldData->IsWordAligned )
                SwapBytesInTypedValueArray(GetFieldDataType(FieldData),
                                           GetFieldDataVoidPtr(FieldData),
                                           0,
                                           NumValues-1,
                                           1);
              else
                SwapBytesInUnalignedTypedValueArray(GetFieldDataType(FieldData),
                                                    GetFieldDataVoidPtr(FieldData),
                                                    0,
                                                    NumValues-1,
                                                    1);
            }
        }
      else
        {
          LgIndex_t N;
          for (N = 0; *IsOk && (N < NumValues); N++)
            {
              double D = GetNextValue(FileStream,FieldDataTypeInFile,-LARGEDOUBLE,LARGEDOUBLE,IsOk);
              if (DoRead)
                SetFieldValue(FieldData,N+StartIndex,D);
            }
        }
    }
}

/**
 */
void ReadClassicOrderedCCBlock(FileStream_s         *DataFileStream,
                               FieldData_pa          FieldData,
                               FieldDataType_e       FieldDataTypeInFile,
                               LgIndex_t             NumIPtsInFile,
                               LgIndex_t             NumJPtsInFile,
                               LgIndex_t             NumKPtsInFile,
                               Boolean_t            *IsOk)
{
  REQUIRE(IMPLICATION(*IsOk, VALID_REF(DataFileStream)));
  REQUIRE(IMPLICATION(*IsOk, VALID_FIELD_DATA_TYPE(FieldDataTypeInFile)));
  REQUIRE(VALID_REF(FieldData));
  REQUIRE(NumIPtsInFile >= 0);
  REQUIRE(NumJPtsInFile >= 0);
  REQUIRE(NumKPtsInFile >= 0);
  REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));
  
  if (*IsOk)
    {
      LgIndex_t J, K;
      LgIndex_t NumIJPts  = NumIPtsInFile * NumJPtsInFile;
      LgIndex_t IEnd      = MAX(NumIPtsInFile-1, 1);
      LgIndex_t JEnd      = MAX(NumJPtsInFile-1, 1);
      LgIndex_t KEnd      = MAX(NumKPtsInFile-1, 1);
      LgIndex_t NumValues = (IEnd * JEnd * KEnd);
      Boolean_t IsLinear  = ((NumJPtsInFile == 1 && NumKPtsInFile == 1) ||
                             (NumIPtsInFile == 1 && NumKPtsInFile == 1) ||
                             (NumIPtsInFile == 1 && NumJPtsInFile == 1));
      if (IsLinear)
        ReadBlock(DataFileStream, FieldData, TRUE, FieldDataTypeInFile,
                  0, NumValues-1, IsOk);
      else
        for (K = 0; K < KEnd && IsOk; K++)
        for (J = 0; J < JEnd && IsOk; J++)
          {
            LgIndex_t CellIndex = 0 + (J*NumIPtsInFile) + (K*NumIJPts);
            ReadBlock(DataFileStream, FieldData, TRUE, FieldDataTypeInFile,
                      CellIndex, CellIndex+IEnd-1, IsOk);
          }
    }

  ENSURE(VALID_BOOLEAN(*IsOk));
}

/**
 */
Boolean_t ReadConnectivityList(FileStream_s *FileStream,
                               short         IVersion,
                               NodeMap_pa    NodeMap,
                               LgIndex_t     TotValues)
{
  Boolean_t InputIsOk = TRUE;
  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(0 < IVersion && IVersion <= TecplotBinaryFileVersion);
  REQUIRE(VALID_REF(NodeMap));
  REQUIRE(TotValues > 0);

  #if defined CRAY
    {
      LgIndex_t ii;
      for (ii = 0; ii < TotValues && InputIsOk; ii++)
        {
          NodeMap_t NodeValue = GetNextI(FileStream,&InputIsOk);
          if (IVersion < 103)
            NodeValue--; /* versions less than 103 are ones based */
          SetNodeValueAtOffset(NodeMap,ii,NodeValue);
        }
    }
  #else
    {
      Int32_t *NMRawPtr = GetNMPtr(NodeMap);
      ReadInt32Block(FileStream,
                     TRUE,
                     NMRawPtr,
                     0,
                     TotValues,
                     &InputIsOk);
      if (InputIsOk)
        {
          if (IVersion < 103)
            {
              /* versions less than 103 are ones based */
              LgIndex_t II;
              for (II = 0; II < TotValues; II++)
                {
                  CHECK(NMRawPtr[II] >= 1);
                  NMRawPtr[II]--;                  
                }
            }
        }
    }
  #endif

  ENSURE(VALID_BOOLEAN(InputIsOk));
  return (InputIsOk);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 */
static void AdjustCustomColor(short         IVersion,
                              ColorIndex_t *BColor)
{
  REQUIRE(0 < IVersion && IVersion <= TecplotBinaryFileVersion);
  REQUIRE(VALID_REF(BColor));

  if ((IVersion < 70) && (*BColor >= 15) && (*BColor <= 22))
    *BColor -= 7;
}







Boolean_t ReadInDataFileTitleAndVarNames(FileStream_s  *FileStream,
                                         short          IVersion,
                                         char         **DataSetTitle,
                                         int           *NumVars,
                                         StringList_pa *VarNames)
{
  EntIndex_t   CurVar;
  Boolean_t    IsOk = TRUE;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(0 < IVersion && IVersion <= TecplotBinaryFileVersion);
  REQUIRE(VALID_REF(DataSetTitle) || (DataSetTitle == NULL));
  REQUIRE(VALID_REF(NumVars));
  REQUIRE(VALID_REF(VarNames) || (VarNames == NULL));

  *NumVars      = 0;
  if (DataSetTitle)
    *DataSetTitle = NULL;

  if (ReadInString(FileStream,
                   IVersion,
                   ((IVersion < 63) ? 80 : MaxChrsDatasetTitle),
                   DataSetTitle,
                   (Boolean_t)(DataSetTitle != NULL)))
    {
      if (DataSetTitle)
        TrimLeadAndTrailSpaces(*DataSetTitle);
      *NumVars = GetIoFileInt(FileStream,IVersion,0,MAXZONEMAP,&IsOk);
    }
  else
    IsOk = FALSE;

  if (IsOk && (*NumVars > MaxNumZonesOrVars))
    {
      ErrMsg(_ts("Too many variables"));
      IsOk = FALSE;
    }

  if (IsOk && VarNames)
    {
      if (*NumVars > 0)
        {
          /* allocate a string list filled with NULL's */
          *VarNames = StringListAlloc();
          IsOk = (*VarNames != NULL);
          if (IsOk)
            IsOk = StringListSetString(*VarNames, *NumVars-1, NULL);

          if (!IsOk)
            {
              if (*VarNames != NULL)
                StringListDealloc(VarNames);
              ErrMsg(_ts("Out of space while allocating var names"));
            }
        }
    }

  for (CurVar = 0; IsOk && (CurVar < *NumVars); CurVar++)
    {
      char *VName = NULL;

      IsOk = ReadInString(FileStream,
                          IVersion,
                          ((IVersion < 63) ? 5 : MaxChrsVarName),
                          VarNames ? &VName : NULL,
                          (Boolean_t)(VarNames != NULL));
      if (IsOk && VarNames)
        {
          if (VName == NULL)
            {
              /* NULL variable names are converted to empty names */
              VName = ALLOC_ARRAY(1, char, "empty variable name");
              strcpy(VName, "");
            }
          TrimLeadAndTrailSpaces(VName);

          /*
           * variables are not allowed to have litteral new line characters
           * within them but they can sneek in from ASCII data files so
           * convert any to their appropriate two character representation
           */
          IsOk = ReplaceNewlineWithBackslashN(&VName);

          IsOk = IsOk && StringListSetString(*VarNames, CurVar, VName);
          if (VName != NULL)
            FREE_ARRAY(VName, "variable name");
        }

      if (!IsOk)
        {
          if (VarNames && *VarNames)
            StringListDealloc(VarNames);
          ErrMsg(_ts("Out of space while allocating variable names"));
        }
    }
  ENSURE(VALID_BOOLEAN(IsOk));
  return (IsOk);
}








/**
 */
static Boolean_t ReadInPresetZoneColor(FileStream_s *FileStream,
                                       short         IVersion,
                                       ZoneSpec_s   *ZoneSpec)
{
  Boolean_t IsOk = TRUE;
  LgIndex_t ZoneColor;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(0 < IVersion && IVersion <= TecplotBinaryFileVersion);
  REQUIRE(VALID_REF(ZoneSpec));

  ZoneColor = GetIoFileInt(FileStream,IVersion,-1,LastBasicColor,&IsOk);
  if (IsOk)
    {
      if ( VALID_BASIC_COLOR(ZoneColor) )
        {
          ZoneSpec->ZoneLoadInfo.PresetZoneColor = (EntIndex_t)ZoneColor;
          AdjustCustomColor(IVersion,&ZoneSpec->ZoneLoadInfo.PresetZoneColor);
        }
      else if (ZoneColor != -1)
        IsOk = FALSE;
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}

/**
 */
static void ConvertCommonTimeToSolutionTime(ZoneSpec_s *ZoneSpec)
{
  REQUIRE(VALID_REF(ZoneSpec));
  REQUIRE(ZoneSpec->AuxData == NULL || VALID_REF(ZoneSpec->AuxData));

  LgIndex_t ItemIndex;
  if (ZoneSpec->AuxData != NULL &&
      AuxDataGetItemIndex(ZoneSpec->AuxData, AuxData_Common_Time, &ItemIndex))
    {
      const char    *SameName;
      ArbParam_t     Value;
      AuxDataType_e  Type;
      Boolean_t      Retain;

      AuxDataGetItemByIndex(ZoneSpec->AuxData, ItemIndex,
                            &SameName, &Value, &Type, &Retain);
      CHECK(ustrcmp(AuxData_Common_Time, SameName) == 0);
      CHECK(Type == AuxDataType_String);

      char *EndPtr = NULL;
      double SolutionTime = strtod((const char *)Value, &EndPtr);
      if (EndPtr != (char *)Value)
        {
          /* we only allow white space to trail a value */
          while (isspace(*EndPtr))
            EndPtr++;
        }
      if (EndPtr != (char *)Value && *EndPtr == '\0')
        {
          ZoneSpec->SolutionTime = SolutionTime;
          ZoneSpec->StrandID     = STRAND_ID_PENDING;
          AuxDataDeleteItemByIndex(ZoneSpec->AuxData, ItemIndex);
        }
    }
}

/*
 * Pass1 for a zone reads in and initializes the zone structures.
 * These structures are released later if the user elects to not read them
 * in. 
 */
Boolean_t ReadInZoneHeader(FileStream_s *FileStream,
                           short         IVersion,
                           ZoneSpec_s   *ZoneSpec,
                           Set_pa        IsVarCellCentered,
                           EntIndex_t    NumVars,
                           Boolean_t    *IsRawFNAvailable)
{
  EntIndex_t Var;
  Boolean_t  IsOk = TRUE;
  LgIndex_t  I1;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(0 < IVersion && IVersion <= TecplotBinaryFileVersion);
  REQUIRE(VALID_REF(ZoneSpec));
  REQUIRE(IsVarCellCentered == NULL || VALID_REF(IsVarCellCentered));
  REQUIRE(NumVars >= 0);
  REQUIRE(VALID_REF(IsRawFNAvailable));

  SetZoneSpecDefaults(ZoneSpec);

  if (IsVarCellCentered != NULL)
    {
      /* assign default variable value location: nodal */
      ClearSet(IsVarCellCentered);
      IsOk = ExpandSet(IsVarCellCentered, NumVars, FALSE);
    }

  if (IsOk)
    IsOk = ReadInString(FileStream,IVersion,
                        ((IVersion < 63) ? 10 : MaxChrsZnTitle),
                        &ZoneSpec->Name,
                        TRUE);

  if (IsOk && ZoneSpec->Name == NULL)
    {
      /* NULL zone names are converted to empty names */
      ZoneSpec->Name = ALLOC_ARRAY(1, char, "empty zone name");
      IsOk = (ZoneSpec->Name != NULL);
      if (IsOk)
        strcpy(ZoneSpec->Name, "");
    }

  if (IsOk)
    TrimLeadAndTrailSpaces(ZoneSpec->Name);

  if (IVersion < 101)
    {
      Boolean_t    IsZoneFinite;
      DataFormat_e ZoneDataFormat;

      I1 = GetIoFileInt(FileStream,IVersion,0,3,&IsOk);

      if ((I1 < 0) || (I1 > 3))
        {
          return (FALSE);
        }

      ZoneDataFormat = (DataFormat_e)I1;

      IsZoneFinite = (ZoneDataFormat == DataFormat_FEPoint ||
                      ZoneDataFormat == DataFormat_FEBlock);

      ZoneSpec->ZoneLoadInfo.IsInBlockFormat = (ZoneDataFormat == DataFormat_IJKBlock ||
                                                ZoneDataFormat == DataFormat_FEBlock);

      if (IVersion > 62)
        IsOk = ReadInPresetZoneColor(FileStream, IVersion, ZoneSpec);

      if (IVersion < 60)
        GetNextValue(FileStream,FieldDataType_Float,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);  /* Old ZPlane Value */

      if (IsOk)
        {
          ZoneSpec->NumPtsI = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          ZoneSpec->NumPtsJ = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          if (IVersion >= 60)
            ZoneSpec->NumPtsK = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          else
            ZoneSpec->NumPtsK = 1;
        }

      if (IsOk)
        {
          /* If IMax,JMax, & KMax are all zero then this zone was "zombied" by
             a partial read to make layout files align. */

          if (!((ZoneSpec->NumPtsI == 0) && 
                (ZoneSpec->NumPtsJ == 0) && 
                (ZoneSpec->NumPtsK == 0)) &&
              ((ZoneSpec->NumPtsI <= 0) || 
               (ZoneSpec->NumPtsJ <= 0) || 
               (ZoneSpec->NumPtsK < 0)  || 
               ((!IsZoneFinite && (ZoneSpec->NumPtsK == 0)))))
            {
              ErrMsg(_ts("Datafile is corrupted"));
              IsOk = FALSE;
            }

          if (IsZoneFinite)
            {
              if (IVersion >= 61)
                {
                  ZoneSpec->Type = (ZoneType_e)(ZoneSpec->NumPtsK+1);
                  switch (ZoneSpec->Type)
                    {
                      case ZoneType_FETriangle : 
                        {
                          ZoneSpec->NumPtsK = 3;
                        } break;
                      case ZoneType_FEQuad :
                        {
                          ZoneSpec->NumPtsK = 4;
                        } break;
                      case ZoneType_FETetra :
                        {
                          ZoneSpec->NumPtsK = 4;
                        } break;
                      case ZoneType_FEBrick :     
                        {
                          ZoneSpec->NumPtsK = 8;
                        } break;
                      case ZoneType_FELineSeg :     
                        {
                          ZoneSpec->NumPtsK = 2;
                        } break;
                      default :
                        {
                          ErrMsg(_ts("Datafile corrupted: Invalid element type for FE DataSet"));
                          IsOk = FALSE;
                        }
                    }
                }
              else
                {
                  ZoneSpec->Type = ZoneType_FEQuad;
                  ZoneSpec->NumPtsK = 4;
                }
            }
          else 
            {
              ZoneSpec->Type = ZoneType_Ordered;

              ZoneSpec->ICellDim = ZoneSpec->NumPtsI - 1;
              ZoneSpec->JCellDim = ZoneSpec->NumPtsJ - 1;
              ZoneSpec->KCellDim = ZoneSpec->NumPtsK - 1;
            }
        }
    }
  else
    {
      if (IVersion >= 107)
        {
          ZoneSpec->ParentZone = GetIoFileInt(FileStream,IVersion,-1,MAXZONEMAP-1,&IsOk);
          if (!IsOk)
            ErrMsg(_ts("Invalid datafile: parent zone assignment must be to an existing zone within the same datafile."));
        }

      if (IVersion >= 106)
        {
          /* Strand ID and solution time. Strand ID's of STRAND_ID_PENDING, -2, instruct Tecplot to generate strand ID's */
          ZoneSpec->StrandID     = GetIoFileInt(FileStream,IVersion,-2,MAXZONEMAP-1,&IsOk);
          ZoneSpec->SolutionTime = GetNextValue(FileStream,FieldDataType_Double,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
          if (!IsOk)
            ErrMsg(_ts("Invalid datafile: bad StrandID or SolutionTime"));
        }

      /* preset zone color */
      IsOk = IsOk && ReadInPresetZoneColor(FileStream, IVersion, ZoneSpec);

      /* ZoneType */
      I1 = (ZoneType_e)GetIoFileInt(FileStream,IVersion,0,5,&IsOk);
      switch (I1)
        {
          case 0: ZoneSpec->Type = ZoneType_Ordered; break;
          case 1: ZoneSpec->Type = ZoneType_FELineSeg; break;
          case 2: ZoneSpec->Type = ZoneType_FETriangle; break;
          case 3: ZoneSpec->Type = ZoneType_FEQuad; break;
          case 4: ZoneSpec->Type = ZoneType_FETetra; break;
          case 5: ZoneSpec->Type = ZoneType_FEBrick; break;
          default: ErrMsg(_ts("Invalid datafile: unknown zone type.")); IsOk = FALSE; break;
        }

      /* DataPacking */
      ZoneSpec->ZoneLoadInfo.IsInBlockFormat =
          ((DataPacking_e)GetIoFileInt(FileStream,IVersion,0,1,&IsOk) == DataPacking_Block);

      /* is the variable value location specified? */
      if ((Boolean_t)GetIoFileInt(FileStream,IVersion,0,1,&IsOk) && IsOk)
        {
          /* Variable Value Location foreach Var */
          for (Var = 0; Var < NumVars && IsOk; Var++)
            {
              if ((Boolean_t)GetIoFileInt(FileStream,IVersion,0,1,&IsOk) && IsOk &&
                  IsVarCellCentered != NULL)
                {
                  IsOk = (ZoneSpec->ZoneLoadInfo.IsInBlockFormat);
                  if (IsOk)
                    IsOk = AddToSet(IsVarCellCentered, Var, FALSE);
                  else
                    ErrMsg(_ts("Invalid datafile: cell centered "
                           "variable must be in block format.",
                           "See the Tecplot User's Manual for a definition of 'block format'"));
                }
            }
        }

      /* are raw face neighbors supplied in the data section? */
      if (IVersion == 108 && IsOk)
        {
          *IsRawFNAvailable = GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
          if (*IsRawFNAvailable &&
              (ZoneSpec->Type == ZoneType_Ordered ||
               ZoneSpec->Type == ZoneType_FELineSeg))
            {
              IsOk = FALSE;
              ErrMsg(_ts("Invalid datafile: raw face neighbors may not be "
                         "supplied for ordered or FE line segment zones."));
            }
        }
      else
        *IsRawFNAvailable = FALSE;

      /* miscellaneous face neighbor info */
      ZoneSpec->UserDefFNConnections = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
      if (ZoneSpec->UserDefFNConnections != 0)
        ZoneSpec->FNMode = (FaceNeighborMode_e)GetIoFileInt(FileStream,IVersion,0,3,&IsOk);

      if (IVersion == 108 && IsOk)
        {
          Boolean_t FaceNeighborsComplete = FALSE;
          if (ZoneSpec->UserDefFNConnections != 0 &&
              ZoneSpec->Type != ZoneType_Ordered)
            FaceNeighborsComplete = (Boolean_t)GetIoFileInt(FileStream,IVersion,0,1,&IsOk);

          /*
           * If miscellaneous face neighbors completely specify the set then we
           * don't want Tecplot to auto-assign the face neighbors.
           */
          ZoneSpec->AutoAssignFN = !FaceNeighborsComplete;
        }

      if (ZoneSpec->Type == ZoneType_Ordered)
        {
          /* IMax, JMax, KMax */
          ZoneSpec->NumPtsI = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          ZoneSpec->NumPtsJ = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          ZoneSpec->NumPtsK = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          /*
           * if not a zombie zone (zombie zone: points in all dimensions are
           * zero) then points in each direction must be specified
           */
          if (IsOk &&
              !(ZoneSpec->NumPtsI == 0 && 
                ZoneSpec->NumPtsJ == 0 && 
                ZoneSpec->NumPtsK == 0)  &&
              (ZoneSpec->NumPtsI == 0 || 
               ZoneSpec->NumPtsJ == 0 || 
               ZoneSpec->NumPtsK == 0))
            {
              ErrMsg(_ts("Invalid data file: incorrect specification of "
                     "I, J, or K points for ordered data set."));
              IsOk = FALSE;
            }
        }
      else
        {
          ZoneSpec->NumPtsI = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          ZoneSpec->NumPtsJ = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          switch (ZoneSpec->Type)
            {
              case ZoneType_FETriangle : 
                {
                  ZoneSpec->NumPtsK = 3;
                } break;
              case ZoneType_FEQuad :
                {
                  ZoneSpec->NumPtsK = 4;
                } break;
              case ZoneType_FETetra :
                {
                  ZoneSpec->NumPtsK = 4;
                } break;
              case ZoneType_FEBrick :     
                {
                  ZoneSpec->NumPtsK = 8;
                } break;
              case ZoneType_FELineSeg :     
                {
                  ZoneSpec->NumPtsK = 2;
                } break;
              default :
                {
                  ErrMsg(_ts("Invalid data file: invalid element type for FE data set."));
                  IsOk = FALSE;
                }
            }

          ZoneSpec->ICellDim = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          ZoneSpec->JCellDim = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
          ZoneSpec->KCellDim = GetIoFileInt(FileStream,IVersion,0,MAXINDEX,&IsOk);
        }

      /* Zone Auxiliary Data indicator followed by Zone Auxiliary Data */
      for (I1 = GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
           IsOk && I1 != 0;
           I1 = GetIoFileInt(FileStream,IVersion,0,1,&IsOk))
        {
          if (ZoneSpec->AuxData == NULL)
            ZoneSpec->AuxData = AuxDataAlloc();
          IsOk = (ZoneSpec->AuxData != NULL);
          if (IsOk)
            IsOk = ReadInAuxData(FileStream, IVersion, ZoneSpec->AuxData);
        }
    }

  /*
   * Convert AuxZone's Common.Time from non-time aware data files to zone
   * solution time if it exists.
   */
  if (IVersion < 106 && IsOk)
    ConvertCommonTimeToSolutionTime(ZoneSpec);

  ENSURE(VALID_BOOLEAN(IsOk));
  return (IsOk);
}




/*
 *  Pass1 for Custom labels simply acknowledges that a custom label was 
 *  parsed and skips over the labels.  
 */

Boolean_t ReadInCustomLabels(FileStream_s  *FileStream,
                             short          IVersion,
                             Boolean_t      OkToLoad,
                             StringList_pa *CustomLabelBase)
{
  LgIndex_t NumLabels;
  short     I;
  Boolean_t IsOk = TRUE;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(IVersion > 0);
  REQUIRE(VALID_BOOLEAN(OkToLoad));
  REQUIRE(!(OkToLoad) || VALID_REF(CustomLabelBase));

  NumLabels = (short)GetIoFileInt(FileStream,IVersion,1,MAXINDEX,&IsOk);
  if (IsOk && NumLabels != 0 && OkToLoad)
    {
      *CustomLabelBase = StringListAlloc();
      IsOk = (*CustomLabelBase != NULL);
      if (!IsOk)
        ErrMsg(_ts("Cannot allocate memory for Custom Labels."));
    }

  for (I = 0; IsOk && (I < NumLabels); I++)
    {
      char *TLabel = NULL;

      IsOk = ReadInString(FileStream,IVersion,
                          49,
                          &TLabel,
                          OkToLoad);
      TrimLeadAndTrailSpaces(TLabel);

      if (IsOk && OkToLoad)
        {
          #if defined TECPLOTKERNEL
            if (GeneralBase.Interface.PrintDebug)
              fprintf(StdError,"Label %d: %s\n",I+1,TLabel);
          #endif
          IsOk = StringListAppendString(*CustomLabelBase, TLabel);
          if (TLabel != NULL)
            FREE_ARRAY(TLabel, "custom label");
          if (!IsOk)
            ErrMsg(_ts("Cannot allocate memory for Custom Label."));
        }
    }
  if (!IsOk)
    ErrMsg(_ts("Invalid custom axis label record in binary datafile"));

  ENSURE(VALID_BOOLEAN(IsOk));
  ENSURE(!(IsOk && NumLabels != 0 && OkToLoad) ||
         StringListValid(*CustomLabelBase));
  return IsOk;
}


Boolean_t ReadInUserRec(FileStream_s  *FileStream,
                        short          IVersion,
                        int            MaxCharactersAllowed,
                        char         **UserRec) /* NULL if to ignore */
{
  if (!ReadInString(FileStream,IVersion,
                    MaxCharactersAllowed,
                    UserRec,
                    (Boolean_t)(UserRec != NULL)))
    {
      ErrMsg(_ts("Invalid USERREC record in binary datafile"));
      return (FALSE);
    }
  return (TRUE);
}


/**
 */
Boolean_t ReadInAuxData(FileStream_s *FileStream,
                        short         IVersion,
                        AuxData_pa    AuxData)
{
  Boolean_t IsOk;
  Boolean_t DoCollectData;
  char      *AuxName  = NULL;
  LgIndex_t AuxValueType;
  char      *AuxValue = NULL;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(0 < IVersion && IVersion <= TecplotBinaryFileVersion);
  REQUIRE(AuxData == NULL || VALID_REF(AuxData));

  DoCollectData = (AuxData != NULL);
  IsOk = ReadInString(FileStream,
                      IVersion,
                      MaxChrsVarName, /* ... seems reasonable */
                      &AuxName,
                      DoCollectData);
  if (IsOk && DoCollectData && !AuxDataIsValidName(AuxName))
    {
      ErrMsg(_ts("Invalid auxiliary data name."));
      IsOk = FALSE;
    }

  /*
   * currently only one value type is supported
   *   0: AuxiliaryValueFormat_String
   */
  if (IsOk)
    {
      AuxValueType = GetIoFileInt(FileStream, IVersion, 0, 0, &IsOk);
      if (IsOk && (AuxValueType != (LgIndex_t)AuxDataType_String))
        {
          ErrMsg(_ts("Unsupported auxiliary data type"));
          IsOk = FALSE;
        }
    }

  if (IsOk)
    IsOk = ReadInString(FileStream,
                        IVersion,
                        MaxChrsAuxValueString,
                        &AuxValue,
                        DoCollectData);
  if (IsOk && DoCollectData)
    IsOk = AuxDataSetItem(AuxData,
                          AuxName, (ArbParam_t)AuxValue,
                          AuxDataType_String,
                          TRUE); /* Retain */

  /* cleanup: auxiliary data made a copy of the name and value */
  if (AuxName != NULL)
    FREE_ARRAY(AuxName, "data set auxiliary data item name");
  if (AuxValue != NULL)
    FREE_ARRAY(AuxValue, "data set auxiliary data item value");

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


static void GetZoneAttachment(FileStream_s *FileStream,
                              short         IVersion,
                              EntIndex_t   *Z,
                              Boolean_t    *IsAttached,
                              Boolean_t    *IsOk)
{
  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(Z));
  REQUIRE(VALID_REF(IsAttached));
  REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));

  if (IVersion >= 47)
    *Z = (EntIndex_t)GetIoFileInt(FileStream,IVersion,-1,MAXZONEMAP,IsOk);
  else
    *Z = 0;

  if (IVersion < 70)
    (*Z)--;

  if (*Z == -1)
    {
      *Z          = 0;
      *IsAttached = FALSE;
    }
  else
    *IsAttached = TRUE;

  ENSURE(VALID_BOOLEAN(*IsAttached));
  ENSURE(VALID_BOOLEAN(*IsOk));
  ENSURE(*Z >= 0);
}


static Boolean_t ReadMacroFunctionCommand(FileStream_s  *FileStream,
                                          short          IVersion,
                                          Boolean_t      OkToLoad,
                                          char         **MacroFunctionCommand)
{
  Boolean_t Result = FALSE;
  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(IVersion > 0);
  REQUIRE(VALID_BOOLEAN(OkToLoad));
  REQUIRE(VALID_REF(MacroFunctionCommand));

  Result = ReadInString(FileStream,IVersion,MAXASCIICOMMANDLENGTH,MacroFunctionCommand,OkToLoad);

  ENSURE(VALID_BOOLEAN(Result));
  return (Result);
}


/*
 *  Pass1 for Geometries simply acknowledges that a geometry was
 *  parsed and skips over the geometry.
 */
Boolean_t ReadInGeometry(FileStream_s *FileStream,
                         short         IVersion,
                         Boolean_t     OkToLoad,
                         Geom_s       *Geom,
                         LgIndex_t     MaxDataPts)
{
  LgIndex_t        I;
  LgIndex_t        S;
  FieldDataType_e  FFT;
  Boolean_t        IsOk = TRUE;
  const char      *ErrMsgString = _ts2char(_ts("Invalid geometry record"));

  REQUIRE(VALID_REF(Geom));

  if (IVersion < 70)
    FFT = FieldDataType_Float;
  else
    FFT = FieldDataType_Double;

  if (IVersion < 101)
    I = GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
  else
    I = GetIoFileInt(FileStream,IVersion,0,4,&IsOk);

  if ( I == 0 )
    Geom->PositionCoordSys = CoordSys_Grid;
  else if ( I == 1 )
    Geom->PositionCoordSys = CoordSys_Frame;
  /*
   * I == 2 is for CoordSys_FrameOffset and is not used currently
   *
   * I == 3 is for the old window coordinate system
   */
  else if ( I == 4 )
    Geom->PositionCoordSys = CoordSys_Grid3D;
  else
    {
      ErrMsgString = _ts2char(_ts("Invalid geometry coordinate system"));
      IsOk = FALSE;
    }

  Geom->Scope = (Scope_e)GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
  if (IVersion >= 102)
    Geom->DrawOrder = (DrawOrder_e)GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
  Geom->AnchorPos.Generic.V1 = GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  Geom->AnchorPos.Generic.V2 = GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  if (IVersion >= 45)
    Geom->AnchorPos.Generic.V3 = GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  else
    Geom->AnchorPos.Generic.V3 = 0.0;

  GetZoneAttachment(FileStream,IVersion,&Geom->Zone,&Geom->AttachToZone,&IsOk);

  Geom->BColor = (SmInteger_t)GetIoFileInt(FileStream,IVersion,0,255,&IsOk);

  AdjustCustomColor(IVersion,&Geom->BColor);

  if (IVersion > 47)
    {
      Geom->FillBColor = (SmInteger_t)GetIoFileInt(FileStream,IVersion,0,255,&IsOk);
      Geom->IsFilled  = (Boolean_t)GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
      AdjustCustomColor(IVersion,&Geom->FillBColor);
    }
  else
    {
      Geom->FillBColor = Geom->BColor;
      Geom->IsFilled  = FALSE;
    }

  if (IVersion < 101)
    {
      Geom->GeomType = (GeomType_e)GetIoFileInt(FileStream,IVersion,0,5,&IsOk);
      if (Geom->GeomType == GeomType_LineSegs3D)
        {
          /*
           * GeomType_LineSegs3D is deprecated, converter to GeomType_LineSegs
           * with CoordSys_Grid3D instead
           */
          Geom->GeomType         = GeomType_LineSegs;
          Geom->PositionCoordSys = CoordSys_Grid3D; /*...should have been anyway */
        }
    }
  else
    {
      Geom->GeomType = (GeomType_e)GetIoFileInt(FileStream,IVersion,0,4,&IsOk);
    }

  /*
   * Check geom coord sys versus geom type
   */
  if ( Geom->PositionCoordSys == CoordSys_Grid3D &&
       Geom->GeomType != GeomType_LineSegs )
    {
      ErrMsgString = _ts2char(_ts("Mismatch between geometry coordinate system and geometry type"));
      IsOk = FALSE;
    }

  if (IVersion > 41)
    {
      Geom->LinePattern  = (LinePattern_e)GetIoFileInt(FileStream,IVersion,0,(LgIndex_t)LinePattern_DashDotDot,&IsOk);
    }
  else
    {
      Geom->LinePattern  = (LinePattern_e)((int)Geom->GeomType % 2);
      Geom->GeomType     = (GeomType_e)((int)Geom->GeomType / 10);
    }

  if ((IVersion < 49) && ((short)(Geom->GeomType) == 2))
    {
      Geom->GeomType = GeomType_Rectangle;
      Geom->IsFilled = TRUE;
    }

  if ((IVersion < 70) && ((short)(Geom->GeomType) > 1))
    Geom->GeomType = (GeomType_e)((short)Geom->GeomType + 1);

  ResetString(&Geom->MacroFunctionCommand,NULL,TRUE);

  Geom->ImageResizeFilter = ImageResizeFilter_Texture;

  if (IVersion >= 70)
    {
      Geom->PatternLength       = GetNextValue(FileStream,FFT,
                                            PatternLengthInputSpec.Min,
                                            PatternLengthInputSpec.Max,
                                            &IsOk);
      Geom->LineThickness       = GetNextValue(FileStream,FFT,
                                            LineThicknessInputSpec.Min,
                                            LineThicknessInputSpec.Max,
                                            &IsOk);
      Geom->NumEllipsePts       = (SmInteger_t)GetIoFileInt(FileStream,IVersion,2,MaxPtsCircleOrEllipse,&IsOk);
      Geom->ArrowheadStyle      = (ArrowheadStyle_e)GetIoFileInt(FileStream,IVersion,0,(LgIndex_t)Arrowhead_Hollow,&IsOk);
      Geom->ArrowheadAttachment = (ArrowheadAttachment_e)GetIoFileInt(FileStream,IVersion,
                                                                      0,
                                                                      (LgIndex_t)ArrowheadAttach_AtBothEnds,
                                                                      &IsOk);

      Geom->ArrowheadSize  = GetNextValue(FileStream,FFT,
                                       ArrowheadSizeInputSpec.Min,
                                       ArrowheadSizeInputSpec.Max,
                                       &IsOk);
      Geom->ArrowheadAngle = GetNextValue(FileStream,FFT,
                                       ArrowheadAngleInputSpec.Min,
                                       ArrowheadAngleInputSpec.Max,
                                       &IsOk);

      if (IVersion >= 75)
        {
          IsOk = ReadMacroFunctionCommand(FileStream,
                                          IVersion,
                                          OkToLoad,
                                          &Geom->MacroFunctionCommand);
        } /* version >= 75 */
    } /* version >= 70 */
  else
    {
      Geom->LineThickness        = 0.001;
      Geom->PatternLength        = 0.02;
      Geom->ArrowheadStyle       = Arrowhead_Plain;
      Geom->ArrowheadAttachment  = ArrowheadAttach_None;
      Geom->ArrowheadSize        = 0.05;
      Geom->ArrowheadAngle       = 12.0/DEGPERRADIANS;
    }

  if (IVersion < 41)
    {
      GetNextValue(FileStream,FieldDataType_Float,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
      GetNextValue(FileStream,FieldDataType_Float,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
      GetNextValue(FileStream,FieldDataType_Float,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
    }

  if (IVersion < 70)
    Geom->DataType = FieldDataType_Float;
  else
    Geom->DataType = (FieldDataType_e)GetIoFileInt(FileStream,IVersion,1,2,&IsOk);
  CHECK(VALID_GEOM_FIELD_DATA_TYPE(Geom->DataType));

  Geom->Clipping = Clipping_ClipToViewport; /* default value for pre 101 versions */
  if (IVersion >= 101)
    {
      Geom->Clipping = (Clipping_e)GetIoFileInt(FileStream,IVersion,0,2,&IsOk);
      /*
       * The second clipping value was deprecated during v10 development and thus removed.
       * This moved Clipping_ClipToFrame to the 2nd position from the 3rd, so we convert
       * value 2 to ClipToFrame to support files made during v10 developement.
       */
      if ( Geom->Clipping == (Clipping_e)2 )
        Geom->Clipping = Clipping_ClipToFrame;
    }
   
  if (IVersion < 50 || 
      Geom->GeomType == GeomType_LineSegs)
    {
      Geom->NumSegments = (SmInteger_t)GetIoFileInt(FileStream,IVersion,1,MaxGeoSegments,&IsOk);
      #if defined TECPLOTKERNEL
        if (GeneralBase.Interface.PrintDebug && OkToLoad)
          fprintf(StdError,"GEOMETRY: NumSegs = %d\n",Geom->NumSegments);
      #endif
      S = -1;
      I = 0;
      while ((S+1 < Geom->NumSegments) && 
             !feof(FileStream->File) && IsOk)
        {
          S++;
          Geom->NumSegPts[S] = GetIoFileInt(FileStream,IVersion,1,MAXINDEX,&IsOk);
          if ((I+Geom->NumSegPts[S] > MaxDataPts) && OkToLoad)
            {
              ErrMsgString = _ts2char(_ts("Geometry is too big"));
              IsOk = FALSE;
            }
          else
            {
              ReadBlock(FileStream,Geom->GeomData.Generic.V1Base,OkToLoad,Geom->DataType,I,I+Geom->NumSegPts[S]-1,&IsOk);
              ReadBlock(FileStream,Geom->GeomData.Generic.V2Base,OkToLoad,Geom->DataType,I,I+Geom->NumSegPts[S]-1,&IsOk);
              if (Geom->PositionCoordSys == CoordSys_Grid3D)
                ReadBlock(FileStream,Geom->GeomData.Generic.V3Base,OkToLoad,Geom->DataType,I,I+Geom->NumSegPts[S]-1,&IsOk);
              I += Geom->NumSegPts[S];
            }
        }
      if (IsOk && (Geom->GeomType == GeomType_Rectangle))     /* IVersion < 50 */
        {
          if (OkToLoad)
            {
              CopyFieldValue(Geom->GeomData.Generic.V1Base,0,Geom->GeomData.Generic.V1Base,2);
              CopyFieldValue(Geom->GeomData.Generic.V2Base,0,Geom->GeomData.Generic.V2Base,2);
            }
        }
    }
  else if (Geom->GeomType == GeomType_Rectangle ||
           Geom->GeomType == GeomType_Ellipse)
    {
      double XX,YY;
      XX = GetNextValue(FileStream,Geom->DataType,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
      YY = GetNextValue(FileStream,Geom->DataType,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
      if (OkToLoad)
        {
          SetFieldValue(Geom->GeomData.XYZ.XBase,0,XX);
          SetFieldValue(Geom->GeomData.XYZ.YBase,0,YY);
        }
      Geom->NumSegments = 1;
      Geom->NumSegPts[0]   = 1;
    }
  else
    {
      double XX;
      CHECK((Geom->GeomType == GeomType_Square) ||
            (Geom->GeomType == GeomType_Circle));
      XX = GetNextValue(FileStream,Geom->DataType,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
      if (OkToLoad)
        {
          SetFieldValue(Geom->GeomData.XYZ.XBase,0,XX);
        }
      Geom->NumSegments  = 1;
      Geom->NumSegPts[0] = 1;
    }

  if (!IsOk)
    ErrMsg(_nts(ErrMsgString));

  return (IsOk);
} 

/*
 *  Pass1 for text simply acknowledges that a text was
 *  parsed and skips over the text.  
 */
Boolean_t ReadInText(FileStream_s *FileStream,
                     short         IVersion,
                     Boolean_t     OkToLoad,
                     Text_s       *Text,
                     LgIndex_t     MaxTextLen)
{
  LgIndex_t       I;
  FieldDataType_e FFT;
  SmInteger_t     TextLength = 0;
  Boolean_t       IsOk = TRUE;
  const char     *ErrMsgString = _ts2char(_ts("Invalid text record"));

  REQUIRE(VALID_REF(Text));

  if (IVersion < 70)
    FFT = FieldDataType_Float;
  else
    FFT = FieldDataType_Double;

  if (IVersion < 101)
    I = GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
  else
    I = GetIoFileInt(FileStream,IVersion,0,4,&IsOk);

  if ( I == 0 )
    Text->PositionCoordSys = CoordSys_Grid;
  else if ( I == 1 )
    Text->PositionCoordSys = CoordSys_Frame;
  /*
   * I == 2 is for CoordSys_FrameOffset and is not used currently
   *
   * I == 3 is for the old window coordinate system
   */
  else if ( I == 4 )
    Text->PositionCoordSys = CoordSys_Grid3D;
  else
    {
      ErrMsgString = _ts2char(_ts("Invalid text coordinate system."));
      IsOk = FALSE;
    }

  Text->Scope   = (Scope_e)GetIoFileInt(FileStream,IVersion,0,1,&IsOk);
  Text->AnchorPos.Generic.V1 = GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  Text->AnchorPos.Generic.V2 = GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  if (IVersion >= 101)
    Text->AnchorPos.Generic.V3 = GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  else
    Text->AnchorPos.Generic.V3 = 0.0; /* default value for pre 101 versions */

  if (IVersion > 40)
    Text->TextShape.Font = (Font_e)GetIoFileInt(FileStream,IVersion,0,(LgIndex_t)Font_CourierBold,&IsOk);
  else
    Text->TextShape.Font = Font_Helvetica;
  if (IVersion < 43)
    GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  if (IVersion < 70)
    {
      if (Text->PositionCoordSys == CoordSys_Grid)
        Text->TextShape.SizeUnits = Units_Grid;
      else
        Text->TextShape.SizeUnits = Units_Frame;
    }
  else
    Text->TextShape.SizeUnits = (Units_e)GetIoFileInt(FileStream,IVersion,0,(LgIndex_t)Units_Point,&IsOk);
    
  Text->TextShape.Height = GetNextValue(FileStream,FFT,-LARGEDOUBLE,LARGEDOUBLE,&IsOk);
  if (IVersion > 47)
    {
      Text->Box.BoxType = (TextBox_e)GetIoFileInt(FileStream,IVersion,0,(LgIndex_t)TextBox_Hollow,&IsOk);
      if (IVersion < 70)
        {
          if (Text->Box.BoxType == TextBox_Hollow)
            Text->Box.BoxType = TextBox_Filled;
          else if (Text->Box.BoxType == TextBox_Filled)
            Text->Box.BoxType = TextBox_Hollow;
        }
      Text->Box.Margin     = GetNextValue(FileStream,FFT,
                                       TextBoxMarginInputSpec.Min,
                                       TextBoxMarginInputSpec.Max,
                                       &IsOk);
      if (IVersion >= 70)
        Text->Box.LineThickness = GetNextValue(FileStream,FFT,
                                            LineThicknessInputSpec.Min,
                                            LineThicknessInputSpec.Max,
                                            &IsOk);
      else
        Text->Box.LineThickness = 0.01;
      Text->Box.BColor     = (ColorIndex_t)GetIoFileInt(FileStream,IVersion,0,255,&IsOk);
      Text->Box.FillBColor = (ColorIndex_t)GetIoFileInt(FileStream,IVersion,0,255,&IsOk);
      AdjustCustomColor(IVersion,&Text->Box.BColor);
      AdjustCustomColor(IVersion,&Text->Box.FillBColor);
    }
  else
    {
      Text->Box.BoxType    = TextBox_None;
      Text->Box.Margin     = 0.0;
      Text->Box.BColor     = White_C;
      Text->Box.FillBColor = Black_C;
    }
  if (IVersion < 70)
    {
      Text->Angle       = GetIoFileInt(FileStream,IVersion,-720,720,&IsOk)/DEGPERRADIANS;
      Text->LineSpacing = 1;
      Text->Anchor      = TextAnchor_Left;
    }
  else
    {
      Text->Angle       = GetNextValue(FileStream,FFT,
                                    TextAngleInputSpec.Min,
                                    TextAngleInputSpec.Max,
                                    &IsOk);
      Text->LineSpacing = GetNextValue(FileStream,FFT,
                                    TextLineSpacingInputSpec.Min,
                                    TextLineSpacingInputSpec.Max,
                                    &IsOk);
      Text->Anchor      = (TextAnchor_e)GetIoFileInt(FileStream,IVersion,0,(LgIndex_t)TextAnchor_HeadRight,&IsOk);
    }

  GetZoneAttachment(FileStream,IVersion,&Text->Zone,&Text->AttachToZone,&IsOk);

  Text->BColor   = (ColorIndex_t)GetIoFileInt(FileStream,IVersion,0,255,&IsOk);
  AdjustCustomColor(IVersion,&Text->BColor);
  if (IVersion < 70)
    TextLength = (short)GetIoFileInt(FileStream,IVersion,0,5000,&IsOk);

  ResetString(&Text->MacroFunctionCommand,NULL,TRUE);

  Text->Clipping = Clipping_ClipToViewport; /* default value for pre 101 versions */
 
  if (IVersion < 70)
    {
      short I,S;
      for (I = 0; I < TextLength; I++)
        {
          S = (short)GetIoFileInt(FileStream,IVersion,0,1000,&IsOk);
          if (OkToLoad && (I <= MaxTextLen))
            Text->Text[I] = (char)S;
        }
      if (OkToLoad)
        Text->Text[MIN(TextLength,MaxTextLen)] = '\0';
    }
  else
    {
      char *S = NULL;

      if (IVersion >= 75)
        {
          IsOk = ReadMacroFunctionCommand(FileStream,
                                          IVersion,
                                          OkToLoad,
                                          &Text->MacroFunctionCommand);
        } /* IVersion >= 75 */

      if (IVersion >= 101)
        {
          /*
           * The second clipping value was deprecated during v10 development and thus removed.
           * This moved Clipping_ClipToFrame to the 2nd position from the 3rd, so we convert
           * value 2 to ClipToFrame to support files made during v10 developement.
           */
          Text->Clipping = (Clipping_e)GetIoFileInt(FileStream,IVersion,0,2,&IsOk);
          if ( Text->Clipping == (Clipping_e)2 )
            Text->Clipping = Clipping_ClipToFrame;
        }

      if (ReadInString(FileStream,
                       IVersion,
                       MaxTextLen,
                       &S,
                       OkToLoad))
        {
          REQUIRE( !(S || OkToLoad) || VALID_REF(Text->Text));
          if (S)
            {
              #if defined TECPLOTKERNEL
                /* These symbols aren't available for the tecio library (yet, anyway). */
                if (IVersion < 101)
                  IsOk = ConvertOldTextFormattingToNew(&S,
                                                       Text->TextShape.Font,
                                                       MaxTextLen);
              #endif /* TECPLOTKERNEL */
              if (IsOk)
                {
                  strcpy(Text->Text,S);
                }
              FREE_ARRAY(S,"Release temp string for new text label");
            }
          else if (OkToLoad)
            Text->Text[0] = '\0';
        }
      else
        IsOk = FALSE;
    }

  if (!IsOk)
    ErrMsg(_nts(ErrMsgString));

  return (IsOk);
} 


static Boolean_t CompareVersion(float      Version,
                                char      *VersionString,
                                Boolean_t  IsByteOrderNative)
{
  char *VersionBuf = (char *)&Version;

  REQUIRE(VALID_REF(VersionString));

  if (IsByteOrderNative)
    return ((VersionString[0] == VersionBuf[0]) &&
            (VersionString[1] == VersionBuf[1]) &&
            (VersionString[2] == VersionBuf[2]) &&
            (VersionString[3] == VersionBuf[3]));
  else
    return ((VersionString[3] == VersionBuf[0]) &&
            (VersionString[2] == VersionBuf[1]) &&
            (VersionString[1] == VersionBuf[2]) &&
            (VersionString[0] == VersionBuf[3]));
}

static float ValidVersions[] = {7.0F,
                                6.3F, 6.2F, 6.1F, 6.0F,
                                5.0F,
                                4.7F, 4.6F, 4.5F, 4.4F, 4.3F, 4.2F, 4.1F, 4.0F};
#define NUMVALIDVERSIONS ((int)(sizeof(ValidVersions)/sizeof(ValidVersions[0])))


/* 
 * Extra caution taken here in case value read is invalid float 
 */
static Boolean_t GetDoubleVersion(char      *VersionString,
                                  float     *FInputVersion,
                                  Boolean_t  IsByteOrderNative)
{
  int  I;
  REQUIRE(VALID_REF(FInputVersion));
  
  for (I = 0; I < NUMVALIDVERSIONS; I++)
    if (CompareVersion(ValidVersions[I],VersionString,IsByteOrderNative))
      {
        *FInputVersion = ValidVersions[I];
        return (TRUE);
      }
  return (FALSE);
}


static short GetNewInputVersion(FileStream_s *FileStream)
{
  /*
   * 
   */
  char       Buf[4];
  short      IVersion = 0;
  short      I;
  LgIndex_t  OneValue;
  Boolean_t  IsOk = TRUE;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(FileStream->IsByteOrderNative);

  if (FREAD(Buf,4,1,FileStream->File) != 1)
    return (0);

  if (strncmp(Buf,"#!TD",4))
    return (0);

  if (FREAD(Buf,4,1,FileStream->File) != 1)
    return (0);

  if (Buf[0] != 'V')
    return (0);

  I = 1;
  while ((I < 4) && isdigit(Buf[I]))
    IVersion = IVersion*10 + Buf[I++] - '0';

  if (IVersion < 70)
    return (0);
  else if (IVersion > TecplotBinaryFileVersion)
    {
      ErrMsg(_ts("Binary file version newer than Tecplot version. "
             "Upgrade Tecplot or use an older Preplot to produce "
             "the datafile."));
      return (IVersion);
    }

  /*
   * Determine Byte Order.
   */

  OneValue = GetIoFileInt(FileStream,
                          IVersion,
                          -MAXINDEX,
                          MAXINDEX,
                          &IsOk);
 
  if (!IsOk)
    return (0);

  FileStream->IsByteOrderNative = (OneValue == 1);

  return (IVersion);
}

/**
 * Return value of zero is to be considered as an invalid 
 * tecplot binary datafile header.  Actually binary files
 * older than version 4.0 (return value of 40) are not supported
 * (See notes in preplot.c).
 */
short GetInputVersion(FileStream_s *FileStream)
{
  Boolean_t    IsOk = TRUE;
  float        FInputVersion;
  short        IVersion;
  char         VersionString[4];
  FileOffset_t StartOffset = 0;

  /*
   * First check to see if file uses new
   * input version format.
   */

  /* keep track of our start offset */
  StartOffset = FTELL(FileStream->File);

  IVersion = GetNewInputVersion(FileStream);

  if ( IVersion > TecplotBinaryFileVersion )
    return IVersion; /* unsupported version */
  else if (IVersion == 0)
    {
      /* rewind to clear any errors and seek to the start offset */
      rewind(FileStream->File);
      IsOk = (FSEEK(FileStream->File, StartOffset, SEEK_SET) == 0);

      if (IsOk && FREAD(VersionString,4,1,FileStream->File) == 1)
        {
          /* try both native and foreign versions numbers */
          if (!GetDoubleVersion(VersionString,&FInputVersion, FileStream->IsByteOrderNative))
            {
              FileStream->IsByteOrderNative = !FileStream->IsByteOrderNative; /* ...reverse the byte order */
              IsOk = GetDoubleVersion(VersionString,&FInputVersion, FileStream->IsByteOrderNative);
            }
          if (IsOk)
            IVersion = ROUNDS(FInputVersion*10);
        }
    }

  if (IsOk)
    return (IVersion);
  else
    return ((short)0);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE /* TOD(RMS)-H 12/12/2005: ENGINE: refactor to use just the Interrupted flag as-is */
#endif
#endif



/**********************************************************************
 **********************************************************************
 **********************        OUTPUT        **************************
 **********************************************************************
 **********************************************************************/

/**
 * Type Byte_t cannot be unaligned or reversed in byte order
 */
Boolean_t WriteBinaryByte(FileStream_s *FileStream,
                          Byte_t        ByteValue)
{
  Boolean_t IsOk;
  size_t    NumValuesWritten;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));

  NumValuesWritten = FWRITE(&ByteValue,sizeof(Byte_t),1,FileStream->File);
  IsOk = (Boolean_t)(NumValuesWritten == 1);

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 * Use Byte_t instead of Int16_t to support unaligned values
 */
static Boolean_t WriteBinaryInt16Unaligned(FileStream_s *FileStream,
                                           Byte_t       *ValueBuffer,
                                           Boolean_t     ValueInNativeOrder)
{
  Boolean_t IsOk;
  size_t    NumValuesWritten;
  Int16_t   Int16Value;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(ValueBuffer));
  REQUIRE(VALID_BOOLEAN(ValueInNativeOrder));

  if (ValueInNativeOrder != FileStream->IsByteOrderNative)
    COPY_AND_REVERSE_2_UNALIGNED_BYTES(&Int16Value, ValueBuffer);
  else
    COPY_2_UNALIGNED_BYTES(&Int16Value, ValueBuffer);

  NumValuesWritten = FWRITE(&Int16Value, sizeof(Int16_t), 1 ,FileStream->File);
  IsOk = (Boolean_t)(NumValuesWritten == 1);

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 * Use Byte_t instead of Int32_t to support unaligned values
 */
static Boolean_t WriteBinaryInt32Unaligned(FileStream_s *FileStream,
                                           Byte_t       *ValueBuffer,
                                           Boolean_t     ValueInNativeOrder)
{
  Boolean_t IsOk;
  size_t    NumValuesWritten;
  Int32_t   Int32Value;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(ValueBuffer));
  REQUIRE(VALID_BOOLEAN(ValueInNativeOrder));

  if (ValueInNativeOrder != FileStream->IsByteOrderNative)
    COPY_AND_REVERSE_4_UNALIGNED_BYTES(&Int32Value, ValueBuffer);
  else
    COPY_4_UNALIGNED_BYTES(&Int32Value, ValueBuffer);

  #if defined CRAY
    {
      int    TempI  = Int32Value;
      int    i      = 1;
      int    type   = 1; 
      int    bitoff = 0;
      char   Int4Buffer[4];

      CRAY2IEG(&type,&i,Int4Buffer,&bitoff,&TempI);
      NumValuesWritten = FWRITE(Int4Buffer, sizeof(Int32_t), 1, FileStream->File);
    }
  #else
    {
      NumValuesWritten = FWRITE(&Int32Value, sizeof(Int32_t), 1, FileStream->File);
    }
  #endif

  IsOk = (Boolean_t)(NumValuesWritten == 1);

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/*
 * This is used in many places and requires the value be in proper order.
 */
Boolean_t WriteBinaryInt32(FileStream_s *FileStream,
                           Int32_t       Value)
{
  Boolean_t IsOk;
  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE("Value can be any Int32_t");
  IsOk = WriteBinaryInt32Unaligned(FileStream, (Byte_t *)&Value, TRUE/*ValueInNativeOrder*/);
  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 * Use Byte_t instead of float to support unaligned values
 */
static Boolean_t WriteBinaryFloatUnaligned(FileStream_s *FileStream,
                                           Byte_t       *ValueBuffer,
                                           Boolean_t     ValueInNativeOrder)
{
  Boolean_t IsOk;
  size_t    NumValuesWritten;
  float     FloatValue;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(ValueBuffer));
  REQUIRE(VALID_BOOLEAN(ValueInNativeOrder));

  if (ValueInNativeOrder != FileStream->IsByteOrderNative)
    COPY_AND_REVERSE_4_UNALIGNED_BYTES(&FloatValue, ValueBuffer);
  else
    COPY_4_UNALIGNED_BYTES(&FloatValue, ValueBuffer);

  #if defined CRAY
    {
      double TempD  = (double)FloatValue;
      int    i      = 1;
      int    type   = 2; 
      int    bitoff = 0;
      char   Real4Buffer[4];

      CRAY2IEG(&type,&i,Real4Buffer,&bitoff,&TempD);
      NumValuesWritten = FWRITE(Real4Buffer,sizeof(float),1,FileStream->File);
    }
  #else
    {
      NumValuesWritten = FWRITE(&FloatValue,sizeof(float),1,FileStream->File);
    }
  #endif

  IsOk = (Boolean_t)(NumValuesWritten == 1);

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 */
static Boolean_t WriteBinaryDoubleUnaligned(FileStream_s *FileStream,
                                            Byte_t       *ValueBuffer,
                                            Boolean_t     ValueInNativeOrder)
{
  Boolean_t IsOk;
  size_t    NumValuesWritten;
  double    DoubleValue;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(ValueBuffer));
  REQUIRE(VALID_BOOLEAN(ValueInNativeOrder));

  if (ValueInNativeOrder != FileStream->IsByteOrderNative)
    COPY_AND_REVERSE_8_UNALIGNED_BYTES(&DoubleValue, ValueBuffer);
  else
    COPY_8_UNALIGNED_BYTES(&DoubleValue, ValueBuffer);

  #if defined CRAY
    {
      double TempD  = DoubleValue;
      int    i      = 1;
      int    type   = 8; 
      int    bitoff = 0;
      char   Real8Buffer[8];

      CRAY2IEG(&type,&i,Real8Buffer,&bitoff,&TempD);
      NumValuesWritten = FWRITE(Real8Buffer,sizeof(double),1,FileStream->File);
    }
  #else
    {
      NumValuesWritten = FWRITE(&DoubleValue,sizeof(double),1,FileStream->File);
    }
  #endif

  IsOk = (Boolean_t)(NumValuesWritten == 1);

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 * Byte blocks cannot be unaligned or reversed in bytes
 */
Boolean_t WriteBinaryByteBlock(FileStream_s *FileStream,
                               Byte_t       *ByteValues,
                               LgIndex_t     NumValues)
{
  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(ByteValues));
  REQUIRE(NumValues >= 0);

  size_t    NumValuesWritten = FWRITE(ByteValues,1,NumValues,FileStream->File);
  Boolean_t IsOk             = (NumValuesWritten == (size_t)NumValues);

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 * Use Byte_t instead of Int16_t to support unaligned values
 */
static Boolean_t WriteBinaryInt16BlockUnaligned(FileStream_s *FileStream,
                                                Byte_t       *Int16Values,
                                                LgIndex_t     NumValues,
                                                Boolean_t     ValuesInNativeOrdering)
{
  Boolean_t IsOk = TRUE;
  Boolean_t WriteEachValueSeparately;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(Int16Values));
  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_BOOLEAN(ValuesInNativeOrdering));

  #if defined CRAY
    WriteEachValueSeparately = TRUE;
  #else
    WriteEachValueSeparately = (ValuesInNativeOrdering != FileStream->IsByteOrderNative);
  #endif

  if (WriteEachValueSeparately)
    {
      LgIndex_t NIndex;
      for ( NIndex = 0; NIndex < NumValues; NIndex++ )
        if ( !WriteBinaryInt16Unaligned(FileStream, Int16Values+NIndex*sizeof(Int16_t), ValuesInNativeOrdering) )
          {
            IsOk = FALSE;
            break;
          }
    }
  else
    {
      size_t NumBytesToWrite = NumValues*sizeof(Int16_t);
      size_t NumBytesWritten = FWRITE(Int16Values,1,NumBytesToWrite,FileStream->File);
      IsOk = (NumBytesToWrite==NumBytesWritten);
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 * Use Byte_t instead of Int32_t to support unaligned values
 */
static Boolean_t WriteBinaryInt32BlockUnaligned(FileStream_s *FileStream,
                                                Byte_t       *Int32Values,
                                                LgIndex_t     NumValues,
                                                Boolean_t     ValuesInNativeOrdering)
{
  Boolean_t IsOk = TRUE;
  Boolean_t WriteEachValueSeparately;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(Int32Values));
  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_BOOLEAN(ValuesInNativeOrdering));

  #if defined CRAY
    WriteEachValueSeparately = TRUE;
  #else
    WriteEachValueSeparately = (ValuesInNativeOrdering != FileStream->IsByteOrderNative);
  #endif

  if (WriteEachValueSeparately)
    {
      LgIndex_t NIndex;
      for ( NIndex = 0; NIndex < NumValues; NIndex++ )
        if ( !WriteBinaryInt32Unaligned(FileStream, Int32Values+NIndex*sizeof(Int32_t), ValuesInNativeOrdering) )
          {
            IsOk = FALSE;
            break;
          }
    }
  else
    {
      size_t NumBytesToWrite = NumValues*sizeof(Int32_t);
      size_t NumBytesWritten = FWRITE(Int32Values,1,NumBytesToWrite,FileStream->File);
      IsOk = (NumBytesToWrite==NumBytesWritten);
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


/**
 * Use Byte_t instead of float to support unaligned values
 */
static Boolean_t WriteBinaryFloatBlockUnaligned(FileStream_s *FileStream,
                                                Byte_t       *FloatValues,
                                                LgIndex_t     NumValues,
                                                Boolean_t     ValuesInNativeOrdering)
{
  Boolean_t IsOk = TRUE;
  Boolean_t WriteEachValueSeparately;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(FloatValues));
  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_BOOLEAN(ValuesInNativeOrdering));

  #if defined CRAY
    WriteEachValueSeparately = TRUE;
  #else
    WriteEachValueSeparately = (ValuesInNativeOrdering != FileStream->IsByteOrderNative);
  #endif

  if (WriteEachValueSeparately)
    {
      LgIndex_t NIndex;
      for ( NIndex = 0; NIndex < NumValues; NIndex++ )
        if ( !WriteBinaryFloatUnaligned(FileStream, FloatValues+NIndex*sizeof(float), ValuesInNativeOrdering) )
          {
            IsOk = FALSE;
            break;
          }
    }
  else
    {
      size_t NumBytesToWrite = NumValues*sizeof(float);
      size_t NumBytesWritten = FWRITE(FloatValues,1,NumBytesToWrite,FileStream->File);
      IsOk = (NumBytesToWrite==NumBytesWritten);
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}

/**
 * Use Byte_t instead of double to support unaligned values
 */
static Boolean_t WriteBinaryDoubleBlockUnaligned(FileStream_s *FileStream,
                                                 Byte_t       *DoubleValues,
                                                 LgIndex_t     NumValues,
                                                 Boolean_t     ValuesInNativeOrdering)
{
  Boolean_t IsOk = TRUE;
  Boolean_t WriteEachValueSeparately;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(DoubleValues));
  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_BOOLEAN(ValuesInNativeOrdering));

  #if defined CRAY
    WriteEachValueSeparately = TRUE;
  #else
    WriteEachValueSeparately = (ValuesInNativeOrdering != FileStream->IsByteOrderNative);
  #endif

  if (WriteEachValueSeparately)
    {
      LgIndex_t NIndex;
      for ( NIndex = 0; NIndex < NumValues; NIndex++ )
        if ( !WriteBinaryDoubleUnaligned(FileStream, DoubleValues+NIndex*sizeof(double), ValuesInNativeOrdering) )
          {
            IsOk = FALSE;
            break;
          }
    }
  else
    {
      size_t NumBytesToWrite = NumValues*sizeof(double);
      size_t NumBytesWritten = FWRITE(DoubleValues,1,NumBytesToWrite,FileStream->File);
      IsOk = (NumBytesToWrite==NumBytesWritten);
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}

/**
 */
Boolean_t WriteBinaryReal(FileStream_s    *FileStream,
                          double           RR,
                          FieldDataType_e  FieldDataType)
{
  Boolean_t IsOk = FALSE; /* ...quite compiler */

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE((FieldDataType == FieldDataType_Float)  ||
          (FieldDataType == FieldDataType_Double) ||
          (FieldDataType == FieldDataType_Byte));

  switch (FieldDataType)
    {
      case FieldDataType_Float :
        {
          float FloatVal = CONVERT_DOUBLE_TO_FLOAT(RR);
          IsOk = WriteBinaryFloatUnaligned(FileStream, (Byte_t *)&FloatVal, TRUE/*NativeOrdering*/);
        } break;
      case FieldDataType_Double :
        {
          double DoubleVal = CLAMP_DOUBLE(RR);
          IsOk = WriteBinaryDoubleUnaligned(FileStream, (Byte_t *)&DoubleVal, TRUE/*NativeOrdering*/);
        } break;
      case FieldDataType_Byte :
        {
          /* Note: type Byte cannot be unaligned or reversed in bytes */
          Byte_t B;
          if (RR > 255)
            B = 255;
          else if (RR < 0)
            B = 0;
          else
            B = (Byte_t)RR;
          IsOk = WriteBinaryByte(FileStream, B);
        } break;
      default: CHECK(FALSE); break;
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


Boolean_t WriteFieldDataType(FileStream_s    *FileStream,
                             FieldDataType_e  FDT,
                             Boolean_t        WriteBinary)
{
  if (WriteBinary)
    return (WriteBinaryInt32(FileStream,(LgIndex_t)FDT));
  else
    {
      short S = 0;
      switch (FDT)
        {
          case FieldDataType_Float :  S = fprintf(FileStream->File,"SINGLE "); break;
          case FieldDataType_Double : S = fprintf(FileStream->File,"DOUBLE "); break;
          case FieldDataType_Int32 :  S = fprintf(FileStream->File,"LONGINT "); break;
          case FieldDataType_Int16 :  S = fprintf(FileStream->File,"SHORTINT "); break;
          case FieldDataType_Byte :   S = fprintf(FileStream->File,"BYTE "); break;
          case FieldDataType_Bit :    S = fprintf(FileStream->File,"BIT "); break;
          default: CHECK(FALSE);
        }
      return (FPRINTFOK(S));
    }
}


static Boolean_t WriteFieldDataValue(FileStream_s *FileStream,
                                     FieldData_pa  FieldData,
                                     LgIndex_t     Offset,
                                     Boolean_t     IsBinary,
                                     SmInteger_t   AsciiPrecision)
{
  Boolean_t IsOk = FALSE; /* ...quiet compiler */

  #if defined TECPLOTKERNEL
    REQUIRE(VALID_REF(FileStream) && VALID_FILE_HANDLE(FileStream->File));
    REQUIRE(VALID_REF(FieldData));
    REQUIRE(IsFieldDataUsable(FieldData));
    REQUIRE(Offset >= 0);
    REQUIRE(VALID_BOOLEAN(IsBinary));
    REQUIRE(AsciiPrecision >= 0);
    REQUIRE(VALID_BOOLEAN(IsOemLaunch));
    /* sanity check */
    INVARIANT(IMPLICATION(!FieldData->IsNativeOrdering, FieldData->MemMap!=NULL));
    INVARIANT(IMPLICATION(!FieldData->IsWordAligned, FieldData->MemMap!=NULL));
  #endif

  if (IsBinary)
    {
      size_t  NumBytes    = 0;
      Byte_t *ByteArray   = NULL;
      double  ValueBuffer = 0.0;

      /*
       * Use Byte_t instead of float, etc. to support unaligned values
       */
      switch ( GetFieldDataType(FieldData) )
        {
          case FieldDataType_Float :
            {
              if (FieldDataHasContiguousAddressSpace(FieldData))
                ByteArray = GetFieldDataBytePtr(FieldData) + Offset*sizeof(float);
              else
                {
                  /* data managed by an add-on has to be fetched via the callback */
                  ValueBuffer = GetFieldValue(FieldData, Offset);
                  ByteArray = (Byte_t *)&ValueBuffer;
                }
              IsOk = WriteBinaryFloatUnaligned(FileStream, ByteArray, FieldData->IsNativeOrdering);
              NumBytes = sizeof(float);
            } break;
          case FieldDataType_Double :
            {
              if (FieldDataHasContiguousAddressSpace(FieldData))
                ByteArray = GetFieldDataBytePtr(FieldData) + Offset*sizeof(double);
              else
                {
                  /* data managed by an add-on has to be fetched via the callback */
                  ValueBuffer = GetFieldValue(FieldData, Offset);
                  ByteArray = (Byte_t *)&ValueBuffer;
                }
              IsOk = WriteBinaryDoubleUnaligned(FileStream, ByteArray, FieldData->IsNativeOrdering);
              NumBytes = sizeof(double);
            } break;
          case FieldDataType_Int32 :
            {
              if (FieldDataHasContiguousAddressSpace(FieldData))
                ByteArray = GetFieldDataBytePtr(FieldData) + Offset*sizeof(Int32_t);
              else
                {
                  /* data managed by an add-on has to be fetched via the callback */
                  ValueBuffer = GetFieldValue(FieldData, Offset);
                  ByteArray = (Byte_t *)&ValueBuffer;
                }
              IsOk = WriteBinaryFloatUnaligned(FileStream, ByteArray, FieldData->IsNativeOrdering);
              NumBytes = sizeof(Int32_t);
            } break;
          case FieldDataType_Int16 :
            {
              if (FieldDataHasContiguousAddressSpace(FieldData))
                ByteArray = GetFieldDataBytePtr(FieldData) + Offset*sizeof(Int16_t);
              else
                {
                  /* data managed by an add-on has to be fetched via the callback */
                  ValueBuffer = GetFieldValue(FieldData, Offset);
                  ByteArray = (Byte_t *)&ValueBuffer;
                }
              IsOk = WriteBinaryFloatUnaligned(FileStream, ByteArray, FieldData->IsNativeOrdering);
              NumBytes = sizeof(Int16_t);
            } break;
          case FieldDataType_Byte :
            {
              /* type Byte cannot be unaligned or reversed in byte order */
              if (FieldDataHasContiguousAddressSpace(FieldData))
                ByteArray = GetFieldDataBytePtr(FieldData) + Offset;
              else
                {
                  /* data managed by an add-on has to be fetched via the callback */
                  ValueBuffer = GetFieldValue(FieldData, Offset);
                  ByteArray = (Byte_t *)&ValueBuffer;
                }
              IsOk = WriteBinaryByte(FileStream, *ByteArray);
              NumBytes = sizeof(Byte_t);
            } break;
          case FieldDataType_Bit :
            {
              /* type Bit cannot be unaligned or reversed in byte order */
              if ((Offset % 8) == 0)
                {
                  if (FieldDataHasContiguousAddressSpace(FieldData))
                    ByteArray = GetFieldDataBytePtr(FieldData) + Offset/8;
                  else
                    {
                      /* data managed by an add-on has to be fetched via the callback */
                      ValueBuffer = GetFieldValue(FieldData, Offset);
                      ByteArray = (Byte_t *)&ValueBuffer;
                    }
                  IsOk = WriteBinaryByte(FileStream, *ByteArray);
                  NumBytes = 1;
                }
              else
                {
                  IsOk = TRUE;
                }
            } break;
          default: CHECK(FALSE); break;
        }

      #if defined TECPLOTKERNEL
        {
          if ( IsOk && IsOemLaunch && ByteArray != NULL )
           {
              UPDATEOEMCHECKSUM(&GlobalOemChecksum,
                                ByteArray,
                                NumBytes,
                                1,
                                FALSE); /* Reverse bytes */       
           }
        }
      #endif
    }
  else
    {
      double V = 0.0;
      char buffer[100*MAX_SIZEOFUTF8CHAR];
      char *AsciiValue = buffer;

      #ifdef TECPLOTKERNEL
        size_t AsciiLen;
      #endif

      V = GetFieldValue(FieldData,Offset);

      switch ( GetFieldDataType(FieldData) )
        {
          case FieldDataType_Float :
          case FieldDataType_Double :
            /*IsOk = FPRINTFOK(fprintf(FileStream->File," %.*E",(int)AsciiPrecision,V)); */
            sprintf(buffer," %.*E",(int)AsciiPrecision,V);
            break;
          case FieldDataType_Int32 :
            /* IsOk = FPRINTFOK(fprintf(FileStream->File," %*d",(int)AsciiPrecision,ROUNDL(V))); */
            sprintf(buffer," %*d",(int)AsciiPrecision,ROUNDL(V));
            break;
          case FieldDataType_Int16 :
            /* IsOk = FPRINTFOK(fprintf(FileStream->File," %6d",ROUND2(V))); */
            sprintf(buffer," %6d",ROUND2(V));
            break;
          case FieldDataType_Byte :
            /* IsOk = FPRINTFOK(fprintf(FileStream->File," %3d",ROUNDS(V))); */
            sprintf(buffer," %3d",ROUNDS(V));
            break;
          case FieldDataType_Bit :
            /* IsOk = FPRINTFOK(fprintf(FileStream->File," %c",((V == 0) ? '0' : '1'))); */
            sprintf(buffer," %c",((V == 0) ? '0' : '1'));
            break;
          default: CHECK(FALSE); break;
        }
      IsOk = FPRINTFOK(fprintf(FileStream->File,buffer));
      #if defined TECPLOTKERNEL
        /* We have to strip leading spaces from the ascii numbers */
        AsciiLen = strlen(buffer);
        while ( *AsciiValue && isspace(*AsciiValue) )
          {
            AsciiValue++;
            AsciiLen--;
          }

        UPDATEOEMCHECKSUM(&GlobalOemChecksum,
                          AsciiValue,
                          AsciiLen,
                          1,
                          FALSE); /* Reverse Bytes */
      #endif

    }
  return (IsOk);
} 


/*
 */
Boolean_t WriteFieldDataBlock(FileStream_s *FileStream,
                              FieldData_pa  FieldData,
                              LgIndex_t     StartI,
                              LgIndex_t     NumValues)
{
  Boolean_t IsOk = FALSE; /* ...quite compiler */
  Byte_t *ByteArray = NULL;
  size_t  size = 0; /* for OEM checksums */
  size_t  nobj = NumValues; /* "" */

  #if defined TECPLOTKERNEL
    REQUIRE(VALID_REF(FileStream) && VALID_FILE_HANDLE(FileStream->File));
    REQUIRE(VALID_REF(FieldData));
    REQUIRE(IsFieldDataUsable(FieldData));
    REQUIRE(NumValues >= 0);
    REQUIRE(0 <= StartI && StartI + NumValues - 1 < GetFieldDataNumValues(FieldData));
    REQUIRE(VALID_BOOLEAN(IsOemLaunch));
    /* sanity check */
    INVARIANT(IMPLICATION(!FieldData->IsNativeOrdering, FieldData->MemMap!=NULL));
    INVARIANT(IMPLICATION(!FieldData->IsWordAligned, FieldData->MemMap!=NULL));
  #endif

  if (FieldDataHasContiguousAddressSpace(FieldData))
    {
      /*
       * Use Byte_t instead of float, etc. to support unaligned values
       */
      switch ( GetFieldDataType(FieldData) )
        {
          case FieldDataType_Float :
            {
              ByteArray = GetFieldDataBytePtr(FieldData) + StartI*sizeof(float);
              IsOk = WriteBinaryFloatBlockUnaligned(FileStream, ByteArray, NumValues, FieldData->IsNativeOrdering);
              size = sizeof(float);
            } break;
          case FieldDataType_Double :
            {
              ByteArray = GetFieldDataBytePtr(FieldData) + StartI*sizeof(double);
              IsOk = WriteBinaryDoubleBlockUnaligned(FileStream, ByteArray, NumValues, FieldData->IsNativeOrdering);
              size = sizeof(double);
            } break;
          case FieldDataType_Int32 :
            {
              ByteArray = GetFieldDataBytePtr(FieldData) + StartI*sizeof(Int32_t);
              IsOk = WriteBinaryInt32BlockUnaligned(FileStream, ByteArray, NumValues, FieldData->IsNativeOrdering);
              size = sizeof(Int32_t);
            } break;
          case FieldDataType_Int16 :
            {
              ByteArray = GetFieldDataBytePtr(FieldData) + StartI*sizeof(Int16_t);
              IsOk = WriteBinaryInt16BlockUnaligned(FileStream, ByteArray, NumValues, FieldData->IsNativeOrdering);
              size = sizeof(Int16_t);
            } break;
          case FieldDataType_Byte :
            {
              /* Note: Byte type cannot be unaligned or reversed in bytes */
              ByteArray = GetFieldDataBytePtr(FieldData) + StartI;
              IsOk = WriteBinaryByteBlock(FileStream, ByteArray, NumValues);
              size = sizeof(Byte_t);
            } break;
          case FieldDataType_Bit :
            {
              /* Note: Bit type cannot be unaligned or reversed in bytes */
              IsOk = (StartI == 0); /* ...only valid for StartI == 0 */
              if (IsOk)
                {
                  size_t NumBytes = 1+(NumValues-1)/8;
                  ByteArray = GetFieldDataBytePtr(FieldData);
                  IsOk = WriteBinaryByteBlock(FileStream, ByteArray, (LgIndex_t)NumBytes);
                  size = sizeof(Byte_t);
                  nobj = NumBytes;
                }
              else
                ErrMsg(_ts("Internal Error: Trying to dump block of bit values with start offset != 0"));
            } break;
          default: CHECK(FALSE); break;
        }

      #if defined TECPLOTKERNEL
        if (IsOemLaunch && ByteArray != NULL)
          {
            UPDATEOEMCHECKSUM(&GlobalOemChecksum,
                               ByteArray,
                               size,
                               nobj,
                               FALSE); /* Reverse Bytes */


          }
      #endif
    }
  else
    {
      /*
       * For data managed by an add-on we can't perform a block operation since
       * we don't have the raw data so we have to simulate it one value at a
       * time.
       */
      LgIndex_t II;
      for (II = 0, IsOk = TRUE; II < NumValues && IsOk; II++)
        IsOk = WriteFieldDataValue(FileStream,
                                   FieldData,
                                   StartI+II,
                                   TRUE, /* IsBinary */
                                   0);   /* AsciiPrecision */
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
} 


/**
 */
Boolean_t WriteCCFieldDataBlock(FileStream_s *FileStream,
                                FieldData_pa  FieldData,
                                Boolean_t     IsOrderedData,
                                LgIndex_t     NumIPts,
                                LgIndex_t     NumJPts,
                                LgIndex_t     NumKPts,
                                Boolean_t     WriteBinary,
                                SmInteger_t   AsciiPrecision)
{
  Boolean_t IsOk = TRUE;
  LgIndex_t NumValues;
  LgIndex_t I, J, K;
  LgIndex_t NumIJPts = -1;
  LgIndex_t IEnd     = -1;
  LgIndex_t JEnd     = -1;
  LgIndex_t KEnd     = -1;
  Boolean_t IsLinear = -1;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(FieldData));
  REQUIRE(VALID_BOOLEAN(IsOrderedData));
  REQUIRE(NumIPts >= 0);
  REQUIRE(NumJPts >= 0);
  REQUIRE(NumKPts >= 0);
  REQUIRE(VALID_BOOLEAN(WriteBinary));
  REQUIRE(IMPLICATION(!WriteBinary, AsciiPrecision >= 0));

  /*
   * As of version 103 Tecplot writes binary data files so that ordered cell
   * centered field data includes the ghost cells. This makes it much easier
   * for Tecplot to map the data when reading by simply writing out
   * FieldData->NumValues. As of version 104 the ghost cells of the slowest
   * moving index are not included but that does effect the output as it is
   * still FieldData->NumValues.
   */
  if (IsOrderedData && !WriteBinary)
    {
      /*
       * Ordered ASCII output is always layed out using
       * DataValueStructure_Classic format.
       */
      NumIJPts  = NumIPts * NumJPts;
      IEnd      = MAX(NumIPts-1, 1);
      JEnd      = MAX(NumJPts-1, 1);
      KEnd      = MAX(NumKPts-1, 1);
      NumValues = (IEnd * JEnd * KEnd);
      IsLinear  = ((NumJPts == 1 && NumKPts == 1) ||
                   (NumIPts == 1 && NumKPts == 1) ||
                   (NumIPts == 1 && NumJPts == 1));
    }
  else
    {
      NumValues = GetFieldDataNumValues(FieldData);
    }

  if (WriteBinary)
    {
      IsOk = WriteFieldDataBlock(FileStream, FieldData, 0, NumValues);
    }
  else
    {
      LgIndex_t NumValuesPerLine = 80/(AsciiPrecision+5);
      if (IsOrderedData && !IsLinear)
        {
          LgIndex_t ValueIndex = 0;
          for (K = 0; K < KEnd && IsOk; K++)
          for (J = 0; J < JEnd && IsOk; J++)
          for (I = 0; I < IEnd && IsOk; I++)
            {
              LgIndex_t CellIndex = I + (J*NumIPts) + (K*NumIJPts);
              IsOk = WriteFieldDataValue(FileStream, FieldData, CellIndex,
                                         WriteBinary, AsciiPrecision);
              if ((ValueIndex+1) % NumValuesPerLine == 0 || ValueIndex == NumValues-1)
                IsOk = (fputc('\n', FileStream->File) != EOF);
              ValueIndex++;
            }
        }
      else
        {
          for (I = 0; I < NumValues && IsOk; I++)
            {
              IsOk = WriteFieldDataValue(FileStream, FieldData, I,
                                         WriteBinary, AsciiPrecision);
              if ((I+1) % NumValuesPerLine == 0 || I == NumValues-1)
                IsOk = (fputc('\n', FileStream->File) != EOF);
            }
        }
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}


Boolean_t DumpDatafileString(FileStream_s *FileStream,
                             const char   *S,
                             Boolean_t     WriteBinary)
{
  Boolean_t IsOk = TRUE;
  const char *CPtr = S;
  if (WriteBinary)
    {
      const char *CPtr = S;
      while (IsOk && CPtr && *CPtr)
        IsOk = WriteBinaryInt32(FileStream,(LgIndex_t)(unsigned char)*CPtr++);
      if (IsOk)
        IsOk = WriteBinaryInt32(FileStream,0);
    }
  else
    {
      fputc('"',FileStream->File);
      while (CPtr && *CPtr)
        {
          if (*CPtr == '\n')
            {
              CPtr++;
              fputc('\\',FileStream->File);
              fputc('\\',FileStream->File);
              fputc('n',FileStream->File);
            }
          else
            {
              if ((*CPtr == '"') || (*CPtr == '\\'))
                fputc('\\',FileStream->File);
              fputc(*CPtr++,FileStream->File);
            }
        }
      fputc('"',FileStream->File);
      IsOk = (fputc('\n',FileStream->File) != EOF);
    }
  return (IsOk);
}


static void WriteAsciiColor(FILE        *File,
                            ColorIndex_t Color)
{
  if (Color >= FirstCustomColor && Color <= LastCustomColor)
    fprintf(File,"CUST%1d ",Color-FirstCustomColor+1);
  else
    {
      switch (Color)
        {
          case Black_C  : fprintf(File,"BLACK "); break;
          case Red_C    : fprintf(File,"RED "); break;
          case Green_C  : fprintf(File,"GREEN "); break;
          case Blue_C   : fprintf(File,"BLUE "); break;
          case Cyan_C   : fprintf(File,"CYAN "); break;
          case Yellow_C : fprintf(File,"YELLOW "); break;
          case Purple_C : fprintf(File,"PURPLE "); break;
          case White_C  : fprintf(File,"WHITE "); break;
        }
    }
}

static void WriteAsciiTextGeomBasics(FILE        *File,
                                     CoordSys_e   CoordSys,
                                     Boolean_t    AttachToZone,
                                     EntIndex_t   Zone,
                                     ColorIndex_t Color,
                                     Scope_e      Scope,
                                     Boolean_t    IncludeZ,
                                     Boolean_t    WriteGridDataAsPolar,
                                     AnchorPos_u *AnchorPos,
                                     double       ScaleFact)
{
  REQUIRE(VALID_REF(File));
  REQUIRE(VALID_TEXT_COORDSYS(CoordSys) || VALID_GEOM_COORDSYS(CoordSys));
  REQUIRE(VALID_BOOLEAN(AttachToZone));
  REQUIRE(IMPLICATION(AttachToZone, Zone>=0));
  REQUIRE(VALID_ENUM(Scope, Scope_e));
  REQUIRE(VALID_BOOLEAN(IncludeZ));
  REQUIRE(VALID_BOOLEAN(WriteGridDataAsPolar));
  REQUIRE(VALID_REF(AnchorPos));

  fprintf(File,"CS=");
  if (CoordSys == CoordSys_Frame)
    fprintf(File,"FRAME");
  else if (CoordSys == CoordSys_Grid)
    fprintf(File,"GRID");
  /*
   * Not currently used
   *
  else if (CoordSys == CoordSys_FrameOffset)
    fprintf(File,"FRAMEOFFSET");
   */
  else if (CoordSys == CoordSys_Grid3D)
    fprintf(File,"GRID3D");
  else
    CHECK(FALSE);

  if ( CoordSys==CoordSys_Grid && !IncludeZ && WriteGridDataAsPolar )
    {
      fprintf(File, "\nTHETA=%.12G,R=%.12G",
              ScaleFact*AnchorPos->ThetaR.Theta,
              ScaleFact*AnchorPos->ThetaR.R);
      CHECK(!IncludeZ);
    }
  else
    {
      fprintf(File,"\nX=%.12G,Y=%.12G",
              ScaleFact*AnchorPos->XYZ.X,
              ScaleFact*AnchorPos->XYZ.Y);
      if (IncludeZ)
        fprintf(File,",Z=%.12G", ScaleFact*AnchorPos->XYZ.Z);
    }

  if (AttachToZone)
    fprintf(File,"\nZN=%d",Zone+1);

  fprintf(File,"\nC=");
  WriteAsciiColor(File,Color);

  fprintf(File,"\nS=");
  if (Scope == Scope_Global)
    fprintf(File,"GLOBAL");
  else if (Scope == Scope_Local)
    fprintf(File,"LOCAL");
  else
    CHECK(FALSE);

  fputc('\n',File);
}


Boolean_t DumpGeometry(FileStream_s *FileStream,
                       Geom_s       *Geom,
                       Boolean_t     WriteBinary,
                       Boolean_t     WriteGridDataAsPolar)
{
  LgIndex_t       I,Index;
  LgIndex_t       SegIndex;
  Boolean_t       IsOk = TRUE;
  FieldDataType_e FDT;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(Geom));
  REQUIRE(Geom->GeomType != GeomType_Image);

  if (WriteBinary)
    {
      WriteBinaryReal(FileStream,GeomMarker,FieldDataType_Float);
      if (Geom->PositionCoordSys == CoordSys_Grid)
        WriteBinaryInt32(FileStream,0);
      else if (Geom->PositionCoordSys == CoordSys_Frame)
        WriteBinaryInt32(FileStream,1);
#if 0 /*
       * Not currently used
       */
      else if (Geom->PositionCoordSys == CoordSys_FrameOffset)
        WriteBinaryInt32(FileStream,2);
#endif
      /*
       * PositionCoordSys == 3 is for old window coordinate system
       */
      else if (Geom->PositionCoordSys == CoordSys_Grid3D)
        WriteBinaryInt32(FileStream,4);
      else
        CHECK(FALSE);

      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->Scope);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->DrawOrder);
      WriteBinaryReal(FileStream,Geom->AnchorPos.Generic.V1,FieldDataType_Double);
      WriteBinaryReal(FileStream,Geom->AnchorPos.Generic.V2,FieldDataType_Double);
      WriteBinaryReal(FileStream,Geom->AnchorPos.Generic.V3,FieldDataType_Double);
      if (Geom->AttachToZone)
        WriteBinaryInt32(FileStream,(LgIndex_t)Geom->Zone);
      else
        WriteBinaryInt32(FileStream,(LgIndex_t)-1);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->BColor);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->FillBColor);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->IsFilled);
      CHECK(Geom->GeomType != GeomType_LineSegs3D); /* deprecated */
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->GeomType);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->LinePattern);
      WriteBinaryReal(FileStream,Geom->PatternLength,FieldDataType_Double);
      WriteBinaryReal(FileStream,Geom->LineThickness,FieldDataType_Double);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->NumEllipsePts);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->ArrowheadStyle);
      WriteBinaryInt32(FileStream,(LgIndex_t)Geom->ArrowheadAttachment);
      WriteBinaryReal(FileStream,Geom->ArrowheadSize,FieldDataType_Double);

      WriteBinaryReal(FileStream,Geom->ArrowheadAngle,FieldDataType_Double);
      
      /* MACRO FUNCTION COMMAND */
      DumpDatafileString(FileStream,Geom->MacroFunctionCommand,TRUE);

      /*
       *  Assume geometry has X,Y (and Z) all using same field
       *  data type.
       */
      FDT = GetGeomFieldDataType(Geom);
      WriteFieldDataType(FileStream,FDT,TRUE);
      WriteBinaryInt32(FileStream, (LgIndex_t)Geom->Clipping);

      if (Geom->GeomType == GeomType_LineSegs)
        {
          short S;
          WriteBinaryInt32(FileStream,Geom->NumSegments);
          I = 0;
          for (S = 0; IsOk && (S < Geom->NumSegments); S++)
            {
              WriteBinaryInt32(FileStream,Geom->NumSegPts[S]);
              WriteFieldDataBlock(FileStream,Geom->GeomData.Generic.V1Base,I,Geom->NumSegPts[S]);
              IsOk = WriteFieldDataBlock(FileStream,Geom->GeomData.Generic.V2Base,I,Geom->NumSegPts[S]);
              if (Geom->PositionCoordSys == CoordSys_Grid3D)
                IsOk = WriteFieldDataBlock(FileStream,Geom->GeomData.Generic.V3Base,I,Geom->NumSegPts[S]);
              I += Geom->NumSegPts[S];
            }
        }
      else if (Geom->GeomType == GeomType_Rectangle ||
               Geom->GeomType == GeomType_Ellipse)
        {
          WriteBinaryReal(FileStream,GetFieldValue(Geom->GeomData.XYZ.XBase,0),FDT);
          IsOk = WriteBinaryReal(FileStream,GetFieldValue(Geom->GeomData.XYZ.YBase,0),FDT);
        }
      else
        {
          CHECK((Geom->GeomType == GeomType_Square) ||
                (Geom->GeomType == GeomType_Circle));
          IsOk = WriteBinaryReal(FileStream,GetFieldValue(Geom->GeomData.XYZ.XBase,0),FDT);
        }

    }
  else
    {
      double ScaleFact;
      if (Geom->PositionCoordSys == CoordSys_Frame)
        ScaleFact = 100.0;
      else
        ScaleFact = 1.0;

      fprintf(FileStream->File,"GEOMETRY\nF=POINT\n");
      WriteAsciiTextGeomBasics(FileStream->File,
                               Geom->PositionCoordSys,
                               Geom->AttachToZone,
                               Geom->Zone,
                               Geom->BColor,
                               Geom->Scope,
                               TRUE,
                               WriteGridDataAsPolar,
                               &Geom->AnchorPos,
                               ScaleFact);

      switch (Geom->LinePattern)
        {
          case LinePattern_Solid      : fprintf(FileStream->File,"L=SOLID\n"); break;
          case LinePattern_Dashed     : fprintf(FileStream->File,"L=DASHED\n"); break;
          case LinePattern_DashDot    : fprintf(FileStream->File,"L=DASHDOT\n"); break;
          case LinePattern_Dotted     : fprintf(FileStream->File,"L=DOTTED\n"); break;
          case LinePattern_LongDash   : fprintf(FileStream->File,"L=LONGDASH\n"); break;
          case LinePattern_DashDotDot : fprintf(FileStream->File,"L=DASHDOTDOT\n"); break;
          default: CHECK(FALSE); break;
        }
      fprintf(FileStream->File,"PL=%.12G\n",
              Geom->PatternLength*PatternLengthInputSpec.InterfaceAdjust.ScaleFact);
      fprintf(FileStream->File,"LT=%.12G\n",
              Geom->LineThickness*LineThicknessInputSpec.InterfaceAdjust.ScaleFact);

      if (Geom->IsFilled)
        {
          fprintf(FileStream->File,"FC=");
          WriteAsciiColor(FileStream->File,Geom->FillBColor);
        }

      if (Geom->Clipping == Clipping_ClipToViewport)
        fprintf(FileStream->File,"CLIPPING=CLIPTOVIEWPORT\n");
      else if (Geom->Clipping == Clipping_ClipToFrame)
        fprintf(FileStream->File,"CLIPPING=CLIPTOFRAME\n");
      else
        CHECK(FALSE);

      if (Geom->DrawOrder == DrawOrder_AfterData)
        fprintf(FileStream->File,"DRAWORDER=AFTERDATA\n");
      else if (Geom->DrawOrder == DrawOrder_BeforeData)
        fprintf(FileStream->File,"DRAWORDER=BEFOREDATA\n");
      else
        CHECK(FALSE);

      /* Macro function command */
      fprintf(FileStream->File,"MFC=");
      DumpDatafileString(FileStream,Geom->MacroFunctionCommand,FALSE);

      if ((Geom->GeomType == GeomType_Circle) || (Geom->GeomType == GeomType_Ellipse))
        fprintf(FileStream->File,"EP=%ld\n",(long)Geom->NumEllipsePts);

      if (Geom->GeomType == GeomType_LineSegs && Geom->PositionCoordSys!=CoordSys_Grid3D)
        {
          switch (Geom->ArrowheadStyle)
            {
              case Arrowhead_Plain  : fprintf(FileStream->File,"AST=PLAIN\n"); break;
              case Arrowhead_Filled : fprintf(FileStream->File,"AST=FILLED\n"); break;
              case Arrowhead_Hollow : fprintf(FileStream->File,"AST=HOLLOW\n"); break;
              default: CHECK(FALSE); break;
            }
        
          switch (Geom->ArrowheadAttachment)
            {
              case ArrowheadAttach_None        : break;
              case ArrowheadAttach_AtBeginning : fprintf(FileStream->File,"AAT=BEGINNING\n"); break;
              case ArrowheadAttach_AtEnd       : fprintf(FileStream->File,"AAT=END\n"); break;
              case ArrowheadAttach_AtBothEnds  : fprintf(FileStream->File,"AAT=BOTH\n"); break;
              default: CHECK(FALSE); break;
            }
          if (Geom->ArrowheadAttachment != ArrowheadAttach_None)
            {
              fprintf(FileStream->File,"ASZ=%.12G\n",
                      Geom->ArrowheadSize*ArrowheadSizeInputSpec.InterfaceAdjust.ScaleFact);
              fprintf(FileStream->File,"AAN=%.12G\n",
                      Geom->ArrowheadAngle*ArrowheadAngleInputSpec.InterfaceAdjust.ScaleFact);
            }
        }

      switch (Geom->GeomType)
        {
          case GeomType_LineSegs :
            {
              fprintf(FileStream->File,"T=LINE\n");
              fprintf(FileStream->File,"DT=");
              WriteFieldDataType(FileStream,GetFieldDataType(Geom->GeomData.Generic.V1Base),FALSE);
              fputc('\n',FileStream->File);
              fprintf(FileStream->File,"%d\n",(int)Geom->NumSegments);
              SegIndex = 0;
              for (I = 0; IsOk && (I < Geom->NumSegments); I++)
                {
                  fprintf(FileStream->File,"%ld\n",(long)Geom->NumSegPts[I]);
                  for (Index = 0; Index < Geom->NumSegPts[I]; Index++)
                    {
                      fprintf(FileStream->File,"%.12G ",GetFieldValue(Geom->GeomData.Generic.V1Base,SegIndex+Index)*ScaleFact);
                      fprintf(FileStream->File,"%.12G",GetFieldValue(Geom->GeomData.Generic.V2Base,SegIndex+Index)*ScaleFact);
                      if (Geom->PositionCoordSys == CoordSys_Grid3D)
                        IsOk = FPRINTFOK(fprintf(FileStream->File," %.12G\n",GetFieldValue(Geom->GeomData.Generic.V3Base,SegIndex+Index)));
                      else
                        IsOk = (Boolean_t)(fputc('\n',FileStream->File) != EOF);
                    }
                  SegIndex += Geom->NumSegPts[I];
                }
            } break;
          case GeomType_Rectangle :
            {
              fprintf(FileStream->File,"T=RECTANGLE %.12G %.12G\n",
                      GetFieldValue(Geom->GeomData.XYZ.XBase,0)*ScaleFact,
                      GetFieldValue(Geom->GeomData.XYZ.YBase,0)*ScaleFact);
            } break;
          case GeomType_Square :
            {
              fprintf(FileStream->File,"T=SQUARE %.12G\n",
                      GetFieldValue(Geom->GeomData.XYZ.XBase,0)*ScaleFact);
            } break;
          case GeomType_Circle :
            {
              fprintf(FileStream->File,"T=CIRCLE %.12G\n",
                      GetFieldValue(Geom->GeomData.XYZ.XBase,0)*ScaleFact);
            } break;
          case GeomType_Ellipse :
            {
              fprintf(FileStream->File,"T=ELLIPSE %.12G %.12G\n",
                      GetFieldValue(Geom->GeomData.XYZ.XBase,0)*ScaleFact,
                      GetFieldValue(Geom->GeomData.XYZ.YBase,0)*ScaleFact);
            } break;
          default: CHECK(FALSE);
        }
    }
  return (IsOk);
}

Boolean_t DumpText(FileStream_s *FileStream,
                   Text_s       *T,
                   Boolean_t     WriteBinary,
                   Boolean_t     WriteGridDataAsPolar)
{
  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_REF(T));
  REQUIRE(VALID_BOOLEAN(WriteBinary));
  REQUIRE(VALID_BOOLEAN(WriteGridDataAsPolar));

  if (WriteBinary)
    {
      WriteBinaryReal(FileStream,TextMarker,FieldDataType_Float);
      if (T->PositionCoordSys == CoordSys_Grid)
        WriteBinaryInt32(FileStream,0);
      else if (T->PositionCoordSys == CoordSys_Frame)
        WriteBinaryInt32(FileStream,1);
#if 0 /*
       * Not currently used
       */
      else if (Geom->PositionCoordSys == CoordSys_FrameOffset)
        WriteBinaryInt32(FileStream,2);
#endif
      /*
       * 3 is used for old window coordinate system
       */
      else if (T->PositionCoordSys == CoordSys_Grid3D)
        WriteBinaryInt32(FileStream,4);
      else
        CHECK(FALSE);

      WriteBinaryInt32(FileStream,(LgIndex_t)T->Scope);
      WriteBinaryReal(FileStream,T->AnchorPos.Generic.V1,FieldDataType_Double);
      WriteBinaryReal(FileStream,T->AnchorPos.Generic.V2,FieldDataType_Double);
      WriteBinaryReal(FileStream,T->AnchorPos.Generic.V3,FieldDataType_Double);
      WriteBinaryInt32(FileStream,(LgIndex_t)T->TextShape.Font);
      WriteBinaryInt32(FileStream,(LgIndex_t)T->TextShape.SizeUnits);
      WriteBinaryReal(FileStream,T->TextShape.Height,FieldDataType_Double);
      WriteBinaryInt32(FileStream,(LgIndex_t)T->Box.BoxType);
      WriteBinaryReal(FileStream,T->Box.Margin,FieldDataType_Double);
      WriteBinaryReal(FileStream,T->Box.LineThickness,FieldDataType_Double);
      WriteBinaryInt32(FileStream,(LgIndex_t)T->Box.BColor);
      WriteBinaryInt32(FileStream,(LgIndex_t)T->Box.FillBColor);
      WriteBinaryReal(FileStream,T->Angle,FieldDataType_Double);
      WriteBinaryReal(FileStream,T->LineSpacing,FieldDataType_Double);
      WriteBinaryInt32(FileStream,(LgIndex_t)T->Anchor);
      if (T->AttachToZone)
        WriteBinaryInt32(FileStream,(LgIndex_t)T->Zone);
      else
        WriteBinaryInt32(FileStream,(LgIndex_t)-1);
      WriteBinaryInt32(FileStream,(LgIndex_t)T->BColor);
    }
  else
    {
      double ScaleFact;
      Boolean_t IncludeZ = T->PositionCoordSys == CoordSys_Grid3D;
      if (T->PositionCoordSys == CoordSys_Frame)
        ScaleFact = 100.0;
      else
        ScaleFact = 1.0;
      fprintf(FileStream->File,"TEXT\n");
      WriteAsciiTextGeomBasics(FileStream->File,
                               T->PositionCoordSys,
                               T->AttachToZone,
                               T->Zone,
                               T->BColor,
                               T->Scope,
                               IncludeZ,
                               WriteGridDataAsPolar,
                               &T->AnchorPos,
                               ScaleFact);
      fprintf(FileStream->File,"HU=");
      switch (T->TextShape.SizeUnits)
        {
          case Units_Grid  : fprintf(FileStream->File,"GRID\n"); break;
          case Units_Frame : fprintf(FileStream->File,"FRAME\n"); break;
          case Units_Point : fprintf(FileStream->File,"POINT\n"); break;
          case Units_AxisPercentage : /* Not allowed */
          default: CHECK(FALSE); break;
        }

      fprintf(FileStream->File,"LS=%.4G ",T->LineSpacing);

      fprintf(FileStream->File,"AN=");
      switch (T->Anchor)
        {
          case TextAnchor_Left       : fprintf(FileStream->File,"LEFT\n");        break;
          case TextAnchor_Center     : fprintf(FileStream->File,"CENTER\n");      break;
          case TextAnchor_Right      : fprintf(FileStream->File,"RIGHT\n");       break;
          case TextAnchor_MidLeft    : fprintf(FileStream->File,"MIDLEFT\n");     break;
          case TextAnchor_MidCenter  : fprintf(FileStream->File,"MIDCENTER\n");   break;
          case TextAnchor_MidRight   : fprintf(FileStream->File,"MIDRIGHT\n");    break;
          case TextAnchor_HeadLeft   : fprintf(FileStream->File,"HEADLEFT\n");    break;
          case TextAnchor_HeadCenter : fprintf(FileStream->File,"HEADCENTER\n");  break;
          case TextAnchor_HeadRight  : fprintf(FileStream->File,"HEADRIGHT\n");   break;
          default: CHECK(FALSE); break;
        }

      switch (T->Box.BoxType)
        {
          case TextBox_Hollow : fprintf(FileStream->File,"BX=Hollow "); break;
          case TextBox_Filled : fprintf(FileStream->File,"BX=Filled "); break;
          default :;
        }
      fprintf(FileStream->File,"BXM=%.4G ",T->Box.Margin*100);
      fprintf(FileStream->File,"LT=%.4G ",T->Box.LineThickness*100.0);
      fprintf(FileStream->File,"BXO="); WriteAsciiColor(FileStream->File,T->Box.BColor);
      fprintf(FileStream->File,"BXF="); WriteAsciiColor(FileStream->File,T->Box.FillBColor);

      fprintf(FileStream->File,"\nF=");
      switch (T->TextShape.Font) 
        {
          case Font_Helvetica     :   fprintf(FileStream->File,"HELV");         break;
          case Font_HelveticaBold :   fprintf(FileStream->File,"HELV-BOLD");    break;
          case Font_Times    :        fprintf(FileStream->File,"TIMES");        break;
          case Font_TimesBold:        fprintf(FileStream->File,"TIMES-BOLD");   break;
          case Font_TimesItalic :     fprintf(FileStream->File,"TIMES-ITALIC"); break;
          case Font_TimesItalicBold : fprintf(FileStream->File,"TIMES-ITALIC-BOLD"); break;
          case Font_Courier  :        fprintf(FileStream->File,"COURIER");      break;
          case Font_CourierBold  :    fprintf(FileStream->File,"COURIER-BOLD"); break;
          case Font_Greek    :        fprintf(FileStream->File,"GREEK");        break;
          case Font_Math     :        fprintf(FileStream->File,"MATH");         break;
          case Font_UserDefined  :    fprintf(FileStream->File,"USER-DEF");     break;
          default: CHECK(FALSE); break;
        }
      if (T->TextShape.SizeUnits == Units_Frame)
        ScaleFact = 100.0;
      else
        ScaleFact = 1.0;
      fprintf(FileStream->File,"\nH=%.12G A=%.12G",
                     T->TextShape.Height*ScaleFact,
                     T->Angle*DEGPERRADIANS);
    }


  if (!WriteBinary)
    fprintf(FileStream->File,"\nMFC=");

  DumpDatafileString(FileStream,T->MacroFunctionCommand,WriteBinary);

  if (!WriteBinary)
    {
      if (T->Clipping == Clipping_ClipToViewport)
        fprintf(FileStream->File,"CLIPPING=CLIPTOVIEWPORT\n");
      else if (T->Clipping == Clipping_ClipToFrame)
        fprintf(FileStream->File,"CLIPPING=CLIPTOFRAME\n");
      else
        CHECK(FALSE);
    }
  else
    {
      WriteBinaryInt32(FileStream, (LgIndex_t)T->Clipping);
    }

  if (!WriteBinary)
    fprintf(FileStream->File,"T=");

  return (DumpDatafileString(FileStream,T->Text,WriteBinary));
}

Boolean_t DumpCustomAxisLabels(FileStream_s  *FileStream,
                               Boolean_t      WriteBinary,
                               StringList_pa  LabelBase)
{
  Boolean_t IsOk = TRUE;
  LgIndex_t Index = 0;
  LgIndex_t Count = 0;

  REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
  REQUIRE(VALID_BOOLEAN(WriteBinary));
  REQUIRE(StringListValid(LabelBase));

  Count = StringListCount(LabelBase);
  if (WriteBinary)
    {
      WriteBinaryReal(FileStream,CustomLabelMarker,FieldDataType_Float);
      WriteBinaryInt32(FileStream,Count);
    }
  else
    {
      fprintf(FileStream->File," CUSTOMLABELS = \n");
    }

  for (Index = 0, IsOk = TRUE; Index < Count && IsOk; Index++)
    {
      const char *CurLabel = StringListGetStringRef(LabelBase, Index);
      IsOk = DumpDatafileString(FileStream,CurLabel,WriteBinary);
    }

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}

Boolean_t WriteBinaryMagicAndVersion(FileStream_s *FileStream)
{
  char VerBuf[5];
  sprintf(VerBuf,"V%-3d",TecplotBinaryFileVersion);
  CHECK(strlen(VerBuf) == 4);
  fprintf(FileStream->File,"#!TD%s",VerBuf);
  /*
   * Write an integer value of 1 to the file.  This is used
   * by the reader to determine byte order of the file.
   */
  return (WriteBinaryInt32(FileStream,1));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if 0 /* first discovered commented out 10/31/2005 */
#endif
#if 0 /* first discovered commented out 10/31/2005 */
#endif
#if 0 /* not used yet */
#endif
#endif
