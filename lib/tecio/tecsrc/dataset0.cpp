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

#define DATASET0MODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "DATASET.h"
#include "SET.h"
#include "DATASHR.h"
#include "FILESTREAM.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "Q_MSG.h"
#include "DATASET0.h"

/*
 * Low level dataset functions.  No references to zones, vars or
 * the DataSet_s master structure here.
 */


/*
 */
void OutOfMemoryMsg(void)
{
  ErrMsg(_ts("Cannot allocate enough memory for this operation."));
} /* OutOfMemoryMsg() */


/**
 */
FieldData_pa FieldDataAlloc(void)
{
  FieldData_pa Result;

  Result = (FieldData_pa)ALLOC_ITEM(struct _FieldData_a, "FieldDataPtr");
  if ( Result != NULL )
    {
      Result->GetFunction            = NULL;
      Result->SetFunction            = NULL;
      Result->Data                   = NULL;
      Result->Type                   = FieldDataType_Invalid;
      Result->IsNativeOrdering       = FALSE;
      Result->IsWordAligned          = FALSE;
      Result->ValueLocation          = ValueLocation_Invalid;
      Result->RefCount               = 1; /* self */
      Result->VarShareRefCount       = 1; /* self */
      Result->NumValues              = 0;
      #if defined TECPLOTKERNEL /* TecIO doesn't require these features yet */
      Result->IsMinMaxValid          = FALSE;
      Result->MinValue               = 0.0;
      Result->MaxValue               = 0.0;

      Result->AccessGeneration       = -1;
      Result->MemMap                 = NULL;
      Result->IsFragmented           = FALSE;
      Result->IsRawDataSelfManaged   = FALSE;
      Result->IsReadableLock         = NULL;
      Result->IsWritableLock         = NULL;
      Result->IsLoaded               = FALSE;
      Result->IsReadOnly             = FALSE;
      Result->VariableLoad           = NULL;
      Result->VariableUnload         = NULL;
      Result->VariableCleanup        = NULL;
      Result->VariableClientData     = 0;

      /*
       * Make sure to mark this item as referenced so it doesn't get unloaded while
       * we are using it.
       */
      DataLoadMarkItemAsReferenced(Result);
      #endif
    }

  ENSURE(VALID_REF(Result) || Result == NULL);
  return Result;
}

/**
 * Most clients should not call this function but FieldDataCleanup() instead.
 * An exception to this would be Tecplot's own storable load-on-demand
 * functions.
 */
void FieldDataDeallocData(FieldData_pa FieldData)
{
  REQUIRE(VALID_REF(FieldData));

# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# else
  if (FieldData->Data != NULL)
    {
      FREE_ARRAY(FieldData->Data, "FieldData _Data");
      FieldData->Data = NULL;
    }
# endif
  ENSURE(FieldData->Data == NULL);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 */
void FieldDataCleanup(FieldData_pa FieldData)
{
  REQUIRE(VALID_REF(FieldData));
  REQUIRE("FieldData->IsReadableLock is not locked");
  REQUIRE("FieldData->IsWritableLock is not locked");

# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# else
  FieldDataDeallocData(FieldData);
# endif
}

/**
 */
void FieldDataDealloc(FieldData_pa *FieldData,
                      Boolean_t     DoTrackVarSharing)
{
  REQUIRE(VALID_REF(FieldData));
  REQUIRE(VALID_REF(*FieldData) || *FieldData == NULL);
  REQUIRE(IMPLICATION(*FieldData != NULL,
                      IsStructureReferenced(*FieldData)));
  REQUIRE(IMPLICATION(*FieldData != NULL && DoTrackVarSharing,
                      IsVarStructureReferenced(*FieldData)));
  REQUIRE(VALID_BOOLEAN(DoTrackVarSharing));
  REQUIRE(IMPLICATION(*FieldData != NULL,
                      (*FieldData)->RefCount >= (*FieldData)->VarShareRefCount));

  if (*FieldData != NULL)
    {
      if (DoTrackVarSharing)
        DecVarStructureReference(*FieldData);
      DecStructureReference(*FieldData);
      if (!IsStructureReferenced(*FieldData))
        {
          FieldDataCleanup(*FieldData);
          FREE_ITEM(*FieldData, "field data");
        }
      *FieldData = NULL;
    }

  ENSURE(*FieldData == NULL);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 * DstArray and SrcArray are aligned on proper word boundaries.
 */
void CopyTypedValueArray(FieldDataType_e  ValueType,
                         void            *DstArray,
                         LgIndex_t        DstStart,
                         LgIndex_t        DstSkip,
                         void            *SrcArray,
                         LgIndex_t        SrcStart,
                         LgIndex_t        SrcEnd,
                         LgIndex_t        SrcSkip)
{
  REQUIRE(VALID_FIELD_DATA_TYPE(ValueType) &&
          ValueType != FieldDataType_Bit);
  REQUIRE(VALID_REF(DstArray));
  REQUIRE(DstStart >= 0);
  REQUIRE(DstSkip > 0);
  REQUIRE(VALID_REF(SrcArray));
  REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd);
  REQUIRE(SrcSkip > 0);
  REQUIRE(DstArray != SrcArray);

  switch (ValueType)
    {
      case FieldDataType_Int64 : CHECK(FALSE); /* Future work: remove check */
      case FieldDataType_Double :
        {
          /* copy 8 bytes at a time */
          UInt64_t *SrcPtr = ((UInt64_t *)SrcArray) + SrcStart;
          UInt64_t *SrcPtrEnd = ((UInt64_t *)SrcArray) + SrcEnd;
          UInt64_t *DstPtr = ((UInt64_t *)DstArray) + DstStart;
          CHECK(sizeof(UInt64_t)==8 && sizeof(double)==8);
          while ( SrcPtr <= SrcPtrEnd )
            {
              *DstPtr = *SrcPtr;
              SrcPtr += SrcSkip;
              DstPtr += DstSkip;
            }
        } break;
      case FieldDataType_Float :
      case FieldDataType_Int32 :
       {
          /* copy 4 bytes at a time */
          UInt32_t *SrcPtr = ((UInt32_t *)SrcArray) + SrcStart;
          UInt32_t *SrcPtrEnd = ((UInt32_t *)SrcArray) + SrcEnd;
          UInt32_t *DstPtr = ((UInt32_t *)DstArray) + DstStart;
          CHECK(sizeof(UInt32_t)==4 && sizeof(float)==4);
          while ( SrcPtr <= SrcPtrEnd )
            {
              *DstPtr = *SrcPtr;
              SrcPtr += SrcSkip;
              DstPtr += DstSkip;
            }
        } break;
      case FieldDataType_Int16 :
        {
          /* copy 2 bytes at a time */
          UInt16_t *SrcPtr = ((UInt16_t *)SrcArray) + SrcStart;
          UInt16_t *SrcPtrEnd = ((UInt16_t *)SrcArray) + SrcEnd;
          UInt16_t *DstPtr = ((UInt16_t *)DstArray) + DstStart;
          CHECK(sizeof(UInt16_t)==2);
          while ( SrcPtr <= SrcPtrEnd )
            {
              *DstPtr = *SrcPtr;
              SrcPtr += SrcSkip;
              DstPtr += DstSkip;
            }
        } break;
      case FieldDataType_Byte :
        {
          /* copy 1 byte at a time */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart;
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd;
          Byte_t *DstPtr = ((Byte_t *)DstArray) + DstStart;
          while ( SrcPtr <= SrcPtrEnd )
            {
              *DstPtr = *SrcPtr;
              SrcPtr += SrcSkip;
              DstPtr += DstSkip;
            }
        } break;
      default : CHECK(FALSE);
    }
}

/**
 * Same as CopyTypedValueArray, but does extra work.  Doesn't assume
 * DstArray and SrcArray are aligned on proper word boundaries.
 */
void CopyUnalignedTypedValueArray(FieldDataType_e  ValueType,
                                  void            *DstArray,
                                  LgIndex_t        DstStart,
                                  LgIndex_t        DstSkip,
                                  void            *SrcArray,
                                  LgIndex_t        SrcStart,
                                  LgIndex_t        SrcEnd,
                                  LgIndex_t        SrcSkip)
{
  REQUIRE(VALID_FIELD_DATA_TYPE(ValueType) &&
          ValueType != FieldDataType_Bit);
  REQUIRE(VALID_REF(DstArray));
  REQUIRE(DstStart >= 0);
  REQUIRE(DstSkip > 0);
  REQUIRE(VALID_REF(SrcArray));
  REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd);
  REQUIRE(SrcSkip > 0);
  REQUIRE(DstArray != SrcArray);

  switch (ValueType)
    {
      case FieldDataType_Int64 : CHECK(FALSE); /* Future work: remove check */
      case FieldDataType_Double :
        {
          /* copy 8 bytes at a time */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart*sizeof(UInt64_t);
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd*sizeof(UInt64_t);
          Byte_t *DstPtr = ((Byte_t *)DstArray) + DstStart*sizeof(UInt64_t);
          size_t src_byte_skip = SrcSkip*sizeof(UInt64_t);
          size_t dst_byte_skip = DstSkip*sizeof(UInt64_t);
          CHECK(sizeof(UInt64_t)==8 && sizeof(double)==8);
          while ( SrcPtr <= SrcPtrEnd )
            {
              COPY_8_UNALIGNED_BYTES(DstPtr, SrcPtr);
              SrcPtr += src_byte_skip;
              DstPtr += dst_byte_skip;
            }
        } break;
      case FieldDataType_Float :
      case FieldDataType_Int32 :
       {
          /* copy 4 bytes at a time */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart*sizeof(UInt32_t);
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd*sizeof(UInt32_t);
          Byte_t *DstPtr = ((Byte_t *)DstArray) + DstStart*sizeof(UInt32_t);
          size_t src_byte_skip = SrcSkip*sizeof(UInt32_t);
          size_t dst_byte_skip = DstSkip*sizeof(UInt32_t);
          CHECK(sizeof(UInt32_t)==4 && sizeof(float)==4);
          while ( SrcPtr <= SrcPtrEnd )
            {
              COPY_4_UNALIGNED_BYTES(DstPtr, SrcPtr);
              SrcPtr += src_byte_skip;
              DstPtr += dst_byte_skip;
            }
        } break;
      case FieldDataType_Int16 :
        {
          /* copy 2 bytes at a time */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart*sizeof(UInt16_t);
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd*sizeof(UInt16_t);
          Byte_t *DstPtr = ((Byte_t *)DstArray) + DstStart*sizeof(UInt16_t);
          size_t src_byte_skip = SrcSkip*sizeof(UInt16_t);
          size_t dst_byte_skip = DstSkip*sizeof(UInt16_t);
          CHECK(sizeof(UInt16_t)==2);
          while ( SrcPtr <= SrcPtrEnd )
            {
              COPY_2_UNALIGNED_BYTES(DstPtr, SrcPtr);
              SrcPtr += src_byte_skip;
              DstPtr += dst_byte_skip;
            }
        } break;
      case FieldDataType_Byte :
        {
          /* copy 1 byte at a time */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart;
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd;
          Byte_t *DstPtr = ((Byte_t *)DstArray) + DstStart;
          while ( SrcPtr <= SrcPtrEnd )
            {
              *DstPtr = *SrcPtr;
              SrcPtr += SrcSkip;
              DstPtr += DstSkip;
            }
        } break;
      default : CHECK(FALSE);
    }
}

/**
 * SrcArray is aligned on proper word boundaries.
 */
void SwapBytesInTypedValueArray(FieldDataType_e  ValueType,
                                void            *SrcArray,
                                LgIndex_t        SrcStart,
                                LgIndex_t        SrcEnd,
                                LgIndex_t        SrcSkip)
{
  REQUIRE(VALID_FIELD_DATA_TYPE(ValueType) &&
          ValueType != FieldDataType_Bit);
  REQUIRE(VALID_REF(SrcArray));
  REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd);
  REQUIRE(SrcSkip > 0);

  switch (ValueType)
    {
      case FieldDataType_Int64: CHECK(FALSE); /* Future work: remove CHECK */
      case FieldDataType_Double:
        {
          /* swap 8 bytes blocks */
          UInt64_t *SrcPtr = ((UInt64_t *)SrcArray) + SrcStart;
          UInt64_t *SrcPtrEnd = ((UInt64_t *)SrcArray) + SrcEnd;
          CHECK(sizeof(UInt64_t)==8 && sizeof(double)==8);
          while ( SrcPtr <= SrcPtrEnd )
            {
              REVERSE_8_BYTES(SrcPtr);
              SrcPtr += SrcSkip;
            }
        } break;
      case FieldDataType_Float:
      case FieldDataType_Int32:
        {
          /* swap 4 bytes blocks */
          UInt32_t *SrcPtr = ((UInt32_t *)SrcArray) + SrcStart;
          UInt32_t *SrcPtrEnd = ((UInt32_t *)SrcArray) + SrcEnd;
          CHECK(sizeof(UInt32_t)==4 && sizeof(float)==4);
          while ( SrcPtr <= SrcPtrEnd )
            {
              REVERSE_4_BYTES(SrcPtr);
              SrcPtr += SrcSkip;
            }
        } break;
      case FieldDataType_Int16:
        {
          /* swap 4 bytes blocks */
          UInt16_t *SrcPtr = ((UInt16_t *)SrcArray) + SrcStart;
          UInt16_t *SrcPtrEnd = ((UInt16_t *)SrcArray) + SrcEnd;
          CHECK(sizeof(UInt16_t)==2);
          while ( SrcPtr <= SrcPtrEnd )
            {
              REVERSE_2_BYTES(SrcPtr);
              SrcPtr += SrcSkip;
            }
        } break;
      case FieldDataType_Byte:
      case FieldDataType_Bit:
        {
          /* nothing to do */
        } break;
      default: CHECK(FALSE);
    }
}

/**
 * Same as SwapBytesInTypedValueArray, but does extra work.  Doesn't assume
 * DstArray and SrcArray are aligned on proper word boundaries.
 */
void SwapBytesInUnalignedTypedValueArray(FieldDataType_e  ValueType,
                                         void            *SrcArray,
                                         LgIndex_t        SrcStart,
                                         LgIndex_t        SrcEnd,
                                         LgIndex_t        SrcSkip)
{
  REQUIRE(VALID_FIELD_DATA_TYPE(ValueType) &&
          ValueType != FieldDataType_Bit);
  REQUIRE(VALID_REF(SrcArray));
  REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd);
  REQUIRE(SrcSkip > 0);

  switch (ValueType)
    {
      case FieldDataType_Int64: CHECK(FALSE); /* Future work: remove CHECK */
      case FieldDataType_Double:
        {
          /* swap 8-byte blocks */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart*sizeof(UInt64_t);
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd*sizeof(UInt64_t);
          size_t byte_skip = SrcSkip*sizeof(UInt64_t);
          CHECK(sizeof(UInt64_t)==8 && sizeof(double)==8);
          while ( SrcPtr <= SrcPtrEnd )
            {
              REVERSE_8_BYTES_1_AT_A_TIME(SrcPtr);
              SrcPtr += byte_skip;
            }
        } break;
      case FieldDataType_Float:
      case FieldDataType_Int32:
        {
          /* swap 4-byte blocks */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart*sizeof(UInt32_t);
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd*sizeof(UInt32_t);
          size_t byte_skip = SrcSkip*sizeof(UInt32_t);
          CHECK(sizeof(UInt32_t)==4 && sizeof(float)==4);
          while ( SrcPtr <= SrcPtrEnd )
            {
              REVERSE_4_BYTES_1_AT_A_TIME(SrcPtr);
              SrcPtr += byte_skip;
            }
        } break;
      case FieldDataType_Int16:
        {
          /* swap 2-byte blocks */
          Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart*sizeof(UInt16_t);
          Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd*sizeof(UInt16_t);
          size_t byte_skip = SrcSkip*sizeof(UInt16_t);
          CHECK(sizeof(UInt16_t)==2);
          while ( SrcPtr <= SrcPtrEnd )
            {
              REVERSE_2_BYTES_1_AT_A_TIME(SrcPtr);
              SrcPtr += byte_skip;
            }
        } break;
      case FieldDataType_Byte:
      case FieldDataType_Bit:
        {
          /* No swapping required. */
        } break;
      default: CHECK(FALSE);
    }
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if 0
#else
#endif
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined DEBUG_FIELDVALUES
# define DEBUG_FIELDVALUES_BAD_VALUE 0x11
  static unsigned char BadValueStr[] =
    {
      DEBUG_FIELDVALUES_BAD_VALUE,
      DEBUG_FIELDVALUES_BAD_VALUE,
      DEBUG_FIELDVALUES_BAD_VALUE,
      DEBUG_FIELDVALUES_BAD_VALUE,
      DEBUG_FIELDVALUES_BAD_VALUE,
      DEBUG_FIELDVALUES_BAD_VALUE,
      DEBUG_FIELDVALUES_BAD_VALUE,
      DEBUG_FIELDVALUES_BAD_VALUE
    };
  /**
   * If Tecplot is responsible for managing (allocating and deallocating) the
   * raw data then
   */
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# else
#   define FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, type) \
                   ((sizeof(type) < 4) /* cannot make reliably test with less than four bytes */ || \
                    memcmp(BadValueStr,((char *)((fd)->Data))+sizeof(type)*(pt), sizeof(type)) != 0)
# endif
#else
# define FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, type) TRUE
#endif



/**
 * Used in macros, thus not static
 */
double STDCALL GetFieldValueForFloat(const FieldData_pa fd,
                                     LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));

  Result = (double)GetFieldDataFloatPtr(fd)[pt];

  return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
static double STDCALL GetFieldValueForFloatReversed(const FieldData_pa fd,
                                                    LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));

  /* reverse byte order */
    {
      FloatRev_t float_rev_val;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        float_rev_val = GetFieldDataFloatRevPtr(fd)[pt];
      REVERSE_4_BYTES(&float_rev_val);
      Result = (double)(*(float *)&float_rev_val);
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForFloatUnaligned(const FieldData_pa fd,
                                                     LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));

  /* get bytes into alignment */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(float);

      float float_val;
      COPY_4_UNALIGNED_BYTES(&float_val, byte_array);
      Result = (double)float_val;
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForFloatReversedUnaligned(const FieldData_pa fd,
                                                             LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));

  /* get bytes into alignment and reverse byte order */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(float);

      float float_val;
      COPY_AND_REVERSE_4_UNALIGNED_BYTES(&float_val, byte_array);
      Result = (double)float_val;
    }

  return Result;
}


/**
 * Used in macros, thus not static
 */
double STDCALL GetFieldValueForDouble(const FieldData_pa fd,
                                      LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));

  Result = GetFieldDataDoublePtr(fd)[pt];

  return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
static double STDCALL GetFieldValueForDoubleReversed(const FieldData_pa fd,
                                                     LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));

  /* reverse byte order */
    {
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        *(DoubleRev_t *)&Result = GetFieldDataDoubleRevPtr(fd)[pt];
      REVERSE_8_BYTES(&Result);
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForDoubleUnaligned(const FieldData_pa fd,
                                                      LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));

  /* get bytes into alignment */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(double);

      COPY_8_UNALIGNED_BYTES(&Result, byte_array);
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForDoubleReversedUnaligned(const FieldData_pa fd,
                                                              LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));

  /* get bytes into alignment and reverse byte order */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(double);

      COPY_8_UNALIGNED_BYTES(&Result, byte_array);
      REVERSE_8_BYTES(&Result);
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt32(const FieldData_pa fd,
                                            LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    Result = (double)GetFieldDataInt32Ptr(fd)[pt];

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt32Reversed(const FieldData_pa fd,
                                                    LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));

  /* reverse byte order */
    {
      Int32Rev_t int32_val_rev;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        int32_val_rev = GetFieldDataInt32RevPtr(fd)[pt];

      REVERSE_4_BYTES(&int32_val_rev);
      Result = (double)(*(Int32_t *)&int32_val_rev);
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt32Unaligned(const FieldData_pa fd,
                                                     LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));

  /* get bytes into alignment */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int32_t);

      Int32_t int32_val;
      COPY_4_UNALIGNED_BYTES(&int32_val, byte_array);
      Result = (double)int32_val;
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt32ReversedUnaligned(const FieldData_pa fd,
                                                             LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));

  /* get bytes into alignment and reverse byte order */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int32_t);

      Int32_t int32_val;
      COPY_AND_REVERSE_4_UNALIGNED_BYTES(&int32_val, byte_array);
      Result = (double)int32_val;
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt16(const FieldData_pa fd,
                                            LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    Result = (double)GetFieldDataInt16Ptr(fd)[pt];

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt16Reversed(const FieldData_pa fd,
                                                    LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));

  /* reverse byte order */
    {
      Int16Rev_t int16_val_rev;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        int16_val_rev = GetFieldDataInt16RevPtr(fd)[pt];

      REVERSE_2_BYTES(&int16_val_rev);
      Result = (double)(*(Int16Rev_t *)&int16_val_rev);
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt16Unaligned(const FieldData_pa fd,
                                                     LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));

  /* get bytes into alignment */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int16_t);

      Int16_t int16_val;
      COPY_2_UNALIGNED_BYTES(&int16_val, byte_array);
      Result = (double)int16_val;
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForInt16ReversedUnaligned(const FieldData_pa fd,
                                                             LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));

  /* get bytes into alignment and reverse byte order */
    {
      Byte_t *byte_array;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int16_t);

      Int16_t int16_val;
      COPY_AND_REVERSE_2_UNALIGNED_BYTES(&int16_val, byte_array);
      Result = (double)int16_val;
    }

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForByte(const FieldData_pa fd,
                                           LgIndex_t          pt)
{
  double Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Byte);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Byte_t));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    Result = (double)GetFieldDataBytePtr(fd)[pt];

  return Result;
}


/**
 */
static double STDCALL GetFieldValueForBit(const FieldData_pa fd,
                                          LgIndex_t          pt)
{
  double    Result;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Bit);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt/8, Byte_t));

  LgIndex_t ByteOffset = pt/8;
  Byte_t    BitMask    = (01 << (pt%8));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd);
      if ( byte_array[ByteOffset] & BitMask )
        Result = 1.0;
      else
        Result = 0.0;
    }

  return Result;
}


/**
 */
FieldValueGetFunction_pf DetermineFieldDataGetFunction(FieldDataType_e DataType,
                                                       Boolean_t       IsNativeOrdering,
                                                       Boolean_t       IsWordAligned)
{
  FieldValueGetFunction_pf Result;

  REQUIRE(VALID_FIELD_DATA_TYPE(DataType));
  REQUIRE(VALID_BOOLEAN(IsNativeOrdering));
  REQUIRE(VALID_BOOLEAN(IsWordAligned));

  switch ( DataType )
    {
      case FieldDataType_Float :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = GetFieldValueForFloat;
              else
                Result = GetFieldValueForFloatUnaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = GetFieldValueForFloatReversed;
              else
                Result = GetFieldValueForFloatReversedUnaligned;
            }
        } break;
      case FieldDataType_Double :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = GetFieldValueForDouble;
              else
                Result = GetFieldValueForDoubleUnaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = GetFieldValueForDoubleReversed;
              else
                Result = GetFieldValueForDoubleReversedUnaligned;
            }
        } break;
      case FieldDataType_Int32 :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = GetFieldValueForInt32;
              else
                Result = GetFieldValueForInt32Unaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = GetFieldValueForInt32Reversed;
              else
                Result = GetFieldValueForInt32ReversedUnaligned;
            }
        } break;
      case FieldDataType_Int16 :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = GetFieldValueForInt16;
              else
                Result = GetFieldValueForInt16Unaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = GetFieldValueForInt16Reversed;
              else
                Result = GetFieldValueForInt16ReversedUnaligned;
            }
        } break;
      case FieldDataType_Byte :
        {
          /* Bytes cannot be misaligned or reversed */
          Result = GetFieldValueForByte;
        } break;
      case FieldDataType_Bit :
        {
          /* Bits cannot be misaligned or reversed */
          Result = GetFieldValueForBit;
        } break;
      default :
        {
          CHECK(FALSE);
          Result = NULL; /* satisfy compiler */
        } break;
    }
  return (Result);
}

/**
 */
static void STDCALL SetFieldValueForFloat(FieldData_pa fd,
                                          LgIndex_t    pt,
                                          double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

  GetFieldDataFloatPtr(fd)[pt] = CONVERT_DOUBLE_TO_FLOAT(val);

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/*
 */
static void STDCALL SetFieldValueForFloatReversed(FieldData_pa fd,
                                                  LgIndex_t    pt,
                                                  double       val)
{
  FloatRev_t FloatRevVal;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

  /* We set FloatRevVal because when we switch bytes, the value may no longer be a valid float. */
  *((float *)&FloatRevVal) = CONVERT_DOUBLE_TO_FLOAT(val);

  REVERSE_4_BYTES(&FloatRevVal);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    GetFieldDataFloatRevPtr(fd)[pt] = FloatRevVal;

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, FloatRev_t));
}


/*
 */
static void STDCALL SetFieldValueForFloatUnaligned(FieldData_pa fd,
                                                   LgIndex_t    pt,
                                                   double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(float);
      float float_val = CONVERT_DOUBLE_TO_FLOAT(val);
      COPY_4_UNALIGNED_BYTES(byte_array, &float_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));
}


/*
 */
static void STDCALL SetFieldValueForFloatReversedUnaligned(FieldData_pa fd,
                                                           LgIndex_t    pt,
                                                           double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Float && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(float);
      float float_val = CONVERT_DOUBLE_TO_FLOAT(val);
      REVERSE_4_BYTES(&float_val);
      COPY_4_UNALIGNED_BYTES(byte_array, &float_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));
}


/**
 */
static void STDCALL SetFieldValueForDouble(FieldData_pa fd,
                                           LgIndex_t    pt,
                                           double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

  GetFieldDataDoublePtr(fd)[pt] = CLAMP_DOUBLE(val);

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/*
 */
static void STDCALL SetFieldValueForDoubleReversed(FieldData_pa fd,
                                                   LgIndex_t    pt,
                                                   double       val)
{
  DoubleRev_t DoubleRevVal;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

  /* We set DoubleRevVal because when we switch bytes, the value may no longer be a valid double. */
  *((double *)&DoubleRevVal) = CLAMP_DOUBLE(val);

  REVERSE_8_BYTES(&DoubleRevVal);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    GetFieldDataDoubleRevPtr(fd)[pt] = DoubleRevVal;

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));
}


/*
 */
static void STDCALL SetFieldValueForDoubleUnaligned(FieldData_pa fd,
                                                    LgIndex_t    pt,
                                                    double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(double);
      double double_val = CLAMP_DOUBLE(val);
      COPY_8_UNALIGNED_BYTES(byte_array, &double_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));
}


/*
 */
static void STDCALL SetFieldValueForDoubleReversedUnaligned(FieldData_pa fd,
                                                            LgIndex_t    pt,
                                                            double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Double && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(double);
      double double_val = CLAMP_DOUBLE(val);
      COPY_AND_REVERSE_8_UNALIGNED_BYTES(byte_array, &double_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));
}


/*
 */
static void STDCALL SetFieldValueForInt32(FieldData_pa fd,
                                          LgIndex_t    pt,
                                          double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    GetFieldDataInt32Ptr(fd)[pt] = CONVERT_DOUBLE_TO_INT32(val);

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));
}


/*
 */
static void STDCALL SetFieldValueForInt32Reversed(FieldData_pa fd,
                                                  LgIndex_t    pt,
                                                  double       val)
{
  Int32Rev_t Int32RevVal;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

  /* We set Int32RevVal here to be consistant with float and doubles. */
  *((Int32_t *)&Int32RevVal) = CONVERT_DOUBLE_TO_INT32(val);

  REVERSE_4_BYTES(&Int32RevVal);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    GetFieldDataInt32RevPtr(fd)[pt] = Int32RevVal;

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));
}


/*
 */
static void STDCALL SetFieldValueForInt32Unaligned(FieldData_pa fd,
                                                   LgIndex_t    pt,
                                                   double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int32_t);
      Int32_t int32_val = CONVERT_DOUBLE_TO_INT32(val);
      COPY_4_UNALIGNED_BYTES(byte_array, &int32_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));
}


/*
 */
static void STDCALL SetFieldValueForInt32ReversedUnaligned(FieldData_pa fd,
                                                           LgIndex_t    pt,
                                                           double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int32 && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int32_t);
      Int32_t int32_val = CONVERT_DOUBLE_TO_INT32(val);
      REVERSE_4_BYTES(&int32_val);
      COPY_4_UNALIGNED_BYTES(byte_array, &int32_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));
}


/*
 */
static void STDCALL SetFieldValueForInt16(FieldData_pa fd,
                                          LgIndex_t    pt,
                                          double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    GetFieldDataInt16Ptr(fd)[pt] = CONVERT_DOUBLE_TO_INT16(val);

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));
}


/*
 */
static void STDCALL SetFieldValueForInt16Reversed(FieldData_pa fd,
                                                  LgIndex_t    pt,
                                                  double       val)
{
  Int16Rev_t Int16RevVal;

  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && !fd->IsNativeOrdering && fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

  /* We set Int16RevVal here to be consistant with float and doubles. */
  *((Int16_t *)&Int16RevVal) = CONVERT_DOUBLE_TO_INT16(val);
  REVERSE_2_BYTES(&Int16RevVal);
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    GetFieldDataInt16RevPtr(fd)[pt] = Int16RevVal;

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));
}


/*
 */
static void STDCALL SetFieldValueForInt16Unaligned(FieldData_pa fd,
                                                   LgIndex_t    pt,
                                                   double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int16_t);
      Int16_t int16_val = CONVERT_DOUBLE_TO_INT16(val);
      COPY_2_UNALIGNED_BYTES(byte_array, &int16_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));
}


/*
 */
static void STDCALL SetFieldValueForInt16ReversedUnaligned(FieldData_pa fd,
                                                           LgIndex_t    pt,
                                                           double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Int16 && !fd->IsNativeOrdering && !fd->IsWordAligned);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      Byte_t *byte_array = GetFieldDataBytePtr(fd) + pt*sizeof(Int16_t);
      Int16_t int16_val = CONVERT_DOUBLE_TO_INT16(val);
      COPY_AND_REVERSE_2_UNALIGNED_BYTES(byte_array, &int16_val);
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));
}


/*
 */
static void STDCALL SetFieldValueForByte(FieldData_pa fd,
                                         LgIndex_t    pt,
                                         double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Byte);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      if (val < 1.0)
        GetFieldDataBytePtr(fd)[pt] = 0;
      else if ( val > 255.0 )
        GetFieldDataBytePtr(fd)[pt] = 255;
      else
        GetFieldDataBytePtr(fd)[pt] = (Byte_t)val;
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Byte_t));
}


/*
 */
static void STDCALL SetFieldValueForBit(FieldData_pa fd,
                                        LgIndex_t    pt,
                                        double       val)
{
  REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(fd->Type==FieldDataType_Bit);
  REQUIRE(0<=pt && pt<GetFieldDataNumValues(fd));
  REQUIRE("val can have any value");

  LgIndex_t ByteOffset = pt/8;
  Byte_t    BitMask    = (01 << (pt%8));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
      if (val < 1.0)
        GetFieldDataBytePtr(fd)[ByteOffset] &= ~BitMask;
      else 
        GetFieldDataBytePtr(fd)[ByteOffset] |= BitMask;
    }

  ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt/8, Byte_t));
}


/**
 */
FieldValueSetFunction_pf DetermineFieldDataSetFunction(FieldDataType_e DataType,
                                                       Boolean_t       IsNativeOrdering,
                                                       Boolean_t       IsWordAligned)
{
  FieldValueSetFunction_pf Result;

  REQUIRE(VALID_FIELD_DATA_TYPE(DataType));
  REQUIRE(VALID_BOOLEAN(IsNativeOrdering));
  REQUIRE(VALID_BOOLEAN(IsWordAligned));

  switch ( DataType )
    {
      case FieldDataType_Float :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = SetFieldValueForFloat;
              else
                Result = SetFieldValueForFloatUnaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = SetFieldValueForFloatReversed;
              else
                Result = SetFieldValueForFloatReversedUnaligned;
            }
        } break;
      case FieldDataType_Double :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = SetFieldValueForDouble;
              else
                Result = SetFieldValueForDoubleUnaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = SetFieldValueForDoubleReversed;
              else
                Result = SetFieldValueForDoubleReversedUnaligned;
            }
        } break;
      case FieldDataType_Int32 :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = SetFieldValueForInt32;
              else
                Result = SetFieldValueForInt32Unaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = SetFieldValueForInt32Reversed;
              else
                Result = SetFieldValueForInt32ReversedUnaligned;
            }
        } break;
      case FieldDataType_Int16 :
        {
          if ( IsNativeOrdering )
            {
              if ( IsWordAligned )
                Result = SetFieldValueForInt16;
              else
                Result = SetFieldValueForInt16Unaligned;
            }
          else
            {
              if ( IsWordAligned )
                Result = SetFieldValueForInt16Reversed;
              else
                Result = SetFieldValueForInt16ReversedUnaligned;
            }
        } break;
      case FieldDataType_Byte :
        {
          /* Bytes cannot be misaligned or reversed */
          Result = SetFieldValueForByte;
        } break;
      case FieldDataType_Bit :
        {
          /* Bits cannot be misaligned or reversed */
          Result = SetFieldValueForBit;
        } break;
      default :
        {
          CHECK(FALSE);
          Result = NULL; /* satisfy compiler */
        } break;
    }
  return (Result);
}


/**
 */
Int64_t FieldDataGetBytesNeeded(LgIndex_t       NumValues,
                                FieldDataType_e DataType)
{
  Int64_t Result = 0; /* ...quite compiler */

  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_FIELD_DATA_TYPE(DataType));

  switch (DataType)
    {
      case FieldDataType_Float:  Result = ((Int64_t)NumValues)*sizeof(float);        break;
      case FieldDataType_Double: Result = ((Int64_t)NumValues)*sizeof(double);       break;
      case FieldDataType_Int32:  Result = ((Int64_t)NumValues)*sizeof(LgIndex_t);    break;
      case FieldDataType_Int16:  Result = ((Int64_t)NumValues)*sizeof(SmInteger_t);  break;
      case FieldDataType_Byte:   Result = ((Int64_t)NumValues)*sizeof(Byte_t);       break;
      case FieldDataType_Bit:    Result = ((Int64_t)(NumValues+7)/8)*sizeof(Byte_t); break;
      default: CHECK(FALSE); break;
    }

  ENSURE(Result >= 0);
  return Result;
}

/**
 * On the SGI, HP, and Sun machines 64 bit objects such as doubles must be 8
 * byte aligned while on all other machines 32 bit alignment suffices. Some
 * allow 1 byte alignment but we won't bother with that.
 */
#if defined IRISX || defined HPUX || defined SUNX
# define SIZEOF_LARGEST_OBJECT_TO_ALIGN sizeof(Int64_t)
#else
# define SIZEOF_LARGEST_OBJECT_TO_ALIGN sizeof(Int32_t)
#endif

/**
 */
Boolean_t IsOffsetAlignedForFieldDataType(FieldDataType_e FieldDataType,
                                          Int64_t         Offset)
{
  REQUIRE(VALID_FIELD_DATA_TYPE(FieldDataType));
  REQUIRE(Offset >= 0);

  Int64_t SizeOfType = FieldDataGetBytesNeeded(1, FieldDataType);
  if (SizeOfType > (Int64_t)SIZEOF_LARGEST_OBJECT_TO_ALIGN)
    SizeOfType = SIZEOF_LARGEST_OBJECT_TO_ALIGN;

  Boolean_t HasValidAlignment = (Offset % SizeOfType == 0);

  ENSURE(VALID_BOOLEAN(HasValidAlignment));
  return HasValidAlignment;
}

/**
 */
Int64_t GetAlignedOffsetForFieldDataType(FieldDataType_e FieldDataType,
                                         Int64_t         Offset)
{
  REQUIRE(VALID_FIELD_DATA_TYPE(FieldDataType));
  REQUIRE(Offset >= 0);

  Int64_t SizeOfType = FieldDataGetBytesNeeded(1, FieldDataType);
  if (SizeOfType > (Int64_t)SIZEOF_LARGEST_OBJECT_TO_ALIGN)
    SizeOfType = SIZEOF_LARGEST_OBJECT_TO_ALIGN;

  Int64_t NumBytesPastAlignment = (Offset % SizeOfType);
  Int64_t Result = Offset - NumBytesPastAlignment;

  ENSURE(0 <= Result && Result <= Offset);
  ENSURE(IsOffsetAlignedForFieldDataType(FieldDataType, Result));
  return Result;
}

/**
 */
void FieldDataDefineData(FieldData_pa    FieldData,
                         LgIndex_t       NumValues,
                         FieldDataType_e DataType,
                         ValueLocation_e ValueLocation)
{
  REQUIRE(VALID_REF(FieldData));
  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_FIELD_DATA_TYPE(DataType));
  REQUIRE(VALID_ENUM(ValueLocation, ValueLocation_e));

  /*
   * Remove any old data (transformed UVW is one example that calls this
   * function with a non-null data pointer when switching the value location
   * when style changes the value location and therefore the amount of data
   * allocated.)
   */
  FieldDataCleanup(FieldData);
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif

  /*
   * The reference count is not modified here. This function only allocates the
   * structure and makes adjustments to the some of the members. The reference
   * count was initialized when the structure was initially created and the
   * structure may be shared before the data portion is even allocated.
   */
  FieldData->NumValues        = NumValues;
  FieldData->Type             = DataType;
  FieldData->IsNativeOrdering = TRUE;

  FieldData->IsWordAligned    = TRUE;
  FieldData->ValueLocation    = ValueLocation;
  FieldData->GetFunction      = DetermineFieldDataGetFunction(DataType,
                                                              FieldData->IsNativeOrdering,
                                                              FieldData->IsWordAligned);
  FieldData->SetFunction      = DetermineFieldDataSetFunction(DataType,
                                                              FieldData->IsNativeOrdering,
                                                              FieldData->IsWordAligned);

# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif

  ENSURE(FieldData->Data == NULL);
}

/**
 */
Boolean_t FieldDataAllocData(FieldData_pa FieldData,
                             Boolean_t    ShowErrMsg)
{
  Boolean_t IsOk;
  Int64_t   ActualBytesNeeded;
  size_t    BytesToAllocate;

  REQUIRE(VALID_REF(FieldData));
  REQUIRE(FieldData->Type != FieldDataType_Invalid); /* ...must call FieldDataDefineData first */
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  REQUIRE(FieldData->Data == NULL);
  REQUIRE(FieldData->IsNativeOrdering);
  REQUIRE(FieldData->IsWordAligned);
  REQUIRE(VALID_BOOLEAN(ShowErrMsg));

  /*
   * The size of size_t may be smaller than our unsigned 64 bit integer value
   * so we might have to squeeze it down possibly loosing precision.
   */
  ActualBytesNeeded = FieldDataGetBytesNeeded(FieldData->NumValues,
                                              FieldData->Type);
  BytesToAllocate = (size_t)ActualBytesNeeded;

  /*
   * 64 bit architectures are effectively unlimited in their allocation size
   * while 32 architectures are limited to 4GB (some may limit further to 2GB
   * which will be borne out by the call to malloc).
   */
  CHECK(sizeof(size_t) == 4 || sizeof(size_t) == 8);
  IsOk = (sizeof(size_t) == 8 ||
          ActualBytesNeeded <= (Int64_t)0xffffffff);
  if (IsOk)
    {
      if (BytesToAllocate == 0)
        FieldData->Data = NULL;
      else
        {
          FieldData->Data = (void *)ALLOC_ARRAY(BytesToAllocate,
                                                char,
                                                "FieldData's Data");
          #if defined DEBUG_FIELDVALUES
            {
              if (FieldDataHasContiguousAddressSpace(FieldData))
                memset(FieldData->Data, DEBUG_FIELDVALUES_BAD_VALUE, BytesToAllocate);
            }
          #endif
        }
      IsOk = (FieldData->NumValues == 0 ||
              FieldData->Data != NULL);
      if (!IsOk && ShowErrMsg)
        OutOfMemoryMsg();
    }
  else if (ShowErrMsg)
    ErrMsg(_ts("Storage limit (%ld) exceeded for a single variable."),(long)MAXINDEX);

# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif

  ENSURE(VALID_REF(FieldData->Data) || FieldData->Data == NULL);
  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if !defined NO_ASSERTS
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 * Allocates a field data pointer with space for "num_pts" of field data type
 * "field_data_type" nodal values.
 *
 * IMPORTANT:
 *   This field data may NOT be used for zones but only for things like
 *   geometries or other temporary field data that will never be placed
 *   into a COB or zone.
 */
FieldData_pa AllocScratchNodalFieldDataPtr(LgIndex_t       NumValues,
                                           FieldDataType_e Type,
                                           Boolean_t       ShowErrMsg)
{
  FieldData_pa Result;

  REQUIRE(NumValues >= 0);
  REQUIRE(VALID_FIELD_DATA_TYPE(Type));
  REQUIRE(VALID_BOOLEAN(ShowErrMsg));

  Result = FieldDataAlloc();
  if (Result != NULL)
    {
      FieldDataDefineData(Result, NumValues, Type, ValueLocation_Nodal);
      if (!FieldDataAllocData(Result, ShowErrMsg))
        FieldDataDealloc(&Result, FALSE);
    }
  else if (ShowErrMsg)
    OutOfMemoryMsg();

  ENSURE(VALID_REF(Result) || Result == NULL);
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
  ENSURE(IMPLICATION(Result != NULL, VALID_FN_REF(Result->GetFunction)));
  ENSURE(IMPLICATION(Result != NULL, VALID_FN_REF(Result->SetFunction)));
  ENSURE(IMPLICATION(Result != NULL,
                     (Result->NumValues >= 0                &&
                      IMPLICATION(Result->NumValues != 0,
                                  VALID_REF(Result->Data))  &&
                      VALID_FIELD_DATA_TYPE(Result->Type))));

  return Result;
}


/**
 * Frees memory allocated with AllocScratchNodalFieldDataPtr().
 *
 * @param ScratchFieldData
 *   Scratch field data pointer to deallocate. This should NEVER be a field
 *   data from a zone or COB. See note in AllocScratchNodalFieldDataPtr().
 */
void DeallocScratchNodalFieldDataPtr(FieldData_pa *FieldDataRef)
{
  FieldDataDealloc(FieldDataRef,
                   FALSE); /* DoTrackVarSharing */
}


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
FieldDataType_e GetFieldDataType_FUNC(FieldData_pa FieldData)
{
  FieldDataType_e Result;

  REQUIRE(VALID_REF(FieldData));

  Result = GetFieldDataType_MACRO(FieldData);

  ENSURE(VALID_FIELD_DATA_TYPE(Result));
  return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
FieldValueGetFunction_pf GetFieldDataGetFunction_FUNC(FieldData_pa FieldData)
{
  FieldValueGetFunction_pf Result;

  REQUIRE(VALID_REF(FieldData));

  Result = GetFieldDataGetFunction_MACRO(FieldData);

  ENSURE(VALID_FN_REF(Result));
  return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
FieldValueSetFunction_pf GetFieldDataSetFunction_FUNC(FieldData_pa FieldData)
{
  FieldValueSetFunction_pf Result;

  REQUIRE(VALID_REF(FieldData));

  Result = GetFieldDataSetFunction_MACRO(FieldData);

  ENSURE(VALID_FN_REF(Result));
  return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
LgIndex_t GetFieldDataNumValues_FUNC(FieldData_pa FieldData)
{
  LgIndex_t Result;

  REQUIRE(VALID_REF(FieldData));

  Result = GetFieldDataNumValues_MACRO(FieldData);

  ENSURE(Result >= 0);
  return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
ValueLocation_e GetFieldDataValueLocation_FUNC(FieldData_pa FieldData)
{
  ValueLocation_e Result;

  REQUIRE(VALID_REF(FieldData));

  Result = GetFieldDataValueLocation_MACRO(FieldData);

  ENSURE(Result == ValueLocation_Invalid || /* i.e. pending assignment */
         VALID_ENUM(Result, ValueLocation_e));
  return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
Boolean_t FieldDataHasContiguousAddressSpace_FUNC(FieldData_pa FieldData)
{
  REQUIRE(VALID_REF(FieldData));
  Boolean_t Result = FieldDataHasContiguousAddressSpace_MACRO(FieldData);
  ENSURE(VALID_BOOLEAN(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
float *GetFieldDataFloatPtr_FUNC(FieldData_pa fd)
{
  float *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Float && fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataFloatPtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
FloatRev_t *GetFieldDataFloatRevPtr_FUNC(FieldData_pa fd)
{
  FloatRev_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Float && !fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataFloatRevPtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
double *GetFieldDataDoublePtr_FUNC(FieldData_pa fd)
{
  double *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Double && fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataDoublePtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
DoubleRev_t *GetFieldDataDoubleRevPtr_FUNC(FieldData_pa fd)
{
  DoubleRev_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Double && !fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataDoubleRevPtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
Int64_t *GetFieldDataInt64Ptr_FUNC(FieldData_pa fd)
{
  Int64_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Int64 && fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataInt64Ptr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
Int64Rev_t *GetFieldDataInt64RevPtr_FUNC(FieldData_pa fd)
{
  Int64Rev_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Int64 && !fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataInt64RevPtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
Int32_t *GetFieldDataInt32Ptr_FUNC(FieldData_pa fd)
{
  Int32_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Int32 && fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataInt32Ptr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
Int32Rev_t *GetFieldDataInt32RevPtr_FUNC(FieldData_pa fd)
{
  Int32Rev_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Int32 && !fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataInt32RevPtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
Int16_t *GetFieldDataInt16Ptr_FUNC(FieldData_pa fd)
{
  Int16_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Int16 && fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataInt16Ptr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 */
Int16Rev_t *GetFieldDataInt16RevPtr_FUNC(FieldData_pa fd)
{
  Int16Rev_t *Result;
  REQUIRE(VALID_REF(fd) && fd->Type==FieldDataType_Int16 && !fd->IsNativeOrdering && fd->IsWordAligned);
  Result = GetFieldDataInt16RevPtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * No byte ordering or alignment issues with single bytes (which are also used with the "Bit" type)
 */
Byte_t *GetFieldDataBytePtr_FUNC(FieldData_pa fd)
{
  Byte_t *Result;
  /*
   * This function gets called for Byte and Bit types, but we cannot REQUIRE
   * those types because it is also used for non-aligned values.  We can't
   * check for non-aligned because we might be copying aligned bytes to a
   * non-aligned location.
   */
  REQUIRE(VALID_REF(fd));
  Result = GetFieldDataBytePtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * Gets a ptr to 2-byte blocks regardless of byte ordering, but still has to
 * worry about byte alignment
 */
UInt16_t *GetFieldData2BytePtr_FUNC(FieldData_pa fd)
{
  UInt16_t *Result;
  REQUIRE(VALID_REF(fd) && fd->IsWordAligned);
  REQUIRE(fd->Type==FieldDataType_Int16);
  Result = GetFieldData2BytePtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * Gets a ptr to 4-byte blocks regardless of byte ordering, but still has to
 * worry about byte alignment
 */
UInt32_t *GetFieldData4BytePtr_FUNC(FieldData_pa fd)
{
  UInt32_t *Result;
  REQUIRE(VALID_REF(fd) && fd->IsWordAligned);
  REQUIRE(fd->Type==FieldDataType_Int32 || fd->Type==FieldDataType_Float);
  Result = GetFieldData4BytePtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * Gets a ptr to 8-byte blocks regardless of byte ordering, but still has to
 * worry about byte alignment
 */
UInt64_t *GetFieldData8BytePtr_FUNC(FieldData_pa fd)
{
  UInt64_t *Result;
  REQUIRE(VALID_REF(fd) && fd->IsWordAligned);
  REQUIRE(fd->Type==FieldDataType_Int64 || fd->Type==FieldDataType_Double);
  Result = GetFieldData8BytePtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * WARNING: GetFieldDataVoidPtr checks nothing, and thus should only be
 * used with extreme caution (that is, checking the alignment
 * and byte order by hand).
 */
void *GetFieldDataVoidPtr_FUNC(FieldData_pa fd)
{
  void *Result;
  REQUIRE(VALID_REF(fd));
  Result = GetFieldDataVoidPtr_MACRO(fd);
  ENSURE(VALID_REF(Result));
  return Result;
}
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/*
 */
void CopyFieldValue(FieldData_pa  dst,
                    LgIndex_t     dstindex,
                    FieldData_pa  src,
                    LgIndex_t     srcindex)
{
  REQUIRE(VALID_REF(dst));
  REQUIRE(VALID_REF(src));
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
  REQUIRE( dstindex >= 0 && dstindex < GetFieldDataNumValues(dst) &&
           srcindex >= 0 && srcindex < GetFieldDataNumValues(src) );
  
  Boolean_t DoBruteForceCopy = TRUE;

  if ( FieldDataHasContiguousAddressSpace(src) &&
       FieldDataHasContiguousAddressSpace(dst) &&
       GetFieldDataType(src) == GetFieldDataType(dst) )
    {
      switch ( GetFieldDataType(src) )
        {
          case FieldDataType_Int64 : CHECK(FALSE); /* Future work: remove and let fall through */
          case FieldDataType_Double :
            {
              CHECK(sizeof(UInt64_t)==8 && sizeof(double)==8);
              if ( src->IsWordAligned && dst->IsWordAligned )
                {
                  UInt64_t *dst_ptr = GetFieldData8BytePtr(dst) + dstindex;
                  UInt64_t *src_ptr = GetFieldData8BytePtr(src) + srcindex;
                  *dst_ptr = *src_ptr;
                  if ( src->IsNativeOrdering != dst->IsNativeOrdering )
                    REVERSE_8_BYTES(dst_ptr);
                }
              else
                {
                  Byte_t *dst_ptr = GetFieldDataBytePtr(dst) + dstindex*sizeof(UInt64_t);
                  Byte_t *src_ptr = GetFieldDataBytePtr(src) + srcindex*sizeof(UInt64_t);
                  if ( src->IsNativeOrdering == dst->IsNativeOrdering )
                    COPY_8_UNALIGNED_BYTES(dst_ptr, src_ptr);
                  else
                    COPY_AND_REVERSE_8_UNALIGNED_BYTES(dst_ptr, src_ptr);
                }
              DoBruteForceCopy = FALSE;
            } break;
          case FieldDataType_Float :
          case FieldDataType_Int32 :
            {
              CHECK(sizeof(UInt32_t)==4 && sizeof(float)==4);
              if ( src->IsWordAligned && dst->IsWordAligned )
                {
                  UInt32_t *dst_ptr = GetFieldData4BytePtr(dst) + dstindex;
                  UInt32_t *src_ptr = GetFieldData4BytePtr(src) + srcindex;
                  *dst_ptr = *src_ptr;
                  if ( src->IsNativeOrdering != dst->IsNativeOrdering )
                    REVERSE_4_BYTES(dst_ptr);
                }
              else
                {
                  Byte_t *dst_ptr = GetFieldDataBytePtr(dst) + dstindex*sizeof(UInt32_t);
                  Byte_t *src_ptr = GetFieldDataBytePtr(src) + srcindex*sizeof(UInt32_t);
                  if ( src->IsNativeOrdering == dst->IsNativeOrdering )
                    COPY_4_UNALIGNED_BYTES(dst_ptr, src_ptr);
                  else
                    COPY_AND_REVERSE_4_UNALIGNED_BYTES(dst_ptr, src_ptr);
                }
              DoBruteForceCopy = FALSE;
            } break;
          case FieldDataType_Int16 :
            {
              CHECK(sizeof(UInt16_t)==2);
              if ( src->IsWordAligned && dst->IsWordAligned )
                {
                  UInt16_t *dst_ptr = GetFieldData2BytePtr(dst) + dstindex;
                  UInt16_t *src_ptr = GetFieldData2BytePtr(src) + srcindex;
                  *dst_ptr = *src_ptr;
                  if ( src->IsNativeOrdering != dst->IsNativeOrdering )
                    REVERSE_2_BYTES(dst_ptr);
                }
              else
                {
                  Byte_t *dst_ptr = GetFieldDataBytePtr(dst) + dstindex*sizeof(UInt16_t);
                  Byte_t *src_ptr = GetFieldDataBytePtr(src) + srcindex*sizeof(UInt16_t);
                  if ( src->IsNativeOrdering == dst->IsNativeOrdering )
                    COPY_2_UNALIGNED_BYTES(dst_ptr, src_ptr);
                  else
                    COPY_AND_REVERSE_2_UNALIGNED_BYTES(dst_ptr, src_ptr);
                }
            } break;
          case FieldDataType_Byte :
            {
              /* no byte ordering or alignment worries */
              GetFieldDataBytePtr(dst)[dstindex] = GetFieldDataBytePtr(src)[srcindex];
              DoBruteForceCopy = FALSE;
            } break;
          case FieldDataType_Bit : break; /* handle below */
          default : CHECK(FALSE); /* Future work: when more complex types are added, remove this CHECK */
        }
    }

  if ( DoBruteForceCopy )
    {
      double val = GetFieldValue(src,srcindex);
      SetFieldValue(dst, dstindex, val);
    }
} /* CopyFieldValue() */


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */


/*
 */
void SetFieldDataPtrToAllZeros(FieldData_pa fd)
{
  REQUIRE(VALID_REF(fd));
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

  LgIndex_t NumValues = GetFieldDataNumValues(fd);

  /*
   * memset each byte to 0 works for floats and doubles and works regardless
   * of byte ordering or alignment.
   */
  size_t NumBytesToMemSet = 0;
  if (FieldDataHasContiguousAddressSpace(fd))
    {
      switch ( GetFieldDataType(fd) )
        {
          case FieldDataType_Int64 : CHECK(FALSE); /* Future work: remove CHECK */
          case FieldDataType_Double :
            {
              CHECK(sizeof(UInt64_t)==8 && sizeof(double)==8);
              NumBytesToMemSet = NumValues*sizeof(Int64_t);
            } break;
          case FieldDataType_Int32 :
          case FieldDataType_Float :
            {
              CHECK(sizeof(UInt32_t)==4 && sizeof(float)==4);
              NumBytesToMemSet = NumValues*sizeof(Int32_t);
            } break;
          case FieldDataType_Int16 :
            {
              CHECK(sizeof(UInt16_t)==2);
              NumBytesToMemSet = NumValues*sizeof(Int16_t);
            } break;
          case FieldDataType_Byte :
            {
              NumBytesToMemSet = NumValues*sizeof(Byte_t);
            } break;
          case FieldDataType_Bit :
            {
              NumBytesToMemSet = ((NumValues+7)/8)*sizeof(Byte_t);
            } break;
          default :
            {
              CHECK(FALSE);
            } break;
        }
    }

  if ( NumBytesToMemSet > 0 )
    {
      void *fd_data = GetFieldDataVoidPtr(fd);
      memset(fd_data, 0, NumBytesToMemSet);
    }
  else
    {
      int ii;
      for ( ii = 0; ii < NumValues; ii++ )
        SetFieldValue(fd, ii, 0.0);
    }

} /* SetFieldDataPtrToAllZeros() */



#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif



#if defined USEFUNCTIONSFORNODEVALUES
/*
 */
NodeMap_t GetNodeValue(NodeMap_pa  node_map,
                       LgIndex_t   element,
                       LgIndex_t   corner)
{
  LgIndex_t Result;
  LgIndex_t pos;
  REQUIRE(VALID_REF(node_map));
  REQUIRE(element < GetNMNumElements(node_map));
  REQUIRE(corner < GetNMPtsPerElement(node_map));

  pos = element*GetNMPtsPerElement(node_map) + corner;
  Result = GetNMPtr(node_map)[pos];

  ENSURE(Result>=0);
  return Result;
} /* GetNodeValue() */


/*
 */
NodeMap_t GetNodeValueAtOffset(NodeMap_pa  node_map,
                               LgIndex_t   offset)
{
  LgIndex_t Result;
  REQUIRE(VALID_REF(node_map));
  REQUIRE(offset < GetNMNumElements(node_map)*GetNMPtsPerElement(node_map));

  Result = GetNMPtr(node_map)[offset];

  ENSURE(Result>=0);
  return Result;
} /* GetNodeValueAtOffset() */


/*
 */
void SetNodeValue(NodeMap_pa  node_map,
                  LgIndex_t   element,
                  LgIndex_t   corner,
                  NodeMap_t   node)
{
  LgIndex_t pos;

  REQUIRE(VALID_REF(node_map));
  REQUIRE(element < GetNMNumElements(node_map));
  REQUIRE(corner < GetNMPtsPerElement(node_map));
  REQUIRE(node >= 0);

  pos = element*GetNMPtsPerElement(node_map) + corner;
  GetNMPtr(node_map)[pos] = node;

} /* SetNodeValue() */


/*
 */
void SetNodeValueAtOffset(NodeMap_pa  node_map,
                          LgIndex_t   offset,
                          NodeMap_t   node)
{
  REQUIRE(VALID_REF(node_map));
  REQUIRE(offset < GetNMNumElements(node_map)*GetNMPtsPerElement(node_map));
  REQUIRE(node >= 0);

  GetNMPtr(node_map)[offset] = node;

} /* SetNodeValueAtOffset() */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
