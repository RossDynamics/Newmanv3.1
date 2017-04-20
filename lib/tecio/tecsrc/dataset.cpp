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

#define DATASETMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "STRUTIL.h"
#include "AUXDATA.h"
#include "ARRLIST.h"
#include "STRLIST.h"
#include "ALLOC.h"
#include "SET.h"
#include "DATASET.h"
#include "FILESTREAM.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined USEPAGES
#endif
#endif
#include "DATASET0.h"

#include <float.h>


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined USE_MACROS_FOR_FUNCTIONS
#endif
#if !defined USE_MACROS_FOR_FUNCTIONS
#endif
#endif /* TECPLOTKERNEL */

/**
 * Cleanout the contents of the zone spec item but leaves the zone spec item.
 * This effectively leaves the zone spec structure in the same state as calling
 * ZoneSpecAlloc initially.
 *
 * param ZoneSpec
 *     Zone spec item to cleanup.
 */
void CleanoutZoneSpec(ZoneSpec_s *ZoneSpec)
{
  REQUIRE(VALID_REF(ZoneSpec));

  if (ZoneSpec->Name != NULL)
    FREE_ARRAY(ZoneSpec->Name, "ZoneSpec name");
  if (ZoneSpec->AuxData != NULL)
    AuxDataDealloc(&ZoneSpec->AuxData);
  SetZoneSpecDefaults(ZoneSpec);
}


/**
 */
void ZoneSpecDealloc(ZoneSpec_s **ZoneSpec)
{
  REQUIRE(VALID_REF(ZoneSpec));
  REQUIRE(VALID_REF(*ZoneSpec) || *ZoneSpec == NULL);

  if (*ZoneSpec != NULL)
    {
      CleanoutZoneSpec(*ZoneSpec);

      FREE_ITEM(*ZoneSpec, "ZoneSpec structure");
      *ZoneSpec = NULL;
    }

  ENSURE(*ZoneSpec == NULL);
}

/**
 */
Boolean_t ZoneSpecItemDestructor(void       *ItemRef,
                                 ArbParam_t ClientData)
{
  ZoneSpec_s **ZoneSpecRef = (ZoneSpec_s **)ItemRef;

  REQUIRE(VALID_REF(ZoneSpecRef));
  REQUIRE(VALID_REF(*ZoneSpecRef) || *ZoneSpecRef == NULL);

  if (*ZoneSpecRef != NULL)
    ZoneSpecDealloc(ZoneSpecRef);

  ENSURE(*ZoneSpecRef == NULL);
  return TRUE;
}

/**
 */
void SetZoneSpecDefaults(ZoneSpec_s *ZoneSpec)
{
  REQUIRE(VALID_REF(ZoneSpec));
  ZoneSpec->Name                         = NULL;
  ZoneSpec->UniqueID                     = INVALID_UNIQUE_ID;
  ZoneSpec->ParentZone                   = BAD_SET_VALUE;
  ZoneSpec->StrandID                     = STRAND_ID_STATIC;
  ZoneSpec->SolutionTime                 = 0.0;
  ZoneSpec->NumPtsI                      = 0;
  ZoneSpec->NumPtsJ                      = 0;
  ZoneSpec->NumPtsK                      = 0;
  ZoneSpec->ICellDim                     = 0;
  ZoneSpec->JCellDim                     = 0;
  ZoneSpec->KCellDim                     = 0;
  ZoneSpec->Type                         = ZoneType_Ordered;
  ZoneSpec->ZoneLoadInfo.PresetZoneColor = NoColor_C;
  ZoneSpec->ZoneLoadInfo.IsInBlockFormat = TRUE;
  ZoneSpec->UserDefFNConnections         = 0;
  ZoneSpec->AutoAssignFN                 = TRUE;
  ZoneSpec->FNMode                       = FaceNeighborMode_LocalOneToOne;
  ZoneSpec->AuxData                      = NULL;
  ZoneSpec->BuildZoneOptInfo             = TRUE;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

/**
 */
ZoneSpec_s *ZoneSpecAlloc(void)
{
  ZoneSpec_s *Result;
  
  Result = (ZoneSpec_s *)ALLOC_ITEM(ZoneSpec_s, "ZoneSpec structure");
  if (Result != NULL)
    SetZoneSpecDefaults(Result);

  ENSURE(Result == NULL || VALID_REF(Result));
  return Result;
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

/**
 * Adjusts the capacity request as necessary to minimize memory reallocations
 * for large lists. The adjusted capacity will be at least as big as requested
 * however it may be larger if it is determined that the space requirement is
 * growing faster.
 *
 * param ZoneOrVarArrayList
 *     Array list requesting the change in capacity.
 * param CurrentCapacity
 *     Current capacity of the array list.
 * param RequestedCapacity
 *     Capacity request or zero for default size.
 * param ClientData
 *     Any client data needed for the adjustment.
 *
 * return
 *     Adjusted capacity that is at least as large as the request or zero if
 *     unable to satisfy the requested capacity.
 */
LgIndex_t ZoneOrVarListAdjustCapacityRequest(ArrayList_pa ZoneOrVarArrayList,
                                             LgIndex_t    CurrentCapacity,
                                             LgIndex_t    RequestedCapacity,
                                             ArbParam_t   ClientData)
{
  LgIndex_t Result;

  REQUIRE(ArrayListIsValid(ZoneOrVarArrayList));
  REQUIRE((RequestedCapacity == 0 && CurrentCapacity == 0) ||
          RequestedCapacity > CurrentCapacity);
  REQUIRE(CurrentCapacity < MaxNumZonesOrVars);

  if (RequestedCapacity <= MaxNumZonesOrVars)
    {
      if (RequestedCapacity != 0 && CurrentCapacity == 0)
        {
          /* first allocation; assume the request is the desired capacityy */
          Result = RequestedCapacity;
        }
      else
        {
          const LgIndex_t DEFAULT_CAPACITY = 32;
          LgIndex_t       BlockSize = MAX(DEFAULT_CAPACITY, CurrentCapacity / 2);
          if (RequestedCapacity == 0)
            Result = DEFAULT_CAPACITY;
          else
            Result = ((RequestedCapacity - 1) / BlockSize + 1) * BlockSize;

          /* put a cap on the maximum */
          if (Result > MaxNumZonesOrVars)
            Result = MaxNumZonesOrVars;
        }
    }
  else
    Result = 0; /* request exceeded maximum; unable to satisfy request */

  ENSURE(Result == 0 || Result >= RequestedCapacity);
  ENSURE(Result <= MaxNumZonesOrVars);
  return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#if !defined USE_MACROS_FOR_FUNCTIONS
#endif
#         if defined DEBUGUNIQUE
#         endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#else
#endif
#if defined USEPAGES
#endif
#endif /* TECPLOTKERNEL */
