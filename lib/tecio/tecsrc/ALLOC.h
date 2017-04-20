/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2004 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/
#if !defined ALLOC_H__
#define ALLOC_H__

#if defined EXTERN
  #undef EXTERN
#endif
#if defined ALLOCMODULE
  #define EXTERN
#else
  #define EXTERN extern
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#ifdef NO_DEBUG_ALLOC
  #undef DEBUG_ALLOC
#endif

#if defined DEBUG_ALLOC

  #if defined CRUDEALLOCVERIFY
    size_t GetArraySize(void *Adr);
    void CheckActiveAllocs(void);
  #endif

EXTERN void *newmalloc(size_t      size,
                       const char *file,
                       long        line,
                       const char *desc);

EXTERN void newfree(void       *ptr,
                    const char *file,
                    long        line,
                    const char *desc);

  #define ALLOC_ARRAY(N,Type,S) (Type *)newmalloc((N)*sizeof(Type),\
                                           __FILE__,__LINE__,(S))
  #define ALLOC_ITEM(Type,S)    (Type *)newmalloc(sizeof(Type),\
                                           __FILE__,__LINE__,(S))

  /* NOTE: the pointer is set to 0xFFFF after the free for debug   */
  /*       versions in the hopes of catching invalid pointer usage */
  #define FREE_ARRAY(X,S) (newfree((void *)(X), __FILE__,__LINE__,(S)), \
                          (*((void **)&(X))=(void *)0xFFFF))
  #define FREE_ITEM(X,S)  (newfree((void *)(X), __FILE__,__LINE__,(S)), \
                          (*((void **)&(X))=(void *)0xFFFF))

#elif defined UNIXX || defined DOS || defined VMS

  #if defined TECPLOTKERNEL
    EXTERN void *newmalloc(size_t size);
    EXTERN void newfree(void *ptr);

    #define ALLOC_ARRAY(N,Type,S) (Type *)newmalloc((N)*sizeof(Type))
    #define ALLOC_ITEM(Type,S)    (Type *)newmalloc(sizeof(Type))
    #define FREE_ARRAY(X,S)       newfree((void *)(X))
    #define FREE_ITEM(X,S)        newfree((void *)(X))
  #else
    #define ALLOC_ARRAY(N,Type,S) (Type *)malloc((N)*sizeof(Type))
    #define ALLOC_ITEM(Type,S)    (Type *)malloc(sizeof(Type))
    #define FREE_ARRAY(X,S)       free((void *)(X))
    #define FREE_ITEM(X,S)        free((void *)(X))
  #endif

#elif defined MSWIN

  #if !defined (__cplusplus)
    #define ALLOC_MALLOCFREE /* new/delete not available */
  #else
    #define ALLOC_NEWDELETE
  #endif /* __cplusplus */

  #ifdef ALLOC_NEWDELETE
    #define ALLOC_ARRAY(N,Type,S) new Type[N]
    #define ALLOC_ITEM(Type,S)    new Type
    #ifdef _DEBUG
      /* NOTE: the pointer is set to 0xFFFF after the free for debug   */
      /*       versions in the hopes of catching invalid pointer usage */
      #define FREE_ARRAY(X,S)  do { CHECK(AfxIsValidAddress((X),1)); delete [] (X); *((void **)&(X)) = (void *)0xFFFF; } while (0)
      #define FREE_ITEM(X,S)   do { CHECK(AfxIsValidAddress((X),1)); delete (X); *((void **)&(X)) = (void *)0xFFFF; } while (0)
    #else
      #define FREE_ARRAY(X,S)  delete [] (X)
      #define FREE_ITEM(X,S)   delete (X)
    #endif
  #endif

  #ifdef ALLOC_MALLOCFREE
    #define ALLOC_ARRAY(N,Type,S) (Type *)malloc((N)*sizeof(Type))
    #define ALLOC_ITEM(Type,S)    (Type *)malloc(sizeof(Type))
    #define FREE_ARRAY(X,S)       free((void *)(X))
    #define FREE_ITEM(X,S)        free((void *)(X))
  #endif

  #ifdef ALLOC_HEAP
    #define ALLOC_ARRAY(N,Type,S) (Type *)HeapAlloc(GetProcessHeap(),\
                                             NULL,(N)*sizeof(Type))
    #define ALLOC_ITEM(Type,S)    (Type *)HeapAlloc(GetProcessHeap(),\
                                             NULL,sizeof(Type))
    #define FREE_ARRAY(X,S)       HeapFree(GetProcessHeap(),NULL,(void *)(X))
    #define FREE_ITEM(X,S)        HeapFree(GetProcessHeap(),NULL,(void *)(X))
  #endif

  #ifdef ALLOC_MIXED
    void *MSWinAlloc(DWORD nSize);
    void MSWinFree(void *pMem);
    #define ALLOC_ARRAY(N,Type,S) (Type *)MSWinAlloc((N)*sizeof(Type))
    #define ALLOC_ITEM(Type,S)    (Type *)MSWinAlloc(sizeof(Type))
    #define FREE_ARRAY(X,S)       MSWinFree((void *)(X))
    #define FREE_ITEM(X,S)        MSWinFree((void *)(X))
  #endif

  #if defined MSWIN_DEBUG_ALLOC
    EXTERN unsigned char *mswin_debug_alloc(int N,
                                            const char *Type,
                                            const char *file,
                                            int Line,
                                            const char *desc);

    EXTERN void mswin_debug_free(unsigned char *X,
                                 const char *file,
                                 int Line,
                                 const char *desc);


    #undef ALLOC_ARRAY
    #undef FREE_ARRAY
    
    #define ALLOC_ARRAY(N,Type,S) (Type *) mswin_debug_alloc((N)*sizeof(Type),#Type,__FILE__,__LINE__,S)
    #define FREE_ARRAY(X,S)                mswin_debug_free((unsigned char*)(X),__FILE__,__LINE__,S)

  #endif /* MSWIN_DEBUG_ALLOC */

#endif

#endif /* ALLOC_H__ */





