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
#if defined TEXTMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

/* These macros for checking CoordSys_e and Units of text objects (i.e., those associated with the text tool). */
#define VALID_TEXT_COORDSYS(sys)  (((sys)==CoordSys_Frame)||((sys)==CoordSys_Grid)||((sys)==CoordSys_Grid3D))
#define VALID_TEXT_UNITS(units)  (((units)==Units_Grid)||((units)==Units_Frame)||((units)==Units_Point))
#define VALID_TEXT_COORDSYS_AND_UNITS(pos_sys, size_units) \
           ( VALID_TEXT_COORDSYS((pos_sys)) && \
             VALID_TEXT_UNITS((size_units)) && \
             ! ((pos_sys) == CoordSys_Frame && (size_units) == Units_Grid) )

/* This is for any type of font in Tecplot. */
#define VALID_FONT_SIZEUNITS(units)  (((units)==Units_Grid)||((units)==Units_Frame)||((units)==Units_Point)||(units)==Units_AxisPercentage)

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
