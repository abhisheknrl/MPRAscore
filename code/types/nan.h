/*******************************************************************************
 *
 *   floatnan.h -- Helper functions for checking float validity
 *
 *   Björn Nilsson, 2006
 */

#ifndef FLOAT_NaN
#define FLOAT_NaN (0x7f800000)

/*
// deprecated 2017-08-29, use double. /BN

inline bool IsN(float x) {	return ((*(long *)&x) & 0x7f800000)!=0x7f800000; }
inline bool IsNaN(float x) { return ((*(long *)&x) & 0x7f800000)==0x7f800000; }
*/

#endif
