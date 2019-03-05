/***************************************************************************
 *
 *   vector.cpp -- Reference counted vector class
 *
 *   Erik Persson and Björn Nilsson, 2000-
 */

#include "types/vector.h"

// A global empty CVectorData. 
static int g_NilVectorDataInts[] = { 0x40000000, 0, 0, 0 };

// Pointer into the empty CVectorData, used by all empty vectors
void *g_pNilVectorElements= ((CVectorData*) &g_NilVectorDataInts)+1;

