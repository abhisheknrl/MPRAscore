/*******************************************************************************
 *
 *   assert.h -- Assert wrapper
 *
 *   Björn Nilsson, 2004-2008
 */

#ifndef SYSTEM_ASSERT_H
#define SYSTEM_ASSERT_H
#ifdef MFC_CLASSES
// MFC in use
#include <afx.h>
#else
// MFC not in use
#ifdef _DEBUG
#include <assert.h>
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)
#endif // _DEBUG
#endif // MFC_CLASSES
#endif // SYSTEM_ASSERT_H
