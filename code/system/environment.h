/***********************************************************************************
 *
 *   environment.h -- Interface to global environment variables
 *
 *   Erik Persson and Björn Nilsson, 2000-
 */

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "types/string.h"
#ifdef _WIN32
#ifndef MFC_CLASSES
#include <windows.h> // For LARGE_INTEGER
#endif
#endif

namespace Environment
{
	CString GetUserName();
	CString GetUserDir();
	CString GetEnv(const char *aName);
#ifdef _WIN32
	unsigned int GetTimeDiff(LARGE_INTEGER c0,LARGE_INTEGER c1);
#endif
};

#endif
