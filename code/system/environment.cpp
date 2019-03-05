/***********************************************************************************
 *
 *   environment.cpp -- Interface to global environment variables
 *
 *   Björn Nilsson, 2000-
 */

#ifdef _WIN32
#ifdef MFC_CLASSES
#include <afx.h>
#else
#include <windows.h>
#endif
#endif

#include "environment.h"
#include <cstdint> // for int64_t

#define STRBUFSIZE (0x7fff)
using namespace Environment;

/***********************************************************************************/

CString Environment::GetEnv(const char *aName)
{
	// TODO: UNIX: Implement
#ifdef _WIN32
	char aBuf[STRBUFSIZE];
	int n= ::GetEnvironmentVariable(aName,aBuf,STRBUFSIZE); //Win32 API

	// Fit inside buffer?
	if (n>=0 && n < STRBUFSIZE)
		return CString(aBuf,n);
#endif
	return CString("");
}

/***********************************************************************************/

CString Environment::GetUserName()
{
#ifdef _WIN32
	// Windows
	char aBuf[STRBUFSIZE];
	DWORD n= STRBUFSIZE;
	if (!::GetUserName(aBuf,&n))
		n= -1;

	// Fit inside buffer?
	if (n>=0 && n < STRBUFSIZE)
		return CString(aBuf,n);

	return CString("");
#else
	// Unix
	return GetEnv("USER");
#endif
}

#ifdef _WIN32
unsigned int Environment::GetTimeDiff(LARGE_INTEGER c0,LARGE_INTEGER c1)
{
		LARGE_INTEGER f;
		if(QueryPerformanceFrequency(&f))
		{
			//	LARGE_INTEGER c;
			//	QueryPerformanceCounter(&c);
			int64_t test=
				 (int64_t(c1.HighPart)*(int64_t(1<<16)*int64_t(1<<16))+c1.LowPart)
				-(int64_t(c0.HighPart)*(int64_t(1<<16)*int64_t(1<<16))+c0.LowPart);
			return (unsigned int)((test/(f.LowPart/1000))%(int64_t(1<<16)*int64_t(1<<16)));
		}
		else
			return 0;
}
#endif
