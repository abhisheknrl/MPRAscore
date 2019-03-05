//
// filehelpers.cpp --
//

#include "filehelpers.h"
#include <time.h>
#include <string.h> // for strlen
#include "stdout.h"
#include "stderr.h"

int g_Verbosity= 2;

bool CheckLicence(const int nDays)
{
	time_t ltime, ref_ltime= 1365863093; // Apr 13, 2013
	time( &ltime );
	ltime -= ref_ltime;

	if (ltime > 0 && ltime < nDays*24*60*60)
		return true;
	return false;
}

void ReportError(const char *aMsg, const char *aContext, int nVerbosity)
{
	if (nVerbosity<=0)
		return;
	if (!aMsg)
		return;
	eprintf("Error: %s", aMsg);
	if (aContext && aContext[0])
	{
		eprintf(": %s\n", aContext);
		return;
	}
	
	// No context
	size_t n= strlen(aMsg);
	if (n>0 && aMsg[n-1]!='.' && aMsg[n-1]!='?' && aMsg[n-1]!='!')
		eprintf(".");
	eprintf("\n");
}

void ReportError(const char *aMsg, const char *aContext) 
{ 
	ReportError(aMsg, aContext, g_Verbosity); 
}

void ReportWarning(const char *aMsg, const char *aContext, int nVerbosity)
{
	if (nVerbosity<=1)
		return;
	if (!aMsg)
		return;
	printf("Warning: %s", aMsg);
	if (aContext && aContext[0])
	{
		printf(": %s\n", aContext);
		return;
	}
	
	// No context
	size_t n= strlen(aMsg);
	if (n>0 && aMsg[n-1]!='.' && aMsg[n-1]!='?' && aMsg[n-1]!='!')
		printf(".");
	printf("\n");
}

void ReportWarning(const char *aMsg, const char *aContext)
{ 
	ReportWarning(aMsg, aContext, g_Verbosity);
}

FILE *OpenInputFile(const char *aFilename, int nVerbosity)
{
	FILE *pf= NULL;
	if (aFilename && aFilename[0]!=0)
		pf= fopen(aFilename, "rb");
	if (!pf)
		ReportError("Cannot open input file", aFilename, nVerbosity);
	return pf;
}

FILE *OpenInputFile(const char *aFilename) 
{ 
	return OpenInputFile(aFilename, g_Verbosity); 
}

FILE *OpenOutputFile(const char *aFilename, int nVerbosity)
{
	FILE *pf= NULL;
	if (aFilename && aFilename[0]!=0)
		pf= fopen(aFilename, "wb");
	if (pf)
	{
		if (nVerbosity>=2)
			printf("Creating output file: %s\n", aFilename);
	}
	else
		ReportError("Cannot open output file", aFilename, nVerbosity);
	return pf;
}

FILE *OpenOutputFile(const char *aFilename, int nVerbosity, const char *aMessage)
{
	FILE *pf= NULL;
	if (aFilename && aFilename[0]!=0)
		pf= fopen(aFilename, "wb");
	if (pf)
	{
		if (nVerbosity>=2 && aMessage)
			printf("%s: %s\n", aMessage, aFilename);
	}
	else
		ReportError("Cannot open output file", aFilename, nVerbosity);
	return pf;
}

FILE *OpenOutputFile(const char *aFilename) 
{ 
	return OpenOutputFile(aFilename, g_Verbosity); 
}

bool CloseFile(FILE *pf)
{
	if (!pf)
		return false;
	return fclose(pf)==0;
}

int GetVerbosity()
{
	return g_Verbosity;
}

int SetVerbosity(int v)
{
	int old= g_Verbosity;
	g_Verbosity= v;
	return old;
}

bool IsLittleEndian()
{
	int x= 1;
	unsigned char y= *(const char *)&x;
	return y==1;
}

bool IsBigEndian()
{
	int x= 1;
	unsigned char y= *(const char *)&x;
	return y==0;
}
