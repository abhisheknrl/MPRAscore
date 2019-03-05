/***********************************************************************************
 *
 *   stdout.cpp -- Standard Output using ITextOutput
 *
 */

#include <stdio.h>
#include "stderr.h"
#include "textoutput.h"
#include "string.h"

/***********************************************************************************/
// CMyNilOutput: default output object

class CMyNilErr: public ITextOutput
{
	void OnPut(const char *ach, size_t len) { printf("%s", (const char *)CString(ach, len)); }
};

static CMyNilErr g_NilErr;

/***********************************************************************************/

static ITextOutput *g_pStderr= &g_NilErr;

ITextOutput *SetStderr(ITextOutput *pOut)
{
	// Set the standard error device and return the previous setting
	ASSERT(pOut);
	ITextOutput *pOld= g_pStderr;
	g_pStderr= pOut;
	return pOld;
}

void eprintf(const char *ach, ...)
{
	va_list args;
	va_start(args,ach);
	g_pStderr->WriteFmtVA(ach,args);
	va_end(args);
}

void eputchar(char c)
{
	g_pStderr->Write(c);
}

void eputs(const char *ach)
{
	// puts, differently from fputs, automatically appends newline
	g_pStderr->Write(ach);
	g_pStderr->Write('\n');
}
