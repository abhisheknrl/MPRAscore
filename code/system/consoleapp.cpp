/*******************************************************************************
 *
 *   consoleapp.cpp -- 
 *
 *   Björn Nilsson, 2004-2008
 */

#include "consoleapp.h"
#include <stdio.h>
#include "filehelpers.h"
#include "types/stringfun.h"
#include "filesystem.h"
#include <stdlib.h>

const char *CConsoleApp::ScanSwitch(const char *ach, CString &sSwitch)
{ 
	// Supports - or -- as (synonymous) switch prefixes
	if (ach[0]=='-')
	{
		ach++;
		if (ach[0]=='-')
			ach++;
		const char *ach0= ach;
		while (*ach && ach[0]!=':')
			ach++;
		sSwitch= CString(ach0, ach-ach0);
		if (ach[0]==':')
			ach++;
		return ach;
	}

	return NULL;
}

bool CConsoleApp::ParseCommandLine(int argc, TCHAR* argv[], TCHAR* envp[])
{
	g_FileSystem.SetExecutableFileName(argv[0]);
	m_vNonSwitchArgs.SetSize(0);

	bool bOk= true;
	for (int i=1;bOk && i<argc;i++)
	{
		CString sSwitch;
		const char *arg= argv[i];
		if ((arg= ScanSwitch(arg, sSwitch)))
		{
			if (sSwitch=="verbosity")
				SetVerbosity(atoi(arg));
			else
				bOk= OnSwitch((const char *)arg, sSwitch);
		}
		else 
		{
			CString s= argv[i];
			m_vNonSwitchArgs.Add(g_FileSystem.ToSystemSlash(s));
		}
	}

	return bOk;
}
