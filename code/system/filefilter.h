/*******************************************************************************
 *
 *   filefilter.h -- 
 *
 *   Björn Nilsson, Sketchware, 2002
 */

#ifndef FILEFILTER_H
#define FILEFILTER_H

#include "afxwin.h"
#include "types/vector.h"

class CFileFilter
{
protected:
	char *m_paBuffer;
	CVector<CString> m_vsTitles;
	CVector<CString> m_vsFilters;

	void FreeBuffer();
public:
	const char *GetFilterArray(); // The returned pointer must not be tampered with.
	void AddFilter(CString sTitle, CString sFilter);

	CFileFilter()	{	m_paBuffer = NULL;	}
	~CFileFilter()	{	FreeBuffer();	}
};

#endif
