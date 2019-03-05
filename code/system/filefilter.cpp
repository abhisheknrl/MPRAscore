/*******************************************************************************
 *
 *   filefilter.cpp -- 
 *
 *   Björn Nilsson, 2002
 */

#include "system/filefilter.h"

void CFileFilter::FreeBuffer()
{
	if (m_paBuffer)
		delete [] m_paBuffer;

	m_paBuffer = NULL;
}

void CFileFilter::AddFilter(CString sTitle, CString sFilter)
{
	m_vsTitles.Add(sTitle + " ("+sFilter+")");
	m_vsFilters.Add(sFilter);
}

const char *CFileFilter::GetFilterArray() // The returned pointer must not be tampered with.
{
	FreeBuffer();

	if (m_vsTitles.GetSize() != m_vsFilters.GetSize())
	{
		ASSERT(false);
		return NULL;
	}

	const int N = m_vsTitles.GetSize();

	int nSize = 2;
	for (int i=0;i<N;i++)
		nSize += m_vsTitles[i].GetLength() + 1;
	for (int i=0;i<N;i++)
		nSize += m_vsFilters[i].GetLength() + 1;

	try
	{
		m_paBuffer = new char[nSize];
	}
	catch(...)
	{
		return NULL;
	}

	// Copy strings.
	char *ach = m_paBuffer;

	for (int i=0;i<N;i++)
	{
		CString s;
		int j;

		s = m_vsTitles[i];
		for (j=0;j<s.GetLength();j++)
			*ach++ = s[j];
		*ach++ = 0;

		s = m_vsFilters[i];
		for (j=0;j<s.GetLength();j++)
			*ach++ = s[j];
		*ach++ = 0;
	}

	*ach++ = 0;
	*ach++ = 0;

	return m_paBuffer;
}
