/*******************************************************************************
 *
 *   annotation.cpp -- 
 *
 *   Björn Nilsson, Veta Mera Data HB, 2005
 */

#include "annotation.h"

// CColorMap

CColorMap::CColorMap(const CColorMap& other)
{
	m_vsKey= other.m_vsKey;
	m_vcValue= other.m_vcValue;
}

int CColorMap::FindKey(const char *aProp) const
{
	// Return index of property or -1 for failure
	for (int i= 0;i<m_vcValue.GetSize();i++)
		if (!strcmp(m_vsKey[i],aProp))
			return i;
	return -1;
}

bool CColorMap::AppendAssoc(const char *aProp, COLORREF aCol)
{
	try
	{
		m_vsKey.Add(aProp);
		m_vcValue.Add(aCol);
		return true;
	}
	catch(...) // CMemoryException
	{
		return false;
	}
}

COLORREF CColorMap::GetColor(const char *aProp) const
{
	int n= FindKey(aProp);
	if (n>=0)
	{
		ASSERT( n<m_vcValue.GetSize() );
		return m_vcValue[n];
	}
	return COLORREF(0x808080); // Return neutral color
}

bool CColorMap::SetColor(const char *aProp, COLORREF aCol)
{
	try
	{
		int n= FindKey(aProp);
		if (n>=0)
		{
			ASSERT( n<m_vcValue.GetSize() );
			m_vcValue.SetAt(n, aCol);
			return true;
		}
		else
			return AppendAssoc(aProp, aCol);
	}
	catch(...) // CMemoryException
	{
		return false;
	}
}

int CColorMap::GetCount() const
{
	return m_vcValue.GetSize();
}

CString CColorMap::GetKey(int n) const
{
	ASSERT(n>=0 && n<m_vsKey.GetSize());
	return m_vsKey[n];
}

COLORREF CColorMap::GetValue(int n) const
{
	ASSERT(n>=0 && n<m_vcValue.GetSize());
	return m_vcValue[n];
}

