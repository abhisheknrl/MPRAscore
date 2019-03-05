/***********************************************************************************
 *
 *   CStringRef.h - Reference to a string within a file that is being parsed
 *
 *   Erik Persson and Björn Nilsson, 1999-2008
 */

#ifndef STRINGREF_H
#define STRINGREF_H

#include "types/string.h"

class CStringRef 
{
public:
	const char *m_ach;
	size_t m_Len;
	
	CStringRef() { m_ach = NULL; m_Len=0; } // Constructor introduced because of dangling references in Release mode (uninitialized memory!).
	CStringRef(const char *ach) { m_ach= ach; m_Len= strlen(ach); }
	CStringRef(const char *ach, size_t Len) { m_ach= ach; m_Len= Len; }
	CStringRef(CString s) { m_ach= /*(LPCTSTR)*/ s; m_Len= s.GetLength(); }

	bool operator==(char *aStr) const {	return strncmp(m_ach, aStr, m_Len)==0 && aStr[m_Len]==0; }
	char operator[](size_t nIndex) const {	return m_ach[nIndex]; }
	operator CString() const { return CString(m_ach, m_Len); } 

	void Assign(char *ach, size_t Len) { m_ach= ach; m_Len= Len; }
	void Empty() { m_Len= 0; }
	size_t GetLength() { return m_Len; }
	size_t length() { return m_Len; }	
};

#endif // STRINGREF_H
