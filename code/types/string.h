/***************************************************************************
 *
 *   string.h -- Japanized MFC Style Subset String Class
 *
 *   Erik Persson and Björn Nilsson, 2000-
 */

#ifndef TYPES_STRING_H
#define TYPES_STRING_H

#ifdef MFC_CLASSES
#include <afx.h>
#else

#include <stdarg.h>
#include <string.h>
#include "types/assert.h"
#include <cstdint>

typedef char TCHAR;
typedef unsigned int UINT;

#define g_stringNil ((CString&)*(CString*)&g_pchNil)

struct CStringData
{
	int nRefs;                // reference count
	size_t nDataLength;        // length of data (including terminator)
	size_t nAllocLength;       // length of allocation

  // Elements are placed after this structure
};

extern const TCHAR *g_pchNil;

class CString
{
public:
#ifdef _WIN32
	typedef TCHAR CHAR;
#else
	typedef char CHAR;
#endif

// Constructors
	CString();                                // constructs empty CString
	CString(const CString& stringSrc);        // copy constructor
	CString(CHAR ch, size_t nRepeat = 1);        // from a single character	
	CString(const CHAR *pstr);                // from a string
	CString(const CHAR *pstr, size_t nLength);   // subset of characters from string

// Attributes & Operations
	size_t GetLength() const; 	                  // get data length
	bool IsEmpty() const;                     // TRUE if zero length
	bool IsSpaces() const;
	void Empty();                             // clear contents to empty	
	CHAR GetAt(size_t nIndex) const;             // return single character at zero-based index
	// CHAR operator[](int nIndex) const;        // return single character at zero-based index
	void SetAt(size_t nIndex, CHAR ch);          // set a single character at zero-based index
	operator const CHAR *() const;            // return pointer to const string

	// overloaded assignment	
	const CString& operator=(const CString& stringSrc); // ref-counted copy from another CString
	const CString& operator=(CHAR ch);                  // set string content to single character
	const CString& operator=(const CHAR *pstr);         // copy string content from ANSI string (converts to CHAR)
	
	// string concatenation	
	const CString& operator+=(const CString& string); // concatenate from another CString	
	const CString& operator+=(CHAR ch);               // concatenate a single character
	const CString& operator+=(const CHAR *pstr);      // concatenate a string

	friend CString operator+(const CString& string1, const CString& string2);
	// friend CString operator+(const CString& string, CHAR ch);
	// friend CString operator+(CHAR ch, const CString& string);
	friend CString operator+(const CString& string, const CHAR *pstr);
	friend CString operator+(const CHAR *pstr, const CString& string);

	// substring extraction
	CString Left(size_t nCount) const;
	CString Mid(size_t nFirst, size_t nCount) const;
	CString Right(size_t nCount) const;
	void Append(const CString &sSrc);

	// string comparison	
	int Compare(const CHAR *pstr) const;             // straight character comparison

	// Access to string implementation buffer as "C" character array	
	CHAR *GetBuffer(size_t nMinBufLength);        // get pointer to modifiable buffer at least as long as nMinBufLength
	// void ReleaseBuffer(size_t nNewLength = -1);	 // release buffer, setting length to nNewLength (or to first nul if -1)
	CHAR *GetBufferSetLength(size_t nNewLength);  // get pointer to modifiable buffer exactly as long as nNewLength
	void FreeExtra();                          // release memory allocated to but unused by string

	// Use LockBuffer/UnlockBuffer to turn refcounting off
	CHAR *LockBuffer();  // turn refcounting back on	
	void UnlockBuffer(); // turn refcounting off

	// Substring searching
	// These functions return ints because size_t is unsigned on _W64
	int64_t Find(const CString &sSub) const;
	int64_t Find(const TCHAR c) const;
	int64_t Find(const TCHAR *aSub, int64_t n) const;

// Implementation
public:
	~CString();
	size_t GetAllocLength() const;

protected:
	CHAR *m_pElements;   // pointer to ref counted string data

	// implementation helpers
public:
	CStringData* GetData() const;
protected:
	void Init();
	void AllocBuffer(size_t nLen);
	void AssignCopy(size_t nSrcLen, const CHAR *aSrcData);
	void ConcatCopy(size_t nSrc1Len, const CHAR *aSrc1Data, size_t nSrc2Len, const CHAR *aSrc2Data);
	void ConcatInPlace(size_t nSrcLen, const CHAR *aSrcData);
	void CopyBeforeWrite();
	void AllocBeforeWrite(size_t nLen);
	void Release();
	static void Release(CStringData* pData);
	static void FreeData(CStringData* pData);

};

/*********** Inline definitions ***********/

// Compare helpers
inline bool operator==(const CString& s1, const CString& s2)        { return s1.Compare(s2) == 0; }
inline bool operator==(const CString& s1, const CString::CHAR *a2)  { return s1.Compare(a2) == 0; }
inline bool operator==(const CString::CHAR *a1, const CString& s2)  { return s2.Compare(a1) == 0; }
inline bool operator!=(const CString& s1, const CString& s2)        { return s1.Compare(s2) != 0; }
inline bool operator!=(const CString& s1, const CString::CHAR *a2)  { return s1.Compare(a2) != 0; }
inline bool operator!=(const CString::CHAR *a1, const CString& s2)  { return s2.Compare(a1) != 0; }
inline bool operator<(const CString& s1, const CString& s2)         { return s1.Compare(s2) < 0; }
inline bool operator<(const CString& s1, const CString::CHAR *a2)   { return s1.Compare(a2) < 0; }
inline bool operator<(const CString::CHAR *a1, const CString& s2)   { return s2.Compare(a1) > 0; }
inline bool operator>(const CString& s1, const CString& s2)         { return s1.Compare(s2) > 0; }
inline bool operator>(const CString& s1, const CString::CHAR *a2)   { return s1.Compare(a2) > 0; }
inline bool operator>(const CString::CHAR *a1, const CString& s2)   { return s2.Compare(a1) < 0; }
inline bool operator<=(const CString& s1, const CString& s2)        { return s1.Compare(s2) <= 0; }
inline bool operator<=(const CString& s1, const CString::CHAR *a2)  { return s1.Compare(a2) <= 0; }
inline bool operator<=(const CString::CHAR *a1, const CString& s2)  { return s2.Compare(a1) >= 0; }
inline bool operator>=(const CString& s1, const CString& s2)        { return s1.Compare(s2) >= 0; }
inline bool operator>=(const CString& s1, const CString::CHAR *a2)  { return s1.Compare(a2) >= 0; }
inline bool operator>=(const CString::CHAR *a1, const CString& s2)  { return s2.Compare(a1) <= 0; }

// CString
inline CStringData* CString::GetData() const            { ASSERT(m_pElements); return ((CStringData*)m_pElements)-1; }
inline void CString::Init()                             { m_pElements= (CHAR *) g_pchNil; }
inline CString::CString()                               { m_pElements= (CHAR *) g_pchNil; }

inline size_t CString::GetLength() const                   { return GetData()->nDataLength; }
inline size_t CString::GetAllocLength() const              { return GetData()->nAllocLength; }
inline bool CString::IsEmpty() const                    { return GetData()->nDataLength == 0; }
inline bool CString::IsSpaces() const  { 
	for (int i=0;i<GetLength();i++)
		if (m_pElements[i]!=32 && (unsigned char)(m_pElements[i])!=160 && m_pElements[i]!=9)
			return false;
	return true; // if all chars are spaces or tabs, or string is empty
}
inline CString::operator const CString::CHAR *() const  { return m_pElements; }

inline int CString::Compare(const CHAR *pstr) const     { return strcmp(m_pElements, pstr); }

inline CString::CHAR CString::GetAt(size_t nIndex) const
{	ASSERT(nIndex >= 0 &&	nIndex < GetData()->nDataLength);
	return m_pElements[nIndex];
}

/*
inline CString::CHAR CString::operator[](int nIndex) const
{	// Same as GetAt
	ASSERT(nIndex >= 0);
	ASSERT(size_t(nIndex) < GetData()->nDataLength);
	return m_pElements[nIndex];
}
*/

#endif

#endif // TYPES_STRING_H
