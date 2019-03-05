/***************************************************************************
 *
 *   string.cpp -- Japanized subset of MFC strcore.cpp
 *
 *   Erik Persson and Björn Nilsson, 2000-2008
 */

#ifndef MFC_CLASSES

#include <stdio.h>
#include "types/string.h"
#include "types/minmax.h"

#define TRACE0(x) (1)

#define InterlockedIncrement(plong) (++(*(plong)))
#define InterlockedDecrement(plong) (--(*(plong)))

typedef char BYTE;
typedef CString::CHAR TCHAR;

/***************************************************************************/

// Global data

// g_chNil is left for backward compatibility
TCHAR g_chNil = '\0';

// For an empty string, m_pElements will point here
// (note: avoids special case of checking for NULL m_pElements)
// empty string data (and locked)

// int(-1), size_t(0), size_t(0), followed by at least one char(0)
// works on both 32-bit and 64-bit systems
static int g_InitData[] = { -1, 0, 0, 0, 0, 0, 0 }; 

static CStringData* g_dataNil = (CStringData*)&g_InitData;

// Pointer into the buffer part of g_dataNil
const TCHAR *g_pchNil = (const TCHAR *)(((BYTE*)&g_InitData)+sizeof(CStringData));

/***************************************************************************/

CString::CString(const CString& stringSrc)
{
	ASSERT(stringSrc.GetData()->nRefs != 0);
	if (stringSrc.GetData()->nRefs >= 0)
	{
		ASSERT(stringSrc.GetData() != g_dataNil);
		m_pElements = stringSrc.m_pElements;
		InterlockedIncrement(&GetData()->nRefs);
	}
	else
	{
		Init();
		*this = stringSrc.m_pElements;
	}
}

CString::CString(TCHAR ch, size_t nLength)
{
	Init();
	if (nLength >= 1)
	{
		AllocBuffer(nLength);
		memset(m_pElements, ch, nLength*sizeof(TCHAR));
	}
}

CString::CString(const TCHAR *pstr, size_t nLength)
{
	Init();
	if (nLength != 0)
	{
		//ASSERT(AfxIsValidAddress(pstr, nLength, FALSE));
		AllocBuffer(nLength);
		memcpy(m_pElements, pstr, nLength*sizeof(TCHAR));
	}
}

CString::CString(const TCHAR *pstr)
{
	Init();

  if (pstr) // MFC Compatible 0-tolerance
  {
    // Removed Implicit LoadString removed /Erik
    // ASSERT(((unsigned int int) pstr) >= 0x10000);
	  size_t nLen = strlen(pstr);
	  if (nLen != 0)
	  {
		  AllocBuffer(nLen);
		  memcpy(m_pElements, pstr, nLen*sizeof(TCHAR));
	  }
  }
}

/***************************************************************************/

void CString::AllocBuffer(size_t nLen)
{
  // always allocate one extra character for '\0' termination
  // assumes [optimistically] that data length will equal allocation length

	ASSERT(nLen >= 0);

	if (nLen == 0) 
    Init();
	else
	{
		size_t nAllocSize= sizeof(CStringData) + (nLen+1)*sizeof(TCHAR);
		CStringData* pData= (CStringData*) new char[nAllocSize];
		pData->nAllocLength = nLen;
		pData->nRefs = 1;

		TCHAR *aElements= (TCHAR *) (pData+1);
		aElements[nLen] = '\0';
		pData->nDataLength = nLen;
		m_pElements= aElements;
	}
}

void CString::FreeData(CStringData* pData)
{
	delete[] (BYTE*)pData;
}

void CString::Release()
{
	if (GetData() != g_dataNil)
	{
		ASSERT(GetData()->nRefs != 0);
		if (InterlockedDecrement(&GetData()->nRefs) <= 0)
			FreeData(GetData());
		Init();
	}
}

void CString::Release(CStringData* pData)
{
	if (pData != g_dataNil)
	{
		ASSERT(pData->nRefs != 0);
		if (InterlockedDecrement(&pData->nRefs) <= 0)
			FreeData(pData);
	}
}

void CString::Empty()
{
	if (GetData()->nDataLength == 0)
		return;
	if (GetData()->nRefs >= 0)
		Release();
	else
		*this = &g_chNil;
	ASSERT(GetData()->nDataLength == 0);
	ASSERT(GetData()->nRefs < 0 || GetData()->nAllocLength == 0);
}

void CString::CopyBeforeWrite()
{
	if (GetData()->nRefs > 1)
	{
		CStringData* pData = GetData();
		Release();
		AllocBuffer(pData->nDataLength);
		memcpy(m_pElements, pData+1, (pData->nDataLength+1)*sizeof(TCHAR));
	}
	ASSERT(GetData()->nRefs <= 1);
}

void CString::AllocBeforeWrite(size_t nLen)
{
	CStringData *pData= GetData();
	if (pData->nRefs > 1 || nLen > pData->nAllocLength)
	{
		Release();
		AllocBuffer(nLen);
	}
	ASSERT(GetData()->nRefs <= 1);
}

CString::~CString()
//  free any attached data
{
	if (GetData() != g_dataNil)
	{
		if (InterlockedDecrement(&GetData()->nRefs) <= 0)
			FreeData(GetData());
	}
}

/***************************************************************************/

// Assignment operators
//  All assign a new value to the string
//      (a) first see if the buffer is big enough
//      (b) if enough room, copy on top of old buffer, set size and type
//      (c) otherwise free old string data, and create a new one
//
//  All routines return the new string (but as a 'const CString&' so that
//      assigning it again will cause a copy, eg: s1 = s2 = "hi there".
//

void CString::AssignCopy(size_t nSrcLen, const TCHAR *aSrcData)
{
	AllocBeforeWrite(nSrcLen);
	memcpy(m_pElements, aSrcData, nSrcLen*sizeof(TCHAR));
	GetData()->nDataLength = nSrcLen;
	m_pElements[nSrcLen] = '\0';
}

const CString& CString::operator=(const CString& stringSrc)
{
	if (m_pElements != stringSrc.m_pElements)
	{
		if ((GetData()->nRefs < 0 && GetData() != g_dataNil) ||
			stringSrc.GetData()->nRefs < 0)
		{
			// actual copy necessary since one of the strings is locked
			AssignCopy(stringSrc.GetData()->nDataLength, stringSrc.m_pElements);
		}
		else
		{
			// can just copy references around
			Release();
			ASSERT(stringSrc.GetData() != g_dataNil);
			m_pElements = stringSrc.m_pElements;
			InterlockedIncrement(&GetData()->nRefs);
		}
	}
	return *this;
}

const CString& CString::operator=(const TCHAR *pstr)
{
	// ASSERT(pstr == NULL || AfxIsValidString(pstr));
	size_t nLen= pstr ? strlen(pstr) : 0; // MFC Compatible 0-tolerance
	AssignCopy(nLen, pstr);
	return *this;
}

/***************************************************************************/

// Concatenation

// NOTE: "operator+" is done as friend functions for simplicity
//      There are three variants:
//          CString + CString
// and for ? = TCHAR, const TCHAR *
//          CString + ?
//          ? + CString

void CString::ConcatCopy(size_t nSrc1Len, const TCHAR *aSrc1Data, size_t nSrc2Len, const TCHAR *aSrc2Data)
{
  // -- master concatenation routine
  // Concatenate two sources
  // -- assume that 'this' is a new CString object

	size_t nNewLen = nSrc1Len + nSrc2Len;
	if (nNewLen != 0)
	{
		AllocBuffer(nNewLen);
		memcpy(m_pElements, aSrc1Data, nSrc1Len*sizeof(TCHAR));
		memcpy(m_pElements+nSrc1Len, aSrc2Data, nSrc2Len*sizeof(TCHAR));
	}
}

CString operator+(const CString& string1, const CString& string2)
{
	CString s;
	s.ConcatCopy(string1.GetData()->nDataLength, string1.m_pElements,
		string2.GetData()->nDataLength, string2.m_pElements);
	return s;
}

CString operator+(const CString& string, const TCHAR *pstr)
{
	CString s;
	if (pstr) // MFC Compatible 0-tolerance
	{
		// ASSERT(AfxIsValidString(pstr));
		s.ConcatCopy(string.GetData()->nDataLength, string.m_pElements, strlen(pstr), pstr);
	}
	return s;
}

CString operator+(const TCHAR *pstr, const CString& string)
{
	// ASSERT(pstr == NULL || AfxIsValidString(pstr));
	CString s;
	if (pstr) // MFC Compatible 0-tolerance
		s.ConcatCopy(strlen(pstr), pstr, string.GetData()->nDataLength, string.m_pElements);

	return s;
}

/***************************************************************************/

// Concatenate in-place

void CString::ConcatInPlace(size_t nSrcLen, const TCHAR *aSrcData)
{
	//  -- the main routine for += operators

	// concatenating an empty string is a no-op!
	if (nSrcLen == 0)
		return;

	// if the buffer is too small, or we have a width mis-match, just
	//   allocate a new buffer (slow but sure)
	if (GetData()->nRefs > 1 || GetData()->nDataLength + nSrcLen > GetData()->nAllocLength)
	{
		// we have to grow the buffer, use the ConcatCopy routine
		CStringData* pOldData = GetData();
		ConcatCopy(GetData()->nDataLength, m_pElements, nSrcLen, aSrcData);
		ASSERT(pOldData != NULL);
		CString::Release(pOldData);
	}
	else
	{
		// fast concatenation when buffer big enough
		memcpy(m_pElements+GetData()->nDataLength, aSrcData, nSrcLen*sizeof(TCHAR));
		GetData()->nDataLength += nSrcLen;
		ASSERT(GetData()->nDataLength <= GetData()->nAllocLength);
		m_pElements[GetData()->nDataLength] = '\0';
	}
}

void CString::Append(const CString &sSrc)
{
	ConcatInPlace(sSrc.GetLength(), sSrc);
}

const CString& CString::operator+=(const TCHAR *pstr)
{
	if (pstr) // MFC Compatible 0-tolerance
	{
		// ASSERT(AfxIsValidString(pstr));
		ConcatInPlace(strlen(pstr), pstr);
	}
	return *this;
}

const CString& CString::operator+=(TCHAR ch)
{
	ConcatInPlace(1, &ch);
	return *this;
}

const CString& CString::operator+=(const CString& string)
{
	ConcatInPlace(string.GetData()->nDataLength, string.m_pElements);
	return *this;
}

/***************************************************************************/
// Advanced direct buffer access

TCHAR *CString::GetBuffer(size_t nMinBufLength)
{
	ASSERT(nMinBufLength >= 0);

	if (GetData()->nRefs > 1 || nMinBufLength > GetData()->nAllocLength)
	{
#ifdef _DEBUG
		// give a warning in case locked string becomes unlocked
		if (GetData() != g_dataNil && GetData()->nRefs < 0)
			TRACE0("Warning: GetBuffer on locked CString creates unlocked CString!\n");
#endif
		// we have to grow the buffer
		CStringData* pOldData = GetData();
		size_t nOldLen = GetData()->nDataLength;   // AllocBuffer will tromp it
		if (nMinBufLength < nOldLen)
			nMinBufLength = nOldLen;
		AllocBuffer(nMinBufLength);
		memcpy(m_pElements, pOldData+1, (nOldLen+1)*sizeof(TCHAR));
		GetData()->nDataLength = nOldLen;
		CString::Release(pOldData);
	}
	ASSERT(GetData()->nRefs <= 1);

	// return a pointer to the character storage for this string
	ASSERT(m_pElements != NULL);
	return m_pElements;
}

/*
// removed because size_t can be either unsigned or signed, depending on the target platform, meaning -1 can crash the code
void CString::ReleaseBuffer(size_t nNewLength)
{
	CopyBeforeWrite();  // just in case GetBuffer was not called

	if (nNewLength == -1)
		nNewLength = strlen(m_pElements); // zero terminated

	ASSERT(nNewLength <= GetData()->nAllocLength);
	GetData()->nDataLength = nNewLength;
	m_pElements[nNewLength] = '\0';
}
*/

TCHAR *CString::GetBufferSetLength(size_t nNewLength)
{
	ASSERT(nNewLength >= 0);

	GetBuffer(nNewLength);
	GetData()->nDataLength = nNewLength;
	m_pElements[nNewLength] = '\0';
	return m_pElements;
}

void CString::FreeExtra()
{
	ASSERT(GetData()->nDataLength <= GetData()->nAllocLength);
	if (GetData()->nDataLength != GetData()->nAllocLength)
	{
		CStringData* pOldData = GetData();
		AllocBuffer(GetData()->nDataLength);
		memcpy(m_pElements, pOldData+1, pOldData->nDataLength*sizeof(TCHAR));
		ASSERT(m_pElements[GetData()->nDataLength] == '\0');
		CString::Release(pOldData);
	}
	ASSERT(GetData() != NULL);
}

TCHAR *CString::LockBuffer()
{
	TCHAR *pstr = GetBuffer(0);
	GetData()->nRefs = -1;
	return pstr;
}

void CString::UnlockBuffer()
{
	ASSERT(GetData()->nRefs == -1);
	if (GetData() != g_dataNil)
		GetData()->nRefs = 1;
}

void CString::SetAt(size_t nIndex, TCHAR ch)
{
	ASSERT(nIndex >= 0);
	ASSERT(nIndex < GetData()->nDataLength);

	CopyBeforeWrite();
	m_pElements[nIndex] = ch;
}

/*
int GetFormattedLength(const TCHAR *pszFormat, va_list args) 
{
	return _vscprintf( pszFormat, args );
}

int CString::Format(const TCHAR *pszFormat, ...)
{
	// Returns number of chars written, not including terminating null character
	if (!pszFormat)
	{
		Empty();
		return 0;
	}

	va_list args;
	va_start(args,pszFormat);

	int nLength= GetFormattedLength( pszFormat, args );
	TCHAR *pszBuffer= GetBuffer(nLength);
	int nWritten= vsprintf(pszBuffer, pszFormat, args);
	ReleaseBuffer(nLength);

	va_end(args);
	return nWritten;
}
*/

// Added by bn 2008-05-30
CString CString::Left(size_t nCount) const
{
	if (IsEmpty())
		return "";
	if (nCount<0)
		nCount= 0;
	size_t n= GetLength();
	if (nCount>n)
		nCount= n;
	return CString(m_pElements, nCount);
}

// Added by bn 2008-05-30
CString CString::Right(size_t nCount) const
{
	if (IsEmpty())
		return "";
	if (nCount<0)
		nCount= 0;
	size_t n= GetLength();
	if (nCount>n)
		nCount= n;
	return CString(m_pElements+n-nCount, nCount);
}

// Added by bn 2008-05-30
CString CString::Mid(size_t nFirst, size_t nCount) const
{
	if (IsEmpty())
		return "";
	size_t n= GetLength();
	if (nFirst>=n)
		return "";
	if (nFirst<0)
	{
		nCount= nCount-nFirst;
		nFirst= 0;
	}
	if (nCount<0)
		nCount= 0;
	if (nFirst+nCount>n)
		nCount= n-nFirst;
	return CString(m_pElements+nFirst, nCount);
}

int64_t CString::Find(const TCHAR *aSub, int64_t nSub) const
{
	if (aSub==NULL || *aSub==0)
		return -1;
	if (nSub==0)
		return -1;

	const int64_t n_this= GetLength();
	if (n_this<nSub)
		return -1;

	for (int64_t i=0;i<=n_this-nSub;)
	{
		int64_t j;
		for (j=nSub;--j>=0;)
			if (aSub[j]!=m_pElements[i+j])
				break;
		if (j<0)
			return i;
		i++; // TODO: This is brute force. Improve!
	}
	return -1;
}

int64_t CString::Find(const TCHAR c) const
{
	for (size_t i=0;i<GetLength();i++)
		if (m_pElements[i]==c)
			return i;
	return -1;
}

int64_t CString::Find(const CString &sSub) const
{
	return Find(sSub, sSub.GetLength());
}

#endif // MFC_CLASSES
