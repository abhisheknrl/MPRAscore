/***********************************************************************************
 *   
 *	 stringfun.cpp -- String utility functions
 *
 *   Erik Persson and Björn Nilsson, 1998-
 *
 */

#include <stdio.h>
#include <locale> // For toupper and tolower
// #include <stdlib.h>
// #include <stdarg.h>
#include "string.h"
#include "stringfun.h"

#define FORMAT_BUFFER_SIZE (0x8000)

CString Format(const char *ach, ...)
{
	// NOTE: vsprintf was changed to _vsnprintf because of the classical 
	//		 security breach: no check against buffer size.
	va_list args;
	va_start(args,ach);
	char aBuf[FORMAT_BUFFER_SIZE ];
	vsnprintf(aBuf, FORMAT_BUFFER_SIZE , ach, args);
	va_end(args);
	return CString(aBuf);
}

CString FormatVA(const char *ach, va_list argptr)
{
	// NOTE: vsprintf was changed to _vsnprintf because of the classical 
	//		 security breach: no check against buffer size.
	char aBuf[FORMAT_BUFFER_SIZE];
	vsnprintf(aBuf, FORMAT_BUFFER_SIZE, ach, argptr);
	return CString(aBuf);
}

bool ParseInt(const char *ach, int *pResult)
{
	// TODO: Write a new version of this that checks for garbage at end
	return sscanf(ach,"%d",pResult)==1;
}

const char *strreverse(char *ach)
{
	// Reverses the order of the characters in a string
	// Returns the pointer to the original string (redundant)
	size_t n= strlen(ach) >> 1;
	for (size_t i=0;i<n;i++)
		ach[i]= ach[n-1-i];
	return ach;
}

int strcmp_nocase(const char *ach, const char *ach2)
{
	/* Case-insensitive strcmp */
	while (*ach)
	{
		int c, c2;
		c= tolower(*ach);
		c2= tolower(*ach2);

		if (c>c2) return 1;
		if (c<c2) return -1;
		ach++;
		ach2++;
	}

	if (*ach2) 
		return -1; /* Had a bug: returned 1 */
	else 
		return 0;
}

int strfind(const char *a0, const CVector<CString> &v)
{
	for (int i=0;i<v.GetSize();i++)
		if (strcmp(a0, v[i])==0)
			return i;
	return -1;
}

int strfind(const char *a0, const CVector<const char *> &v)
{
	for (int i=0;i<v.GetSize();i++)
		if (strcmp(a0, v[i])==0)
			return i;
	return -1;
}

int strfind_nocase(const char *a0, const CVector<CString> &v)
{
	for (int i=0;i<v.GetSize();i++)
		if (strcmp_nocase(a0, v[i])==0)
			return i;
	return -1;
}

int strfind_nocase(const char *a0, const CVector<const char *> &v)
{
	for (int i=0;i<v.GetSize();i++)
		if (strcmp_nocase(a0, v[i])==0)
			return i;
	return -1;
}

bool strmatch_recursive(const char *a0, size_t n0, const char *aF, size_t nF)
{
	// a0, n0= string
	// aF, nF= string with * or ? wildcards
	// Also accepts @ as wildcard for letter, and $ as wildcard for a number
	// TODO: Replace with standard regexp syntax
	size_t j= 0;
	for (size_t i=0;i<nF;i++)
	{
		if (aF[i]=='*')
		{
			if (i+1==nF)
				return true;
			for (size_t k=j;k<n0;k++)
				if (strmatch_recursive(&a0[k], n0-k, &aF[i+1], nF-(i+1)))
					return true;
			return false;
		}
		/*
		else if (aF[i]=='[')
		{
			CVector<bool> v_accepted(256);
			v_accepted.Fill(false);

			int i1= i+1;
			int last= 0;
			while (i1<nF && i1!=']')
			{
			}

			if (!v_accepted[a0[j]

			i= i1+1;
		}
		*/
		else if (
			 aF[i]==a0[j] || 
			 aF[i]=='?' || 
			(aF[i]=='@' && ((a0[j]>='a' && a0[j]<='z') || (a0[j]>='A' && a0[j]<='Z'))) ||
			(aF[i]=='$' && a0[j]>='0' && a0[j]<='9'))
		{
			j++;
		}
		else
			return false;
	}
	return j==n0;
}

bool strmatch_recursive_nocase(const char *a0, size_t n0, const char *aF, size_t nF)
{
	// a0, n0= string
	// aF, nF= string with * or ? wildcards
	size_t j= 0;
	for (size_t i=0;i<nF;i++)
	{
		if (aF[i]=='*')
		{
			if (i+1==nF)
				return true;
			for (size_t k=j;k<n0;k++)
				if (strmatch_recursive_nocase(&a0[k], n0-k, &aF[i+1], nF-(i+1)))
					return true;
			return false;
		}
		else if (tolower(aF[i])==tolower(a0[j]) || aF[i]=='?')
			j++;
		else
			return false;
	}
	return j==n0;
}

bool strmatch(const char *a0, const char *aF)
{
	// a0= string
	// aF= string, potentially with * or ? wildcards
	return strmatch_recursive(a0, strlen(a0), aF, strlen(aF));
}

bool strmatch_nocase(const char *a0, const char *aF)
{
	// a0= string
	// aF= string, potentially with * or ? wildcards
	return strmatch_recursive_nocase(a0, strlen(a0), aF, strlen(aF));
}

int strinset(const char *a0, const char *aset, const char csep)
{
	// checks if a given string is included in a string set (multiple strings separated by csep)
	// csep must not be in a0
	// returns index within set, or -1

	size_t n= strlen(a0);
	int p= 0;
	while (*aset)
	{
		size_t i=0;
		for (;i<n;i++)
			if (a0[i]!=aset[i])
				break;
		if (i==n && (aset[i]==0 || aset[i]==csep))
			return p; // found
		while (*aset && *aset!=csep)
			aset++;
		if (*aset)
			aset++;
		p++;
	}

	return -1;
}

bool MatchPrefixNoCase(const char *ach, const char *aPrefix)
{
	/* Case-insensitive prefix match. aPrefix must be lowercase */
	/* Return true if ach matches the prefix */

	char cm;
	while ((cm= *(aPrefix++)))
	{
		if (tolower(*ach) != cm) 
			return false;
		ach++;
	}
	return true;
}

CString ToLower(CString s)
{
	size_t n= s.GetLength();

	bool bDiff= false;
	for (size_t i= 0; i<n; i++)
		if (s[i] != tolower(s[i]))
			bDiff= true;

	if (bDiff)
	{
		// Must modify string
		char *a= s.LockBuffer();
		for (size_t i= 0; i<n; i++)
			a[i]= tolower(a[i]);
		s.UnlockBuffer();
	}
	return s;
}

CString ToUpper(CString s)
{
	size_t n= s.GetLength();

	bool bDiff= false;
	for (size_t i= 0; i<n; i++)
		if (s[i] != toupper(s[i]))
			bDiff= true;

	if (bDiff)
	{
		// Must modify string
		char *a= s.LockBuffer();
		for (size_t i= 0; i<n; i++)
			a[i]= toupper(a[i]);
		s.UnlockBuffer();
	}
	return s;
}

CString BumpName(CString sBaseName)
{
	size_t i= 2;
	size_t baselen= sBaseName.GetLength();

	// Check for decimal number at end
	size_t digs= 0;
	while (baselen && sBaseName[baselen-1] >= '0' && sBaseName[baselen-1] <= '9')
	{
		digs++;
		baselen--;
	}

	if (digs)
	{
		// Convert from decimal and use as starting number
		i= 0;
		for (size_t pos= 0; pos<digs; pos++)
		{
			i *= 10;
			i += sBaseName[baselen+pos]-'0';
		}

		// Increase the number to get a unique name
		i++;

		// Strip off the digits
		sBaseName= sBaseName.Left(baselen);
	}

	// Generate a unique name
	return Format("%s%d", (const char *)sBaseName,i);
}

// NOTE: Some functions rely on GetUniqueName to return a name 
//		 with a digit part that is always higher than or equal to 
//		 the digit part at the end of the base name.
CString GetUniqueName(CString sBaseName, void *pObj, bool (*nameIsUsed)(void *, const char *) )
{
	// Bump name until it gets unique
	while (nameIsUsed(pObj, sBaseName))
		sBaseName= BumpName(sBaseName);

	return sBaseName;
}

CString ToForwardSlash(CString s)
{
	for (size_t i=s.GetLength();--i>=0;)
		if (s[i]=='\\')
			s.SetAt(i, '/');
	return s;
}

CString ToBackSlash(CString s)
{
	for (size_t i=s.GetLength();--i>=0;)
		if (s[i]=='/')
			s.SetAt(i, '\\');
	return s;
}

CString ToLinefeedOnly(CString s)
{
	// Convert from general format EOL ::= LF | CR LF*
	// To Linefeed-only EOL ::= LF
	CString sOut;
	const char *ach= (const char *) s;
	while (char c= *ach++)
	{
		if (c==13)
		{
			while (*ach == 10) // skip LF*
				++ach;
			sOut += char(10);
		}
		else 
			sOut += c; 
	}
	return sOut;
}

CString ToCrlf(CString s)
{
	// Convert from general format EOL ::= LF | CR LF*
	// To CRLF EOL ::= CR LF
	CString sOut;
	const char *ach= (const char *) s;
	while (char c= *ach++)
	{
		if (c==13)
		{
			while (*ach == 10) // skip LF*
				++ach;
			sOut += char(13);
			sOut += char(10);
		}
		else if (c==10)
		{
			sOut += char(13);
			sOut += char(10);
		}
		else
			sOut += c; 
	}
	return sOut;
}


// Added May 2008/BN
void SplitAtChar(const CString &s, const char c, CVector<CString> &vS)
{
	// Identify substrings of s delimited by c. Leading/trailing spaces are removed.

	vS.SetSize(0);

	const char *ach= s;
	if (!(*ach))
		return;
	while (true)
	{
		while (*ach==32) // Discard space prefix
			ach++;

		const char *ach0= ach;
		while (*ach && *ach!=c)
			ach++;

		const char *ach1=ach; // Discard space suffix
		while (ach1>ach0 && *(ach1-1)==32)
			ach1--;

		vS.Add(CString(ach0, ach1-ach0));

		if (!(*ach))
			return;
		ach++;
	}
}

/***********************************************************************************/
// Common parsing functions, moved here from lingo.h/cpp

const char *Scan_UnsignedInt(const char *ach, int *presult) // BN 040123
{
	// Int ::= ("0b" BinDigit BinDigit* | "0x" HexDigit HexDigit* | DecDigit DecDigit*)
	bool bOk= false;
	int d;

	char c= ach[0];
	if (c == '0' && ach[1]=='b' && IsBinDigit(ach[2]))
	{
		// Binary constant
		bOk= true;
		d= GetBinDigit(ach[2]);
		ach += 3;
		while (IsBinDigit(c= ach[0]))
		{
			d= (d<<1) + GetBinDigit(c);
			ach++;
		}
	}
	else if (c == '0' && ach[1]=='x' && IsHexDigit(ach[2]))
	{
		// Binary constant
		bOk= true;
		d= GetHexDigit(ach[2]);
		ach += 3;
		while (IsHexDigit(c= ach[0]))
		{
			d= (d<<4) + GetHexDigit(c);
			ach++;
		}
	}
	else if (IsDecDigit(c))
	{	
		// Dec constant
		bOk= true;
		d= GetDecDigit(c);
		ach++;
		while (IsDecDigit(c= ach[0]))
		{
			d= d*10 + GetDecDigit(c);
			ach++;
		}
	}

	if (bOk)
	{
		*presult= d;
		return ach;
	}

	return 0;
}

const char *Scan_Int(const char *ach, int *presult)
{
	// Int ::= [-] ("0b" BinDigit BinDigit* | "0x" HexDigit HexDigit* | DecDigit DecDigit*)
	bool bNeg= false;
	if (ach[0]=='-') 
	{	
		bNeg= true;
		ach++;
	}

	int d;
	if (!Scan_UnsignedInt(ach, &d))
		return 0;

	if (bNeg) 
		d= -d;
	*presult= d;
	return ach;
}

const char *Scan_Ident(const char *ach)
{
	// Return new position or NULL if failed
	// Lingo uses C identifiers
	char c= *(ach++);
	if ((c>='a' && c<='z') || (c>='A' && c<='Z') || c=='_')
	{
		while (((c= *(ach))>='a' && c<='z') ||
				 (c>='A' && c<='Z') ||
				 (c>='0' && c<='9') ||
				 c=='_')
			++ach;

		return ach;
	}
	return 0; // Not a valid identifier
}

