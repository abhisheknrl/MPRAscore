/***********************************************************************************
 *   
 *	 stringfun.h -- String utility functions
 *
 *   Erik Persson and Björn Nilsson, 1998-
 *
 */

#ifndef STRINGFUN_H
#define STRINGFUN_H

#include <stdlib.h> // for atoi, atof
#include "vector.h"
#include "types/string.h"

CString ToLower(CString);
CString ToUpper(CString);

CString Format(const char *ach, ...);
CString FormatVA(const char *ach, va_list argptr);

const char *strreverse(char *ach);
int strcmp_nocase(const char *a0, const char *ach2);
bool strmatch(const char *a0, const char *aF);
bool strmatch_nocase(const char *a0, const char *aF);
int strfind(const char *a0, const CVector<CString> &v);
int strfind(const char *a0, const CVector<const char *> &v);
int strfind_nocase(const char *a0, const CVector<CString> &v);
int strfind_nocase(const char *a0, const CVector<const char *> &v);
int strinset(const char *a0, const char *aset, const char csep);


CString BumpName(CString sBaseName);
CString GetUniqueName(CString sBaseName, void *pObj, bool (*nameIsUsed)(void *, const char *) );
bool ParseInt(const char *ach, int *pResult);
bool MatchPrefixNoCase(const char *ach, const char *aPrefix);

CString ToForwardSlash(CString s);
CString ToBackSlash(CString s);
CString ToLinefeedOnly(CString s);
CString ToCrlf(CString s);
int GetColumnCount(const char *ach);

// Added 2008-5 
void SplitAtChar(const CString &s, const char c, CVector<CString> &vS);

// Added 2008-05-30, moved from lingo, stx, and tll
const char *Scan_Int(const char *ach, int *presult);
const char *Scan_UnsignedInt(const char *ach, int *presult);
const char *Scan_Ident(const char *ach);
inline bool IsBinDigit(char c) { return c >= '0' && c <= '1'; }
inline int GetBinDigit(char c) { return c-'0'; }
inline bool IsDecDigit(char c) { return c >= '0' && c <= '9'; }
inline int GetDecDigit(char c) { return c-'0'; }
inline bool IsHexDigit(char c) { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'); }
inline int GetHexDigit(char c) 
{
	if (c >= 'a' && c <= 'f')
		return c-'a'+10;
	if (c >= 'A' && c <= 'F')
		return c-'A'+10;
	return c-'0';
}

inline const char *Scan_SkipSpace(const char *ach)
{
	// spc ::= {space_char}
	// space_char ::= { VT | SPC | SSPC }
	char c;
	while ((c= *ach)==32 || c== (char) 160 || c==9)
		ach++;
	return ach;
}

inline const char *Scan_SkipLine(const char *ach)
{
	// Line ::= non_CR_LF_chars (LF | CR {LF})
	// TODO: handle escaped EOL
	char c;
	while ((c= *ach) != 0 && c != 13 && c != 10)
	{	
		ach++;
		c= *ach;
	}

	if (c==13)
	{
		ach++;
		c= *ach;
		while (c==10)
		{	
			ach++;
			c= *ach;
		}
	}
	else if (c==10)
		ach++;

	return ach;
}

inline const char *Scan_Eol(const char *ach)
{
	// Comments are handled here also
	// TODO: handle escaped EOL
	// EOL ::= LF | CR {LF} | "#" chars [EOL]
	char c;
	c= *ach;
	if (/* c=='#' // SCL comment 
		 || */ c==13
		 || c==10)
		 return Scan_SkipLine(ach);
	return NULL; // Not an EOL
}
inline const char *Scan_SkipRichSpace(const char *ach)
{
  // This function considers EOL to be ordinary whitespace
  // space ::= { space_char | EOL }
  // space_char ::= { VT | SPC | SSPC }

	ach= (const char *)Scan_SkipSpace(ach);
	while (const char *a1= (const char *)Scan_Eol(ach))
		ach= (const char *)Scan_SkipSpace(a1);
	return ach;
}


#endif // STRINGFUN_H
