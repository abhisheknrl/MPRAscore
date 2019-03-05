/*******************************************************************************
 *
 *   annotation.h -- 
 *
 *   Björn Nilsson, Medical Mathematics, 2005-2006
 */

#ifndef ANNOTATION_H
#define ANNOTATION_H

#include <afx.h>
#include "types/vector.h"
#include "types/ref.h"

/******************************************************************************/
// CColorMap - maps string values to rgb values

class CColorMap
{
protected:
	CVector<CString> m_vsKey;
	CVector<COLORREF> m_vcValue;

public:
	CColorMap() {};
	CColorMap(const CColorMap& other);
	~CColorMap() {};

	COLORREF GetColor(const char *aProp) const;
	bool SetColor(const char *aProp, COLORREF aVal);	
	void Clear()
	{
		m_vsKey.SetSize(0);
		m_vcValue.SetSize(0);
	}

	int GetCount() const;
	int FindKey(const char *aProp) const;
	CString GetKey(int n) const;
	COLORREF GetValue(int n) const;
	bool AppendAssoc(const char *aProp, COLORREF aCol);
};

/******************************************************************************/
// CAnnotation

// Annotation flags
#define ANF_Text (1)
#define ANF_Colorize (2)
#define ANF_Box (4)
#define ANF_Dot (8)
#define ANF_Default (ANF_Text)

class CAnnotation 
{
public:
	CString m_sTitle;
	int m_nFlags;
	CColorMap m_Color;
	CVector<CString> m_vValues;

	CAnnotation() { m_nFlags= ANF_Default; }
};

#endif