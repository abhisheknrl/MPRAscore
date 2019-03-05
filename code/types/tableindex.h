/*******************************************************************************
 *
 *   tableindex.h -- Indexes memory containing a tab-delimited text dataset.
 *
 *   Björn Nilsson, 2004-2008
 */

#ifndef TABLEINDEX_H
#define TABLEINDEX_H

#include "vector.h"
#include "nan.h"

class CTableIndexRow
{
public:
	CVector<const char *> m_vFirst;
	CVector<const char *> m_vLast;

	CTableIndexRow() {};
	CTableIndexRow(CTableIndexRow &dst) // Copy ctor
	{
		dst.m_vFirst= m_vFirst;
		dst.m_vLast= m_vLast;
	}

	int GetColumnCount() const
	{ 
		ASSERT(m_vFirst.GetSize()==m_vLast.GetSize()); 
		return m_vFirst.GetSize(); 
	}

	void Extend(int nNewCols, const char *aDummy); // Used by MakeRectangular
};


#ifdef _WIN32
typedef void (__stdcall *PARSECALLBACK)(const int);
#else
typedef void __attribute__((stdcall)) *PARSECALLBACK(const int);
#endif

class CTableIndex
{
protected:
	char m_cSeparator;
	char m_cQuote;
public:
	CVector<CTableIndexRow> m_vRows;

	virtual ~CTableIndex() {};
	CTableIndex() {};
	CTableIndex(const CTableIndex &src) // Copy ctor
	{
		m_cSeparator= src.m_cSeparator;
		m_cQuote= src.m_cQuote;
		m_vRows= src.m_vRows;
	}

	void Clear();

	bool ScanTDS(const char *ach, PARSECALLBACK pCallback=NULL);
	bool ScanCSV(const char *ach, const char cComma, PARSECALLBACK pCallback=NULL);
	bool ScanSpaceSeparated(const char *ach, PARSECALLBACK pCallback);

	bool IsEmpty() const;
	bool IsRectangular() const;

	int GetRowCount() const { return m_vRows.GetSize(); }
	void GetColCount(int &min, int &max) const;
	int GetColCount(int nRow) const { return m_vRows[nRow].GetColumnCount(); }
	int GetMinColCount() const { int min, max; GetColCount(min,max); return min; }
	int GetMaxColCount() const { int min, max; GetColCount(min,max); return max; }
	const CTableIndexRow &GetRow(int nRow) const { return m_vRows[nRow]; }
	char GetSeparatorChar() const { return m_cSeparator; }
	char GetQuoteChar() const { return m_cQuote; }
	bool Transpose(CTableIndex &ti1);

	const char *GetFirstPointer(int nRow, int nCol) const
	{ 
		ASSERT (nRow>=0 && nCol>=0 && nRow<m_vRows.GetSize() && nCol<m_vRows[nRow].GetColumnCount());
		return m_vRows[nRow].m_vFirst[nCol];
	}
	const char *GetLastPointer(int nRow, int nCol) const
	{
		ASSERT (nRow>=0 && nCol>=0 && nRow<m_vRows.GetSize() && nCol<m_vRows[nRow].GetColumnCount());
		return m_vRows[nRow].m_vLast[nCol];
	}	

	void MakeRectangular();
	bool CreateResorted(const CTableIndex &ti, int m_nAnnotRows, CVector<int> vOrder);
};

#endif
