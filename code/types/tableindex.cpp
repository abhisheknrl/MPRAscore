/*******************************************************************************
 *
 *   tableindex.cpp -- 
 *
 *   Björn Nilsson, 2004-2008
 */

#include "tableindex.h"
#include "stringfun.h"

/******************************************************************************/

void CTableIndexRow::Extend(int nNewCols, const char *aDummy)
{
	ASSERT(nNewCols>=0);
	if (nNewCols<0)
		return;
	ASSERT(m_vFirst.GetSize()==m_vLast.GetSize());
	int nOldCols= m_vFirst.GetSize();
	for (int i=nOldCols;i<nNewCols;i++)
	{
		m_vFirst.Add(aDummy);
		m_vLast.Add(aDummy);
	}
}

/******************************************************************************/

bool CTableIndex::IsEmpty() const
{
	if (!GetRowCount())
		return true;
	return !m_vRows[0].m_vFirst.GetSize();
}

bool CTableIndex::IsRectangular() const
{
	// Returns true if table has the same number of elements on all rows.
	int min, max;
	GetColCount(min, max);
	return (min==max);
}

void CTableIndex::GetColCount(int &min, int &max) const
{
	// Returns smallest and largest number of elements per indexed row.
	// The two will co-incide for rectangular (most common case) matrices.
	if (IsEmpty())
	{
		min= max= 0;
		return;
	}

	min= max= m_vRows[0].GetColumnCount();
	for (int i=GetRowCount();--i>0;)
	{
		int j= m_vRows[i].GetColumnCount();
		if (j<min)
			min= j;
		else if (j>max)
			max= j;
	}
}

bool CTableIndex::ScanCSV(const char *ach, const char cComma, PARSECALLBACK pCallback)
{
	// Create table index from comma-separated values (CSV) format

	Clear();

	CTableIndexRow tir;
	int r= 0;
	while (*ach)
	{
		int c=0;
		while (true)
		{
			while (*ach==32 || (unsigned char)(*ach)==160 || *ach==9)
				ach++;
			const char *ach0;
			const char *ach1;
			if (*ach=='\"')
			{
				// Quoted element
				ach0= ++ach;
				while (*ach && *ach!='\"') // EOL allowed inside quotes in CSV format
					ach++;
				ach1= ach;

				// TODO: Handle embedded "" 

				if (*ach!='\"')
					return false; // EOF inside quotes
				ach++;
			}
			else
			{
				// Non-quoted element
				ach0= ach;
				while (*ach && *ach!=10 && *ach!=13 && *ach!=cComma)
					ach++;
				ach1= ach;
			}

			// Ignore leading and trailing spaces
			while (ach0<ach && (*ach0==32 || *ach0==9 || (unsigned char)(*ach0)==160))
				ach0++;
			while (ach0<ach1 && (*(ach1-1)==32 || *(ach1-1)==9 || (unsigned char)(*(ach1-1))==160))
				ach1--;

			if (c<tir.m_vFirst.GetSize())
			{
				tir.m_vFirst.SetAt(c, ach0);
				tir.m_vLast.SetAt(c, ach1);
			}
			else
			{
				tir.m_vFirst.Add(ach0);
				tir.m_vLast.Add(ach1);
			}
			c++;

			while (*ach==32 || (unsigned char)(*ach)==160 || *ach==9)
				ach++;
			if (*ach==0 || *ach==10 || *ach==13)
				break;
			if (*ach!=cComma)
				return false;
			ach++;
		}

		tir.m_vFirst.SetSize(c);
		tir.m_vLast.SetSize(c);
		m_vRows.Add(tir);

		if (pCallback)
			pCallback(r);
		r++;
		ach= Scan_SkipLine(ach);
	}

	return true;
}

bool CTableIndex::ScanTDS(const char *ach, PARSECALLBACK pCallback)
{
	// Create table index from tab-delimited, non-quoted data
	// CVector<const char*> vFirst;
	// CVector<const char*> vLast;

	Clear();

	CTableIndexRow tir;

	int r= 0;
	while (*ach)
	{
		int c=0;
		while (true) 
		{
			// Delimit element
			const char *ach0= ach;
			while (*ach && *ach!=9 && *ach!=10 && *ach!=13)
				ach++;
		
			// Ignore leading and trailing spaces
			while (ach0<ach && (*ach0==32 || (unsigned char)(*ach0)==160))
				ach0++;
			const char *ach1= ach;
			while (ach0<ach1 && (*(ach1-1)==32 || (unsigned char)(*(ach1-1))==160))
				ach1--;

			// Remove quotes and spaces immediately inside the quotes
			if (ach0+1<ach1 && *ach0=='\"' && *(ach1-1)=='\"')
			{
				ach0++;	ach1--;
				while (ach0<ach && (*ach0==32 || (unsigned char)(*ach0)==160))
					ach0++;
				const char *ach1= ach;
				while (ach0<ach1 && (*(ach1-1)==32 || (unsigned char)(*(ach1-1))==160))
					ach1--;
			}

			// BN 2018-09-27
			// convert any interior space 160 (sometimes Excel outputs these) to space 32
			char *ach2= (char *)ach0;
			while (ach2<ach1)
			{
				if ((unsigned char)(*ach2)==160)
					*ach2= 32;
				ach2++;
			}

			if (c<tir.m_vFirst.GetSize())
			{
				tir.m_vFirst.SetAt(c, ach0);
				tir.m_vLast.SetAt(c, ach1);
			}
			else
			{
				tir.m_vFirst.Add(ach0);
				tir.m_vLast.Add(ach1);
			}
			c++;

			if (*ach==9)
				ach++;
			else if (*ach==0 || *ach==10 || *ach==13)
				break;
		} 

		tir.m_vFirst.SetSize(c);
		tir.m_vLast.SetSize(c);
		m_vRows.Add(tir);

		if (pCallback)
			pCallback(r);
		r++;
		ach= Scan_SkipLine(ach);
	}

	return true;
}

bool CTableIndex::ScanSpaceSeparated(const char *ach, PARSECALLBACK pCallback)
{
	// Create table index from tab-delimited, non-quoted data
	// CVector<const char*> vFirst;
	// CVector<const char*> vLast;

	Clear();

	CTableIndexRow tir;

	int r= 0;
	while (*ach)
	{
		int c=0;
		while (true) 
		{
			// Delimit element
			const char *ach0= ach;
			while (*ach && *ach!=32 && (unsigned char)(*ach)!=160 && *ach!=10 && *ach!=13)
				ach++;
			const char *ach1= ach;

			/*
			// Remove quotes and spaces immediately inside the quotes
			if (ach0+1<ach1 && *ach0=='\"' && *(ach1-1)=='\"')
			{
				ach0++;	ach1--;
				while (ach0<ach && *ach0==32)
					ach0++;
				const char *ach1= ach;
				while (ach0<ach1 && *(ach1-1)==32)
					ach1--;
			}
			*/

			if (c<tir.m_vFirst.GetSize())
			{
				tir.m_vFirst.SetAt(c, ach0);
				tir.m_vLast.SetAt(c, ach1);
			}
			else
			{
				tir.m_vFirst.Add(ach0);
				tir.m_vLast.Add(ach1);
			}
			c++;

			if (*ach==32 || (unsigned char)(*ach)==160)
				ach++;
			else if (*ach==0 || *ach==10 || *ach==13)
				break;
		} 

		tir.m_vFirst.SetSize(c);
		tir.m_vLast.SetSize(c);
		m_vRows.Add(tir);

		if (pCallback)
			pCallback(r);
		r++;
		ach= Scan_SkipLine(ach);
	}

	return true;
}

bool CTableIndex::Transpose(CTableIndex &dst)
{
	int R= GetRowCount();
	int C= GetMinColCount();

	try
	{
		dst.m_vRows.SetSize(C);
		for (int c=0;c<C;c++)
		{
			CTableIndexRow tir;
			tir.m_vFirst.SetSize(R);
			tir.m_vLast.SetSize(R);

			for (int r=0;r<R;r++)
			{
				tir.m_vFirst.SetAt(r, GetFirstPointer(r, c));
				tir.m_vLast.SetAt(r, GetLastPointer(r, c));
			}

			dst.m_vRows.SetAt(c, tir);
		}
	}
	catch(...)
	{
		return false; // Out of memory.
	}

	return true;
}

void CTableIndex::MakeRectangular()
{
	int cmax= GetMaxColCount();
	for (int r=0;r<GetRowCount();r++)
		if (GetColCount(r)<cmax)
		{
			CTableIndexRow tir;
			tir= m_vRows.GetAt(r);
			tir.Extend(cmax, GetFirstPointer(0,0));
			m_vRows.SetAt(r, tir);
		}

	ASSERT(IsRectangular());
}

void CTableIndex::Clear()
{
	m_vRows.SetSize(0);
}

/******************************************************************************/

bool CTableIndex::CreateResorted(const CTableIndex &ti, int nAnnotRows, CVector<int> vOrder)
{
	// Creates an index the body of consists of the rows in ti 
	// arranged in the order given in vOrder, except for  the 
	// first m_nAnnotRows rows which are copied as they are.

	m_vRows.SetSize(nAnnotRows+vOrder.GetSize());

	// Annot rows
	for (int i=0;i<nAnnotRows;i++)
		m_vRows.SetAt(i, ti.m_vRows[i]);

	// Reordered body
	for (int i=0;i<vOrder.GetSize();i++)
	{
		int src= vOrder[i];
		if (src<0 || src>=ti.GetRowCount()-nAnnotRows)
		{
			ASSERT(false);
			return false;
		}

		m_vRows.SetAt(i + nAnnotRows, ti.m_vRows[src + nAnnotRows]);
	}

	return true;
}
