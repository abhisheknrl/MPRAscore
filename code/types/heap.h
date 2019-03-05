/*******************************************************************************
 *
 *  heap.h -- Template for min sorted heaps
 *
 *  Björn Nilsson, Veta Mera Data HB, May 2004 (original version March 2002)
 *
 *	T must have
 *	
 *	1. A defined less-than operator <
 *	2. A public member named int m_nHeapIndex
 *	3. A well-defined copy constructor
 *
 *	NOTE:
 *	It is required that you provide a "negative
 *	infinity element" that is smaller than any 
 *	other occuring element.
 *
 *	NOTE:
 *	The heap only contains pointers to T elements
 */

#ifndef HEAP_H
#define HEAP_H

#include <afxwin.h>

template<class T> class CHeap
{
protected:
	// Data pointer
	T **m_ppHeap; 	

	// Maximum heap size.
	int m_nHeapSize;
	int m_nMaxHeapSize;

	void Free();
	void Quicksort(T **pLeft, T **pRight);
public:
	CHeap() { m_ppHeap= NULL; m_nHeapSize=0; m_nMaxHeapSize=0; };
	~CHeap() { Free(); }
	T **GetHeap() { return m_ppHeap; }
	bool SetSize(int nSize);
	void SetNegInfinity(T *pT);
	void AssertHeapProperty(); // Debug function

	int GetHeapSize() const { return m_nHeapSize; }
	bool IsEmpty() const { return m_nHeapSize>0; }
	T *GetRoot() { return m_ppHeap[1]; }
	void Sort();

	void Add(T *pT)
	{
		// Add and treckle into heap.
		ASSERT(pT);
		ASSERT(m_nHeapSize < m_nMaxHeapSize);

		int i0= ++m_nHeapSize;
		int i1;
		while ( (*pT) < (*m_ppHeap[i1= i0>>1]) )
		{
			(m_ppHeap[i0]= m_ppHeap[i1])->m_nHeapIndex= i0;
			i0= i1;
		}
		(m_ppHeap[i0]= pT)->m_nHeapIndex= i0;
	}

	void AddFast(T *pT)
	{
		// Add without sort (must be followed by 
		// sorting of the heap before using operations 
		// that assume that the heap property is true).
		int i0= ++m_nHeapSize;
		m_ppHeap[i0]= pT;
		pT->m_nHeapIndex= i0;		
	}

	void BringUp(int i0)
	{
		T *pE= m_ppHeap[i0];

		int i1;
		while ( (*pE) < (*m_ppHeap[i1= i0>>1]) )
		{
			(m_ppHeap[i0]= m_ppHeap[i1])->m_nHeapIndex= i0;
			i0= i1;
		}
		(m_ppHeap[i0]= pE)->m_nHeapIndex= i0;
	}

	void BringDown(int i0)
	{
		T *pT= m_ppHeap[i0];
		int i1= i0;

		// This neat construction works because the heap 
		// is 1-based instead of 0-based. The 0 position
		// contains a pointer to the neg-inf element.
		while ((i1<<=1)<m_nHeapSize) 
		{
			if (*m_ppHeap[i1+1] < *m_ppHeap[i1]) 
				i1++;

			if (*m_ppHeap[i1] < (*pT) )
			{
				(m_ppHeap[i0]= m_ppHeap[i1])->m_nHeapIndex= i0;
				i0= i1;
			}
			else
				break;
		}
		(m_ppHeap[i0]= pT)->m_nHeapIndex= i0;
	}

	void Delete(int i0)
	{
		ASSERT(i0>0 && i0<=m_nHeapSize);

		m_ppHeap[i0]->m_nHeapIndex= 0;
		T *pLast= m_ppHeap[m_nHeapSize--];
		int i1= i0;
		while ((i1<<=1)<m_nHeapSize) 
		{
			if (*m_ppHeap[i1+1] < *m_ppHeap[i1]) 
				i1++;

			if (*m_ppHeap[i1] < *pLast)
			{
				(m_ppHeap[i0]= m_ppHeap[i1])->m_nHeapIndex= i0;
				i0= i1;
			}
			else
				break;
		}
		(m_ppHeap[i0]= pLast)->m_nHeapIndex= i0;
	}
};

template<class T>
void CHeap<T>::Free()
{
	if (m_ppHeap)
		delete [] m_ppHeap;
	m_ppHeap= NULL;
}

template<class T>
bool CHeap<T>::SetSize(int nSize)
{
	Free();
	if (nSize>0)
	{
		try
		{
			// Allocate nSize+1 elements (the extra 
			// element is used for neginf pointer)
			m_ppHeap= new T*[nSize+1];
			m_nMaxHeapSize= nSize;
			m_nHeapSize= 0;
			return true;
		}
		catch(...)
		{
			return false;
		}
	}
	else
		return true;
}

template<class T>
void CHeap<T>::SetNegInfinity(T *pT)
{
	if (!m_ppHeap)
	{
		// Heap must have been allocated (using
		// SetSize) prior to this call
		ASSERT(false); 
		return;
	}

	m_ppHeap[0]= pT;
}

template<class T>
void CHeap<T>::Quicksort(T **pLeft, T **pRight)
{
	if (pRight<=pLeft)
		return;

	T **pi= pLeft-1;
	T **pj= pRight+1;
	T *v= *pLeft;
	T *tmp;
	do 
	{
		while (**(++pi)<(*v));
		while ((*v)<**(--pj));

		tmp= *pi;
		*pi= *pj;
		*pj= tmp;
	}
	while (pi<pj);

	tmp= *pi;
	*pi= *pj;
	*pj= tmp;

	Quicksort(pLeft, pj);
	Quicksort(pj+1, pRight);
}

template<class T>
void CHeap<T>::Sort()
{
	if (!m_ppHeap || m_nHeapSize<2)
		return;

	T **pLeft= &m_ppHeap[1];
	T **pRight= &m_ppHeap[m_nHeapSize];
	Quicksort(pLeft, pRight);
}

template<class T>
void CHeap<T>::AssertHeapProperty()
{
#ifdef _DEBUG
	for (int i=1;i<=m_nHeapSize;i++)
	{
		ASSERT( !(*m_ppHeap[i] < *m_ppHeap[i>>1]) );
		ASSERT( m_ppHeap[i]->m_nHeapIndex==i );
	}
#endif
}

#endif
