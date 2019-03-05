/***********************************************************************************
 *
 *   singlelist.h -- Singly-linked list data type 
 *
 *   Erik Persson, Sketchware, 2002
 */

// Singly-linked list data type
// InsertLast is o(1)

#ifndef TYPES_SINGLELIST_H
#define TYPES_SINGLELIST_H

template <class SNode>
class CSingleList
{
protected:
	SNode *m_pFirst;
	SNode *m_pLast;

public:
	CSingleList(void) { m_pFirst= m_pLast= NULL; }
	~CSingleList(void) { Clear(); }
	
	inline SNode *GetFirst() const;
	inline SNode *GetLast() const;
	inline void InsertFirst(SNode *pNode);
	inline void InsertLast(SNode *pNode);
	inline void Clear();
	inline void DeleteNode(SNode *pDelete);
	inline int GetSize() const;
	inline int Truncate(int nCount);
	inline bool IsEmpty() const { return !m_pFirst; }
};

template <class SNode>
SNode *CSingleList<SNode>::GetFirst() const
{
	return m_pFirst;
}

template <class SNode>
SNode *CSingleList<SNode>::GetLast() const
{
	return m_pLast;
}

template <class SNode>
void CSingleList<SNode>::InsertFirst(SNode *pNode)
{
	pNode->m_pNext= m_pFirst;
	m_pFirst= pNode;

	if (!m_pLast)
		m_pLast= pNode;
}

template <class SNode>
void CSingleList<SNode>::InsertLast(SNode *pNode)
{
	pNode->m_pNext= NULL;
	if (!m_pFirst)
		m_pFirst= pNode;

	if (m_pLast)
		m_pLast->m_pNext= pNode;

	m_pLast= pNode;
}

template <class SNode>
void CSingleList<SNode>::Clear(void)
{
	SNode *pNode;
	while (pNode= m_pFirst)
	{
		m_pFirst= pNode->m_pNext;
		delete pNode;
	}
	m_pLast= NULL;
}

template <class SNode>
void CSingleList<SNode>::DeleteNode(SNode *pDelete)
{
	SNode *pNode, **ppNode, *pPrev;

	pPrev= NULL;
	ppNode= &m_pFirst;
	while (pNode= *ppNode)
	{
		if (pNode==pDelete)
		{
			// Update the pointer that pointed to this
			*ppNode= pNode->m_pNext;

			// Update m_pLast 
			if (m_pLast == pNode)
				m_pLast= pPrev;

			delete pNode;
			return;
		}
		pPrev= pNode;
		ppNode= &pNode->m_pNext;
	}
}

template <class SNode>
int CSingleList<SNode>::GetSize() const
{
	int n=0;
	SNode *pL= m_pFirst;
	while (pL)
	{
		n++;
		pL= pL->m_pNext;
	}
	return n;
}

template <class SNode>
int CSingleList<SNode>::Truncate(int nCount)
{
	// Truncate a list, i.e. delete all nodes 
	// except the nCount first. Returns the number
	// of nodes left (will be less or eq to nCount).
	int n= 0;
	SNode *pNprev= NULL;
	SNode *pN= GetFirst();
	while (pN && n<nCount)
	{
		n++;
		pNprev= pN;
		pN= pN->m_pNext;
	}

	if (pN)
	{
		pNprev->m_pNext= NULL;
		while (pN)
		{
			SNode *pNtmp= pN->m_pNext;
			delete pN;
			pN= pNtmp;
		}
	}

	return n;
}

#endif // TYPES_SINGLELIST_H