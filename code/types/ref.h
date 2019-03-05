/***********************************************************************************
 *
 *   ref.h -- Automatically Claiming Reference
 *
 *   Erik Persson and Björn Nilsson, 2001-2008
 */

#ifndef REF_H
#define REF_H

#ifndef CLAIMABLE_H
#include "claimable.h"
#endif

extern int g_nRefs;

template <class T>
class CRef
{
	mutable T *m_pObj;
public:
	CRef(T *p= 0)
	{
		g_nRefs++;
		m_pObj= p;
		if (m_pObj) m_pObj->Claim();
	}
	/*
	template <class T2>
	CRef(const CRef<T2>& other)
	{
		g_nRefs++;
		m_pObj= (T2 *) (const T2 *) other; // static type check
		if (m_pObj) m_pObj->Claim();
	}*/

	// This must exist and be nontemplate 
	// to prevent default copy construction
	CRef(const CRef<T>& other)
	{
		g_nRefs++;
		m_pObj= other.m_pObj;
		if (m_pObj) m_pObj->Claim();
	}
	
	~CRef()
	{	
		g_nRefs--;
		if (m_pObj) m_pObj->Release();
	}

	// TODO: Få bort /EP
	operator CRef<CClaimable>() const
	{
		return CRef<CClaimable>((CClaimable *) m_pObj);
	}

	operator const T *() const
	{	return m_pObj;
	}

	const T& operator*() const
	{	return *m_pObj;
	}

	T& CreateNew() // Might throw exception
	{
		T *pNew= new T;
		*this= pNew;
		pNew->Release();
		return *m_pObj;
	}

	T& CreateCopy() // Might throw exception
	{
		T *pNew;
		if (m_pObj)
			pNew= new T(*m_pObj);
		else
			pNew= new T;
		*this= pNew;
		pNew->Release();
		return *m_pObj;
	}

	// Works like GetExclusive() if exclusive reference, 
	// works like CreateCopy() otherwise.
	T& CreateExclusive() // Might throw exception
	{
		T *pNew;
		if (m_pObj)
		{
			if (m_pObj->m_nRefs == 1)
				return *m_pObj;
			else
				pNew= new T(*m_pObj);
		}
		else
			pNew= new T;
		*this= pNew;
		pNew->Release();
		return *m_pObj;
	}

	T& GetExclusive() const
	{
		// TODO: Improve /EP
		if (m_pObj)
		{
			ASSERT(m_pObj->m_nRefs==1);
		}
		return *m_pObj;
	}

	bool IsExclusive() const
	{
		return m_pObj->m_nRefs==1;
	}

	const T *operator->() const
	{
		return m_pObj;
	}

	bool operator==(const CRef<T>& other) const
	{	return m_pObj == other.m_pObj;
	}

	bool operator!=(const CRef<T>& other) const
	{	return m_pObj != other.m_pObj;
	}

	void SetPtr(T *p)
	{
		if (m_pObj) m_pObj->Release();
		m_pObj= p;
		if (m_pObj) m_pObj->Claim();
	}

	/* Experimental cast
	template <class T2>
	operator CRef<T2>() const
	{
		return CRef<T2>((T2 *)(const T2 *) m_pObj); // static type check is performed here
	}
	*/

	// This must exist and be nontemplate 
	// to prevent default copy construction
	const CRef<T>& operator=(const CRef<T>& other)
	{
		SetPtr(other.m_pObj);
		return *this;
	}

	template <class T2>
	inline CRef<T>& operator=(T2 *p)
	{
		SetPtr(p); // static type check is performed here
		return *this;
	}


/*	// Templatized version as separate function
	template <class T2>
	inline const CRef<T>& operator=(const CRef<T2>& other)
	{
		SetPtr((T2 *) (const T2 *) other); // static type check is performed here
		return *this;
	}*/
/*
	const CRef<T>& operator=(T *p)
	{
		if (m_pObj) m_pObj->Release();
		m_pObj= p; // static type check is performed here
		if (m_pObj) m_pObj->Claim();
		return *this;
	}
	const CRef<T>& operator=(const CRef<T>& other)
	{
		if (m_pObj) m_pObj->Release();
		m_pObj= other.m_pObj; // static type check is performed here
		if (m_pObj) m_pObj->Claim();
		return *this;
	}
*/

	/*
	const CRef<T>& operator=(const CRef<T>& other)
	{
		return ((*this)= other.m_pObj);
	}*/
};

/* Templatized version as separate function
template <class T, class T2>
inline const CRef<T>& operator=(CRef<T>& ref, const CRef<T2>& other)
{
	ref.SetPtr((T2 *) (const T2 *) other); // static type check is performed here
	return ref;
}*/

// Moved here from EaseParser.cpp
template <class T>
CRef<T> Create(void)
{
	CRef<T> rObj;
	try
	{
		T *pObj;
		pObj= new T;
		rObj= pObj;
		pObj->Release();
	}
	catch(...)
	{
		// new T failed - rObj will be NULL
	}
	return rObj;
}

#endif // REF_H
