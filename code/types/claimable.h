/***********************************************************************************
 *
 *   claimable.h -- CClaimable class
 *
 *   Erik Persson and Björn Nilsson, 1999-2008
 */

#ifndef CLAIMABLE_H
#define CLAIMABLE_H

// This class implements reference count
// It can be used with multiple inheritance

#include "types/assert.h"

class CClaimable
{
public:
	int m_nRefs;

	void Claim()
	{
		ASSERT(m_nRefs>0);
		++m_nRefs;
	}
	void Release() 
	{
		ASSERT(m_nRefs>0);
		if (!--m_nRefs) delete this;
	}

	CClaimable() { m_nRefs= 1; }

protected:
	// Destructor is protected to prevent use of delete operator
	virtual ~CClaimable() { ASSERT( m_nRefs == 0 ); }
};

#endif // CLAIMABLE_H
