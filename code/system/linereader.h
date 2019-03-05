//
// linereader.h -- For reading huge tab-delimited files line-by-line very efficiently
//

#ifndef LINEREADER_H
#define LINEREADER_H

#include <stdio.h>
#include "types/vector.h"

#define LINEREADER_QUANTA (20000000)

class CLineReader
{
protected:
	FILE *m_pf;
	CVector<char> m_vBuf;
	char *m_aBuf;
	size_t m_nRead;
	size_t m_nPos;

// #ifdef _WIN32
public:
	long long m_nRead_sum;
	long long m_nRead_filesize;
// #endif

public:
	CLineReader();
	~CLineReader();
	bool Open(const char *aFilename);
	bool GetLine(char *aLineBuf, size_t nBufSize, CVector<const char *> &vItems, int &nItems);
	bool GetLine(char *aLineBuf, size_t nBufSize, CVector<const char *> &vItems, int &nItems, char sep1, char sep2, bool bMulti);
	bool Close();
	float GetProgress();
};

#endif
