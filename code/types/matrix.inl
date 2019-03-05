/*******************************************************************************
 *
 *   matrix.inl -- 
 *
 */

// #include <cassert>
#include <memory>
#include <strstream>

//-----------------------------------------------------------------------------------------

template<class T>
CMatrix<T>::CMatrix()
: m_pElems(NULL), m_width(0), m_height(0)
{

}
//-----------------------------------------------------------------------------------------

template<class T>
CMatrix<T>::CMatrix(const int nRows, const int nCols)
: m_pElems(NULL), m_width(0), m_height(0)
{
	ASSERT(nRows > 0 && nCols > 0);
	m_pElems = new T[nRows*nCols];
	m_width = nCols;
	m_height = nRows;
	ASSERT(m_pElems);
}
//-----------------------------------------------------------------------------------------

template<class T>
CMatrix<T>::CMatrix(const int nRows, const int nCols, const T initvalue)
: m_pElems(NULL), m_width(0), m_height(0)
{
	ASSERT(nRows > 0 && nCols > 0);
	m_pElems = new T[nRows*nCols];
	m_width = nCols;
	m_height = nRows;
	ASSERT(m_pElems);
	InitValues(initvalue);
}
//-----------------------------------------------------------------------------------------

/*
template<class T>
CMatrix<T>::CMatrix(const int nRows, const int nCols, const char * pszInitValues)
: m_pElems(NULL), m_width(0), m_height(0)
{
	ASSERT(nRows > 0 && nCols > 0);
	int nSize = nRows*nCols;
	m_pElems = new T[nSize];
	m_width = nCols;
	m_height = nRows;
	ASSERT(m_pElems);
//	SetValues(0); // init to zero first ???
	std::istrstream buffer(pszInitValues, strlen(pszInitValues));
	double str2num;
	for(int i = 0; i < nSize; ++i)
	{
		buffer >> str2num; // ISSUE: What happens if the the stream i smaller than the number of elements??
		m_pElems[i] = static_cast<T>(str2num);
	}
}
//-----------------------------------------------------------------------------------------
*/

template<class T>
CMatrix<T>::CMatrix(const CMatrix<T> & other)
: m_pElems(NULL), m_width(0), m_height(0)
{
	if(other.m_pElems)
	{
		m_width = other.m_width;
		m_height = other.m_height;
		m_pElems = new T[m_width*m_height];
		ASSERT(m_pElems);
		memcpy(m_pElems, other.m_pElems, m_width*m_height*sizeof(T));
	}
	// else nothing to copy! - Leave it be
}
//-----------------------------------------------------------------------------------------

template<class T>
CMatrix<T>& CMatrix<T>::operator=(const CMatrix<T>& other)
{
	ASSERT(other.m_pElems && other.m_width>0 && other.m_height>0);
	if(this != &other)
	{
		if(m_width != other.m_width || m_height != other.m_height)
		{
			if(m_pElems)
			{
				delete[] m_pElems;
				m_pElems = NULL;
			}
			m_pElems = new T[other.m_width*other.m_height];
			ASSERT(m_pElems);
			m_width = other.m_width;
			m_height = other.m_height;
		}
		ASSERT(m_pElems);
		memcpy(m_pElems, other.m_pElems, m_width*m_height*sizeof(T));
	}
	return *this;
}
//-----------------------------------------------------------------------------------------

template<class T>
CMatrix<T>::~CMatrix()
{
	if(m_pElems)
	{
		delete[] m_pElems;
		m_pElems = NULL;
	}
}
//-----------------------------------------------------------------------------------------

template<class T>
void CMatrix<T>::ReInit(const int NRows, const int NCols)
{
	ASSERT(NRows > 0 && NCols > 0);
	if(m_pElems)
	{
		delete[] m_pElems;
		m_pElems = NULL;
	}
	m_pElems = new T[NRows*NCols];
	ASSERT(m_pElems);
	m_width = NCols;
	m_height = NRows;
}
//-----------------------------------------------------------------------------------------

template<class T>
void CMatrix<T>::ReInit(const int NRows, const int NCols, const T initValue)
{
	ReInit(NRows, NCols);
	InitValues(initValue);
}
//-----------------------------------------------------------------------------------------

template<class T>
void CMatrix<T>::InitValues(const T value)
{
	ASSERT(m_pElems && m_width > 0 && m_height > 0);
	if(static_cast<double>(value) == 0.0)
	{
		// Then we can use memset (fast!)
		memset(m_pElems, 0, m_width*m_height*sizeof(T));
	}
	else
	{
		// we have to go about it the slow way, since memset is setting each byte to the value, 
		// datatypes using more than one byte will not get the correct value from memset if != 0
		T 
			* pP = m_pElems, 
			* pL = m_pElems + m_width*m_height; // will actually jump w*h*sizeof(T) bytes
		while(pP < pL)
		{
			*pP++ = value;
		} // for each element
	} // if equal to 0	
}
//-----------------------------------------------------------------------------------------

template<class T>
T * CMatrix<T>::GetFirstPtr() const
{
	return m_pElems;
}
//-----------------------------------------------------------------------------------------

template<class T>
void CMatrix<T>::GetFirstLastPtrs(T ** pFirst, T ** pLast) const
{
	*pFirst = m_pElems;
	if(!m_pElems)
	{
		*pLast = NULL; // pFirst == pLast == NULL;
	}
	else
	{
		*pLast = m_pElems + m_width*m_height - 1; // will actually jump (w*h-1)*sizeof(T) bytes
	}
}
//-----------------------------------------------------------------------------------------

template<class T>
int CMatrix<T>::GetWidth() const
{
	return m_width;
}
//-----------------------------------------------------------------------------------------

template<class T>
int CMatrix<T>::GetHeight() const
{
	return m_height;
}
//-----------------------------------------------------------------------------------------

template<class T>
int CMatrix<T>::GetNRows() const
{
	return m_height;
}
//-----------------------------------------------------------------------------------------

template<class T>
int CMatrix<T>::GetNCols() const
{
	return m_width;
}
//-----------------------------------------------------------------------------------------
// Elementwise multiplication
template<class T>
CMatrix<T> CMatrix<T>::Mul(const CMatrix<T> & a)
{
	ASSERT((a.GetNRows()==GetNRows())&&(a.GetNCols()==GetNCols()));
	CMatrix<T> prod(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pProd=prod.GetFirstPtr();
	T * pA=a.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pProd++=(*pMe++) * (*pA++);
	}
	return prod;
}

//-----------------------------------------------------------------------------------------
// Elementwise division
template<class T>
CMatrix<T> CMatrix<T>::Div(const CMatrix<T> & a)
{
	ASSERT((a.GetNRows()==GetNRows())&&(a.GetNCols()==GetNCols()));
	CMatrix<T> prod(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pProd=prod.GetFirstPtr();
	T * pA=a.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pProd++=(*pMe++) / (*pA++);
	}
	return prod;

}

//-----------------------------------------------------------------------------------------

template<class T>
T& CMatrix<T>::operator()(const int row, const int column)
{
	ASSERT(row < m_height && row >= 0 && column < m_width && column >= 0);
	return *(m_pElems + (row * m_width) + column);
}

//-----------------------------------------------------------------------------------------

template<class T>
T CMatrix<T>::operator()(const int row, const int column) const
{
	ASSERT(row < m_height && row >= 0 && column < m_width && column >= 0);
	return *(m_pElems + (row * m_width) + column);
}

//-----------------------------------------------------------------------------------------

template<class T>
T& CMatrix<T>::operator()(const int index)
{
	ASSERT(index < m_height*m_width && index>=0);
	return *(m_pElems + index);
}

//-----------------------------------------------------------------------------------------

template<class T>
T CMatrix<T>::operator()(const int index) const
{
	ASSERT(index < m_height*m_width && index>=0);
	return *(m_pElems + index);
}

//-----------------------------------------------------------------------------------------

template <class T>
CMatrix<T> CMatrix<T>::operator+(const CMatrix<T> & a)
{
	ASSERT((a.GetNRows()==GetNRows())&&(a.GetNCols()==GetNCols()));
	CMatrix<T> sum(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pSum=sum.GetFirstPtr();
	T * pA=a.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pSum++=(*pMe++) + (*pA++);
	}
	return sum;
} // operator+

//-----------------------------------------------------------------------------------------

template <class T>
CMatrix<T> CMatrix<T>::operator+(const T & t)
{
	CMatrix<T> sum(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pSum=sum.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pSum++=(*pMe++) + t;
	}
	return sum;
} // operator+

//-----------------------------------------------------------------------------------------

template <class T>
CMatrix<T> CMatrix<T>::operator-(const CMatrix<T> & a)
{
	ASSERT((a.GetNRows()==GetNRows())&&(a.GetNCols()==GetNCols()));
	CMatrix<T> sum(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pSum=sum.GetFirstPtr();
	T * pA=a.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pSum++=(*pMe++) - (*pA++);
	}
	return sum;
} // operator-

//-----------------------------------------------------------------------------------------

template <class T>
CMatrix<T> CMatrix<T>::operator-(const T & t)
{
	CMatrix<T> sum(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pSum=sum.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pSum++=(*pMe++) - t;
	}
	return sum;
} // operator-

//-----------------------------------------------------------------------------------------

template <class T>
CMatrix<T> CMatrix<T>::operator-(void)
{
	CMatrix<T> neg(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pNeg=neg.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pNeg++=0-(*pMe++);
	}
	return neg;
} // operator-

//-----------------------------------------------------------------------------------------

template <class T>
CMatrix<T> CMatrix<T>::operator*(const CMatrix<T> & b)
{
	ASSERT(GetNCols()==b.GetNRows());
	CMatrix<T> prod(GetNRows(),b.GetNCols());
	for(register int i=0; i<GetNRows(); i++)
		for(register int j=0; j<b.GetNCols(); j++){
			prod(i,j)=T(0);
			for(register int k=0; k<GetNCols(); k++)
				prod(i,j)+=(*this)(i,k)*b(k,j);
		} // for each elt, calc sum of products
	return prod;
} // operator*


//-----------------------------------------------------------------------------------------

template <class T>
CMatrix<T> CMatrix<T>::operator*(const T & t)
{
	CMatrix<T> prod(GetNRows(), GetNCols());
	T * pMe, * pLast;
	GetFirstLastPtrs(pMe,pLast);
	T * pProd=prod.GetFirstPtr();
	while(pMe<=pLast)
	{
		*pProd++=(*pMe++) * t;
	}
	return prod;
} // operator*

//-----------------------------------------------------------------------------------------

/*
template <class T>
bool CMatrix<T>::SaveMatlab(const char *fname, const char *vname )const // filename, variablename
{				
	FILE *fp = fopen(fname, "wb");
	if(!fp)
		return false;
	else
	{
		WriteMatlabVariable(fp, vname);
		fclose(fp);
		return true;
	} // if-else
} // SaveMatlab()
*/

//-----------------------------------------------------------------------------------------

/*
template <class T>
void CMatrix<T>::WriteMatlabVariable(FILE *fp, const char *pname) const
{
	Fmatrix x;
	x.type = 0; // PC data
	x.mrows = this->GetNRows();
	x.ncols = this->GetNCols();
	x.imagf = 0;
	x.namlen = strlen(pname)+1;
	const int mn = x.mrows*x.ncols;
	fwrite(&x, sizeof(Fmatrix), 1, fp);
	fwrite(pname, sizeof(char), (int)x.namlen, fp);
	// Always write doubles
	for(register int c = 0; c < x.ncols; c++)
	{
		for(register int r = 0; r < x.mrows; r++)
		{
			const double d = double((*this)(r, c));
			fwrite(&d, sizeof(double), 1, fp);
		} // for r
	}
} // WriteMatlabVariable()
*/

//-----------------------------------------------------------------------------------------





///////////////////////////////////////////////////////////////////////////////////
//
// Arithmetic operators - some global (friend) operators
//
//

//-----------------------------------------------------------------------------------------
	
template <class T> 
CMatrix<T> operator+(const T & f,const CMatrix<T> & M)
{
	T *pM=M.GetFirstPtr();
	ASSERT(pM);
	CMatrix<T> sum(M.GetNRows(),M.GetNCols());
	T *pSum,*pLast;
	sum.GetFirstLastPtrs(pSum,pLast);
	while(pSum<=pLast)
		*pSum++ = f+ (*pM++);
	return sum;
} // operator+

//-----------------------------------------------------------------------------------------

template <class T> 
CMatrix<T> operator-(const T & f,const CMatrix<T> & M)
{
	T *pM=M.GetFirstPtr();
	ASSERT(pM);
	CMatrix<T> sum(M.GetNRows(),M.GetNCols());
	T *pSum,*pLast;
	sum.GetFirstLastPtrs(pSum,pLast);
	while(pSum<=pLast)
		*pSum++=f-(*pM++);
	return sum;
} // operator-

//-----------------------------------------------------------------------------------------

template <class T> 
CMatrix<T> operator*(const T & f,const CMatrix<T> & M)
{
	T *pM=M.GetFirstPtr();
	ASSERT(pM);
	CMatrix<T> prod(M.GetNRows(),M.GetNCols());
	T * pProd, * pLast;
	prod.GetFirstLastPtrs(pProd,pLast);
	while(pProd<=pLast)
		*pProd++=f*(*pM++);
	return prod;
} // operator*

//-----------------------------------------------------------------------------------------

template <class T>
void CMatrix<T>::Fill(T data)
{
	for (int i=m_width*m_height;--i>=0;)
		m_pElems[i]= data;
}

//-----------------------------------------------------------------------------------------
