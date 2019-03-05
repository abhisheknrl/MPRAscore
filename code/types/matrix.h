/*******************************************************************************
 *
 *   CMatrix.h -- Matrix class
 *
 */

#ifndef CMATRIX_H
#define CMATRIX_H

#include <stdio.h>
#include "types/assert.h"

template <class T> 
class CMatrix;

// Friends 
template <class T> CMatrix<T> operator+(const T & f,const CMatrix<T> & M);
template <class T> CMatrix<T> operator-(const T & f,const CMatrix<T> & M);
template <class T> CMatrix<T> operator*(const T & f,const CMatrix<T> & M);

template<class T>
class  CMatrix
{
	/*
	// These generate mysterious warnings with GCC
	friend CMatrix<T> operator+(const T & f,const CMatrix<T> & M);
	friend CMatrix<T> operator-(const T & f,const CMatrix<T> & M);
	friend CMatrix<T> operator*(const T & f,const CMatrix<T> & M);
	*/

private:
	T		*m_pElems;
	int		m_width;		
	int		m_height;
public:
	// default ctor
	CMatrix();

	// init ctor
	CMatrix(const int nRows, const int nCols);

	// init ctor
	CMatrix(const int nRows, const int nCols, const T initValue);

	// 
	//CMatrix: Init ctor. This version of the ctor takes a string of values to initialize the CMatrix with
	//PARAM: width, height, and a string of numbers
	//PRECOND: Only numbers are valid in this ctor. There must exactly the same amount of numbers as the total 
	// number of elements in the CMatrix.
	//POSTCOND: The CMatrix will be initialized to the specified size and set with the numbers in the string
	//NOTE: The function takes no responsibility for what will happen if the string has more or less numbers
	// than the number of elements.
	CMatrix(const int nRows, const int nCols, const char * pszInitValues);

	// copy ctor
	CMatrix(const CMatrix<T> & other);

	//operator =: Assignment operator overload
	//PARAM: The other CMatrix to copy
	//PRECOND:	The other CMatrix must be initialized.
	//POSTCOND:	This CMatrix will be an exact copy of the the other CMatrix with all its elements.
	//NOTE:
	CMatrix<T>& operator=(const CMatrix<T> & other);

	// dtor
	~CMatrix();

	//ReInit: ReInitializes the CMatrix to a new size, old values are destroyed if they exist
	//PARAM: The new dimensions of the CMatrix (rows/columns) = (height/width);
	//PRECOND: None
	//POSTCOND: The CMatrix will have the new dimensions
	//NOTE:
	void ReInit(const int NRows, const int NCols);

	//ReInit: ReInitializes the CMatrix to a new size and sets the initvalue, old values are destroyed if they exist
	//PARAM: The new dimensions of the CMatrix (rows/columns) = (height/width);
	//PRECOND: None
	//POSTCOND: The CMatrix will have the new dimensions
	//NOTE:
	void ReInit(const int NRows, const int NCols, const T initValue);

	//InitValues: Sets all the elementvalues to value.
	//PARAM: The value to set.
	//PRECOND: The CMatrix exists (i.e. it has dimensions)
	//POSTCOND: All elements values will be set to value
	//NOTE:
	void InitValues(const T value);

	//GetFirstPtr: Returns a pointer to the first element in the matrix.
	//PARAM: None
	//PRECOND:	None.
	//POSTCOND: The pointer to the first element in the data is returned, null if the CMatrix is uninitialized.
	//NOTE:
	T * GetFirstPtr() const;

	//GetFirstLastPtrs: Retrieves pointers to the first element (top-left) and the last element (bottom-right)
	// in the CMatrix for fast iteration.
	//PARAM:	Pointers that should be set sent in by reference.
	//PRECOND:	None.
	//POSTCOND:	If the CMatrix is uninitialized both pointers will be set to NULL, otherwise pFirst
	// will point to the first element, and pLast to the last element (please note: NOT BEYOND the last element)
	//NOTE:
	void GetFirstLastPtrs(T ** pFirst, T ** pLast) const;

	// GetPointer: Retrieves the address of a specific element.
	T *GetPointer(int nRow, int nCol) const { ASSERT(nRow>=0 && nRow<m_height && nCol>=0 && nCol<m_width); return m_pElems + nRow*m_width + nCol; }

	//GetWidth: Gets the width == number of columns of the matrix
	int GetWidth() const;
	int NCols() const { return GetWidth(); }; // For backward compatibility.

	//GetHeight: Gets the height == number of rows of the matrix
	int GetHeight() const;
	int NRows() const { return GetHeight(); } // For backward compatibility;

	bool IsEmpty() const { return GetWidth()==0 || GetHeight()==0; }

	//GetNRows: Gets the number of rows == height of the matrix
	int GetNRows() const;

	//GetNCols: Gets the number of columns == width of the matrix
	int GetNCols() const;

	//operator(): Parenthesis operator for referencing an element in the array.
	// This version can be used to write to the element.
	//PARAM: Row and column of the desired element.
	//PRECOND: 0 <= row < NRows && 0 <= column < NCols
	//POSTCOND: The value of the element will be returned by reference.
	T& operator()(const int row, const int column);

	//operator(): Parenthesis operator for referencing an element in the array.
	// This version can only be used to read the element.
	//PARAM: Row and column of the desired element.
	//PRECOND: 0 <= row < NRows && 0 <= column < NCols
	//POSTCOND: The value of the element will be returned by value.
	T operator()(const int row, const int column) const;

	//operator(): Parenthesis operator for referencing an element in the array.
	// This version can be used to write to the element.
	//PARAM: index (row*width+col) of the desired element.
	//PRECOND: 0 <= row < NRows*NCols
	//POSTCOND: The value of the element will be returned by reference.
	T& operator()(const int index);

	//operator(): Parenthesis operator for referencing an element in the array.
	// This version can only be used to read the element.
	//PARAM: index (row*width+col) of the desired element.
	//PRECOND: 0 <= index < NRows*NCols
	//POSTCOND: The value of the element will be returned by value.
	T operator()(const int index) const;

	T GetAt(int nRow, int nCol) const 
	{	
		ASSERT((nRow<m_height) && (nCol<m_width));
		return *(m_pElems + nRow*m_width + nCol); 
	}
	void SetAt(int nRow, int nCol, T data) 
	{	
		ASSERT((nRow<m_height) && (nCol<m_width));
		*(m_pElems + nRow*m_width + nCol)= data; 
	}

	///////////////////////////////////////////////////////////////////////////////////
	// Elemetwise arithmetic functions
	
	// Elementwise multiplication
	CMatrix<T> Mul(const CMatrix<T> & a);
	
	// Elementwise division
	CMatrix<T> Div(const CMatrix<T> & a);


	///////////////////////////////////////////////////////////////////////////////////
	// Arithmetic operators

	CMatrix<T>	operator+ (const CMatrix<T> &);
	CMatrix<T>	operator+ (const T &);

	CMatrix<T>	operator- (const CMatrix<T> &);
	CMatrix<T>	operator- (const T &);
	CMatrix<T>	operator- (void);
	
	// CMatrix multiplication
	CMatrix<T>	operator* (const CMatrix<T> &);
	CMatrix<T>	operator* (const T &);

	typedef struct // Used by SaveMatlab()
	{
		long type;		// type 0=double, 50=compact???
		long mrows;		// row dimension
		long ncols;		// column dimension
		long imagf;		// flag indicating imag part
		long namlen;	// name length (including NULL)
	} FCMatrix;

	// bool SaveMatlab(const char*, const char*) const; // filename, variablename.

	void Fill(T data);
private:
	// void WriteMatlabVariable(FILE *fp, const char *pname) const;
};

#include "matrix.inl"

#endif

