/***********************************************************************************
 *
 *   filesystem.h -- Interface to file system
 *
 *   Erik Persson and Björn Nilsson, 2000-
 */

#ifndef FILESYSTEM_H
#define FILESYSTEM_H

#include "types/vector.h"
#include "types/string.h"
#include <time.h>
#ifndef _WIN32
#include <utime.h>
#endif

#ifdef MFC_CLASSES
class CBitmap;
#endif

#ifdef _WIN32
// for FILETIME structure
#include <olectl.h>
#endif

#ifdef _WIN32
#define FILENAME_SLASH '\\'
#else
#define FILENAME_SLASH '/'
#endif

class CLoadBuf;

class CFileSystem
{
	CVector<CLoadBuf *> m_vpLoadBuf;
	mutable CString m_sLastCustomError;
public:
	int m_LoadCnt;

public:
	CFileSystem();
	~CFileSystem();

	// File-system error status
	CString GetLastError() const;
	void SetLastError(const char *ach);

	// Current Directory
	CString GetCurrentDir() const;
#ifdef _WIN32
	void SetCurrentDir(const char *Dir);
#endif

	// Current executable
protected:
	CString m_sExeName;
public:
	void SetExecutableFileName(const char *ach);
	CString GetExecutableFileName() const;

	// User home directory
	CString GetUserDir() const;


	// Simple file load/store
public:
	bool CreateBinary(const char *aFileName, const char *aBinary, size_t nBinary);
	bool CreateBinary(const char *aFileName, const char *aStr);
	CLoadBuf *Load(const char *aFileName);
	CLoadBuf *AllocBuf();
	void FreeBuf(CLoadBuf *pLoadBuf);
	char *GetFileData(CLoadBuf *pLoadBuf);
	bool IsEqual(const char *aFileName, const char *aData);
	bool IsEqualFilename(const char *aFilename0, const char *aFilename1);

	bool FileExists(const char *ach) const;
#ifdef _WIN32
	bool IsProtected(const char *ach) const;
#endif
	int GetAttributes(const char *aFilename) const;

	/*
#ifdef _WIN32
	CTime GetModifiedTime(const char *aFileName) const;
#else
#error TODO: Implement CFileSystem::GetModifiedTime under Unix
#endif
	*/

	/*
	CVector<CString> ListFiles(CString sDir, CString sPattern,
		bool bDirs, bool bRecursive);
	*/

	// Directory Creation
#ifdef _WIN32
	bool CreateDirRecursive(const char *ach);
#endif

	// File name handling
public:
	bool IsRoot(const char *ach) const;
	// bool IsSimpleFileName(const char *ach) const;
	bool IsAbsPath(const char *ach) const;
	CString SubtractFullNames(const char *aFullFileName, const char *aFullDirName) const;
	CString AddFileNames(const char *aFullName, const char *aRelName) const;
public:
	CString GetFullName(const char *ach) const;
	CString GetPathPrefix(const char *ach) const;
	CString GetRelativeName(const char *aFileName, const char *aDirName) const;
	CString GetCombinedName(const char *aDirName, const char *aRelName) const;
	CString GetNameWithoutExtension(const char *ach) const;
	CString GetFileExtension(const char *ach) const;
	CString GetStrippedName(const char *ach) const;
	CString GetParentName(const char *ach) const;
#ifdef MFC_CLASSES
	bool SaveBitmap(HDC hDC, CBitmap &bm, CString sFileName);
#endif
	CString ToSystemSlash(const char *ach) const;
	bool ListFiles(const char *aDir, const char *aFilter, CVector<CString> &v_out);
	bool ListSubdirectories(const char *aDir, const char *aFilter, CVector<CString> &v_out);
	bool MatchFilter(const CString &s0, const CString &sF, bool bCaseDep= false);

	CString GetSafeStrippedName(CString s);
	CString GetTmpFilename(const char *aDir) const;

#ifdef _WIN32
	// WINDOWS
	bool Delete(const char *aFilename);
	bool GetFileTime(const char *aFilename, FILETIME *lpCreate, FILETIME *lpAccess, FILETIME *lpWrite);
	bool SetFileTime(const char *aFilename, FILETIME *lpCreate, FILETIME *lpAccess, FILETIME *lpWrite);
	bool HasFileTime(const CString &sFilename, const FILETIME *lpCreate, const FILETIME *lpAccess, const FILETIME *lpWrite);
	bool GetFileSize(const char *aFilename, long long *p_sz);
	bool CopyFile(const CString &sSrc, const CString &sDst, CString *psError);
	int	CompareFiles(const char *aFile0, const char *aFile1, CString *psError);
#else
	// POSIX
	bool GetFileTime(const char *aFilename, time_t &atime, time_t &mtime);	
	bool SetFileTime(const char *aFilename, const time_t atime, const time_t mtime);
	bool HasFileTime(const char *aFilename, const time_t *p_atime, const time_t *p_mtime);
#endif
};


/***********************************************************************************
 *
 *  CLoadBuf class: realloc-based file read buffer
 *
 */

class CLoadBuf
{
protected:
	bool QuickLoad(const char *aFileName);
public:
	char *m_aLoadBuf;
	size_t m_BufSize; // No. of bytes allocated for buffer ( != bytes loaded)
	size_t m_nFileSize; // No. of bytes loaded
	bool m_bUsed;

	CLoadBuf();
	~CLoadBuf();

	bool SetSize(size_t size);
	bool Load(const char *aFileName);

};

extern CFileSystem g_FileSystem;

#endif
