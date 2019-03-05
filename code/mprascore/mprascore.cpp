/*******************************************************************************
 *
 *   mprascore.exe - Calculate differential regulation scores for MPRA data
 *
 *   Björn Nilsson and Abhishek Niroula, 2018-
 *
 */

#include "types/tablefun.h"
#include "system/consoleapp.h"
#include "system/filesystem.h"
#include "types/stringfun.h"
#include "math/statfun.h"
#include "math/variance.h"
#include <float.h>
#include "types/minmax.h"
#include "math/rand.h"

#define SCORE_LOG2_RNA_DNA (0)
#define SCORE_RNA_DNA (1)
#define SCORE_RNA_ONLY (2)
#define SCORE_MAX (2)

int g_nRefs= 0;

/******************************************************************************/

class CApp : public CConsoleApp
{
protected:
	bool m_bHelpShown;

	CTableIndex m_idx_variants;
	CTableIndex m_idx_barcodes;
	CMatrix<double> m_idx_barcodes_RNA;
	CMatrix<double> m_idx_barcodes_DNA;

	CVector<int> m_barcodes_vFirst;
	CVector<int> m_barcodes_vLast;

	CString m_sRNAcols, m_sDNAcols;
	CVector<int> m_vRNAcols;
	CVector<int> m_vDNAcols;

	bool m_bExportTTest;

	// bool m_bSmythRegularization;
	CVector<double> m_vSmyth_DF0;
	CVector<double> m_vSmyth_Var0;

	CVector<Table_SortItem_int> m_vSI_barcodes;

	int m_nScoreType;
	int m_nMinBarcodesPerOligo;
	bool m_bIncludePosStrand, m_bIncludeNegStrand;
	int m_nNormalizationScaleFactor;
	int m_nNullDistribution_permutations;
	bool m_bExportHeader;
	bool m_bExportNull;

	void DisplayHelp();
	void SetDefaultParameters();
	bool CheckParameters();
	bool CreateOutputFile(const char *aFilename_out);
	bool CreateOutputFile_SNP(FILE *pf, const int r0, const int r1);
	bool CreateBarcodeLookup();
	bool GetBarcodeStats();
	bool NormalizeBarcodeCounts();

	double GetBarcodeScore(const double rna, const double dna);
	bool GetBarcodeScores(const int b0, const int b1, CMatrix<double> &m_scores);
	void GetOligoScores(CMatrix<double> &mBarcodeScores_a, CMatrix<double> &mBarcodeScores_r, CMatrix<double> &mB, CMatrix<double> &mW);

	bool GetStatistics_SNP(Table_SortItem_int *pSI, const int i0, const int i1, const int c_RSID, const int c_Oligo,
		const int c_Ref, const int c_Alt, const int c_Allele_in_oligo, const int c_Strand, const int c_Offset,
		CVector<int> &vOut_r_variant,
		CVector<double> &vOut_beta,
		CVector<double> &vOut_beta_var,
		CVector<double> &vOut_beta_n);
	bool GetStatistics_SNP2(Table_SortItem_int *pSI, const int i0, const int i1, const int c_RSID, const int c_Oligo,
		const int c_Ref, const int c_Alt, const int c_Allele_in_oligo, const int c_Strand, const int c_Offset,
		CVector<int> &vOut_r_variant,
		CVector<double> &vOut_fem,
		CVector<double> &vOut_fem_var,
		CVector<double> &vOut_fem_n,
		CVector<double> &vOut_fem_z,
		CVector<double> &vOut_ttest_t,
		CVector<double> &vOut_ttest_df,
		double *p_NULL, int &n_NULL, const int n_NULL_max);

	bool ExportNullDistribution(const char *aFilename, const double *p_NULL, const int n_NULL);

	virtual bool OnSwitch(const char *ach, CString &sSwitch);

public:
	virtual int Main(int argc, TCHAR* argv[], TCHAR* envp[]);
};

//
// CCommand overrides

void CApp::DisplayHelp()
{
	m_bHelpShown=true;

	CString sName= g_FileSystem.GetExecutableFileName();
	sName= g_FileSystem.GetNameWithoutExtension(g_FileSystem.GetStrippedName(sName));
	printf("Usage: %s <oligo file> <barcode count file> <output file>\n\n", (const char *)sName);

	printf("Input files:\n\n");
	printf("Both input files must be tab-delimited with one header row.\n\n");

	printf("The oligo file must contain the following columns:\n\n");
	printf("Oligo\t\tIdentifier of enhancer oligonucleotide sequence.\n");
	printf("RSID\t\tIdentifier of the variant the oligo represents.\n");
	printf("Ref\t\tVariant reference allele.\n");
	printf("Alt\t\tVariant alternative allele.\n");
	printf("Allele_in_oligo\tAllele in oligo (must be the Ref or Alt allele).\n");
	printf("Strand\t\tOligo strand. O (must be + or -).\n");
	printf("baseOffset\tLocation of variant, in number of base pairs from\n\t\tthe sequence center (0 if variant is in the center of the\n\t\tsequence).\n\n");

	printf("Barcode count file must contain a column named Oligo with identifiers of\nenhancer oligonucleotide sequences matching those in the oligo file.\n");
	printf("In addition, it must contain columns with RNA and DNA counts. The names of\nthese columns are specified using -RNAcols and -DNAcols, see below.\n\n");

	printf("Options:\n\n");

	printf("-RNAcols:col,col...\tSpecify columns containing RNA counts.\n");
	printf("-DNAcols:col,col...\tSpecify columns containing DNA counts.\n");
	/*
	printf("-score\t\t\tSpecify barcode score type:\n");
	printf("\t\t\t0\tlog2((1+#RNA)/(1+#DNA)) (default)\n");
	printf("\t\t\t1\t#RNA/(1+#DNA)\n");
	printf("\t\t\t2\t#RNA\n");
	*/
	printf("-minbarcodes:nnn\tSpecify no. of barcodes required to include\n\t\t\tan oligo (default %d).\n", m_nMinBarcodesPerOligo);
	printf("-posonly\t\tLimit to oligos representing the positive strand.\n");
	printf("-negonly\t\tLimit to oligos representing the negative strand.\n");
	printf("-p:nnn\t\t\tSpecify no. of permutations used to calculate null\n\t\t\tdistribution (default %d).\n", m_nNullDistribution_permutations);
	printf("-?\t\t\tDisplay help text.\n\n");
}

bool CApp::OnSwitch(const char *ach, CString &sSwitch)
{
	if (sSwitch=="?")
	{
		DisplayHelp();
		return true;
	}
	else if (sSwitch=="RNAcols")
	{
		m_sRNAcols= ach;
		return true;
	}
	else if (sSwitch=="DNAcols")
	{
		m_sDNAcols= ach;
		return true;
	}
	else if (sSwitch=="score")
	{
		if (Table_CheckIntFormat_complete(ach)!=0)
		{
			m_nScoreType= atof(ach);
			return true;
		}		
	}
	else if (sSwitch=="negonly")
	{
		m_bIncludePosStrand= false;
		m_bIncludeNegStrand= true;
		return true;
	}
	else if (sSwitch=="posonly")
	{
		m_bIncludePosStrand= true;
		m_bIncludeNegStrand= false;
		return true;
	}
	else if (sSwitch=="minbarcodes")
	{
		if (Table_CheckIntFormat_complete(ach)!=0)
		{
			m_nMinBarcodesPerOligo= atoi(ach);
			return true;
		}
	}
	else if (sSwitch=="p")
	{
		if (Table_CheckIntFormat_complete(ach)!=0)
		{
			m_nNullDistribution_permutations= atoi(ach); 
			return true;
		}
	}
	else if (sSwitch=="noheader")
	{
		m_bExportHeader= false;
		return true;
	}
	else if (sSwitch=="exportnull")
	{
		m_bExportNull= true;
		return true;
	}
	else if (sSwitch=="exportttest")
	{
		m_bExportTTest= true;
		return true;
	}

	CString s= sSwitch;
	if (ach && *ach)
		s += ":" + CString(ach);
	ReportError("Illegal switch", (const char *)s);
	return false;
}

void CApp::SetDefaultParameters()
{
	m_bHelpShown= false;
	m_bExportHeader= true;
	m_bIncludePosStrand= true;
	m_bIncludeNegStrand= true;
	m_nMinBarcodesPerOligo= 3;
	m_nNullDistribution_permutations= 100;
	m_nNormalizationScaleFactor= 10000000; 
	m_nScoreType= SCORE_LOG2_RNA_DNA;
	m_bExportTTest= false;
}

bool CApp::CheckParameters()
{
	if (m_bHelpShown)
		return false;

	if (m_nScoreType<0)
	{
		ReportError("Score type not set, use -score switch", 0);
		return false;
	}
	if (m_nScoreType>SCORE_MAX)
	{
		ReportError("Illegal score type, check -score switch", 0);
		return false;
	}

	// resolve RNA cols
	if (m_sRNAcols.IsEmpty())
	{
		ReportError("-RNAcols not set", 0);
		return false;
	}
	printf("Resolving RNA columns...\n");
	
	CVector<CString> vsRNAcols;
	SplitAtChar(m_sRNAcols, ',', vsRNAcols);
	for (int i=0;i<vsRNAcols.GetSize();i++)
	{
		int c= Table_FindCol(m_idx_barcodes, vsRNAcols[i]);
		if (c<0)
		{
			ReportError(::Format("Column '%s' not found", (const char *)vsRNAcols[i]), 0);
			return false;
		}
		if (GetVerbosity()>2)
			printf("%s\t%d\n", (const char *)vsRNAcols[i], c);
		m_vRNAcols.Add(c);
	}

	// resolve DNA cols
	if (m_sDNAcols.IsEmpty())
	{
		ReportError("-DNAcols not set", 0);
		return false;
	}
	printf("Resolving DNA columns...\n");
	
	CVector<CString> vsDNAcols;
	SplitAtChar(m_sDNAcols, ',', vsDNAcols);
	for (int i=0;i<vsDNAcols.GetSize();i++)
	{
		int c= Table_FindCol(m_idx_barcodes, vsDNAcols[i]);
		if (c<0)
		{
			ReportError(::Format("Column '%s' not found", (const char *)vsDNAcols[i]), 0);
			return false;
		}
		if (GetVerbosity()>2)
			printf("%s\t%d\n", (const char *)vsDNAcols[i], c);
		m_vDNAcols.Add(c);
	}

	for (int i=0;i<m_vRNAcols.GetSize();i++)
		for (int j=0;j<m_vDNAcols.GetSize();j++)
			if (m_vRNAcols[i]==m_vDNAcols[j])
			{
				ReportError("RNA and DNA columns overlap, check -RNAcols and -DNAcols switches", 0);
				return false;
			}			

	printf("Parsing barcode counts...\n");

	// import RNA and DNA matrices
	m_idx_barcodes_RNA.ReInit(m_idx_barcodes.GetRowCount()-1, m_vRNAcols.GetSize());
	m_idx_barcodes_DNA.ReInit(m_idx_barcodes.GetRowCount()-1, m_vDNAcols.GetSize());
	for (int r=1;r<m_idx_barcodes.GetRowCount();r++)
	{
		for (int j=0;j<m_vRNAcols.GetSize();j++)
		{
			if (Table_CheckFloatFormat(m_idx_barcodes.GetFirstPointer(r, m_vRNAcols[j]))!=0)
				m_idx_barcodes_RNA.SetAt(r-1, j, atof(m_idx_barcodes.GetFirstPointer(r, m_vRNAcols[j])));
			else
			{
				ReportError("Illegal floating point format", Table_GetAt(m_idx_barcodes, r, m_vRNAcols[j]));
				return false;
			}
		}
		for (int j=0;j<m_vDNAcols.GetSize();j++)
		{
			if (Table_CheckFloatFormat(m_idx_barcodes.GetFirstPointer(r, m_vDNAcols[j]))!=0)
				m_idx_barcodes_DNA.SetAt(r-1, j, atof(m_idx_barcodes.GetFirstPointer(r, m_vDNAcols[j])));
			else
			{
				ReportError("Illegal floating point format", Table_GetAt(m_idx_barcodes, r, m_vDNAcols[j]));
				return false;
			}
		}
	}

	// sanity check no of permut
	if (m_nNullDistribution_permutations<10 || m_nNullDistribution_permutations>1000000000)
	{
		ReportError("Illegal number of permutations, check -p switch", 0);
		return false;
	}

	return true;
}

//
// helper functions

void PQfromZ(const CVector<double> &vZ, double *pZ0, const int &nZ0, CVector<double> &vP, CVector<double> &vQ)
{
	// sort z scores
	const int nZ= vZ.GetSize();
	CVector<Stat_SortItem_double> vSI(nZ);
	Stat_SortItem_double *pSI= vSI.GetBuffer();
	for (int i=0;i<nZ;i++)
	{
		Stat_SortItem_double si;
		si.m_Value= __abs(vZ[i]);
		si.m_nIndex= i;
		vSI.SetAt(i, si);
	}
	Stat_QuicksortIndexed_double(pSI, nZ, true);

	// sort null distribution of z scores
	for (int i=0;i<nZ0;i++)
		pZ0[i]= __abs(pZ0[i]);
	Stat_Quicksort_double(pZ0, nZ0, true);

	// calculate false discovery rate (q-value) for each z
	vQ.SetSize(nZ);
	vP.SetSize(nZ);
	double pi0= 1.0; // no pi0 estimation implemented at the moment (but close to 1.0 for most applications)
	double qmin= pi0; 
	int i0= 0;
	for (int i=0;i<nZ;i++)
	{
		while (i0<nZ0 && pZ0[i0]<=pSI[i].m_Value)
			i0++;

		int j= pSI[i].m_nIndex;
		double A1= double(nZ-i)/nZ;
		double A0= double(nZ0-i0)/nZ0;
		double q= pi0*A0/A1;
		if (q<qmin)
			qmin= q;
		vQ.SetAt(j, qmin);

		double p= double(nZ0-i0+1)/(nZ0+1);
		vP.SetAt(j, p);
	}
}

CVector<double> QfromP(CVector<double> vP, bool bSortQ, bool bEstimatePi0)
{
	// Sort P's
	const int N= vP.GetSize();
	CVector<Stat_SortItem_double> vSI;
	vSI.SetSize(N);
	for (int i=0;i<N;i++)
	{
		Stat_SortItem_double si;
		si.m_Value= vP[i];
		si.m_nIndex= i;
		vSI.SetAt(i, si);
	}
	Stat_SortItem_double *pSI= vSI.GetBuffer(N);
	Stat_QuicksortIndexed_double(pSI, N, true);

	/*
	// Override Benjamini-Yuketeli correction, i.e. assume independent or 
	// positively dependent p's which yields the standard, and slightly 
	// less convervative Benjamini-Hochberg correction.

	double benj_yuke= 0.0;
	for (i=0;i<N;i++)
		benj_yuke += 1.0/(i+1);
	*/

	// Estimate pi0
	double pi0= 1.0;
	if (bEstimatePi0 && N>1000)
	{
		printf("Estimating pi0\n");
		int nHigh= 0;
		double fLambda= 0.6f; // Cutoff
		for (int i=0;i<N;i++)
		{
			if (vP[i]>fLambda)
				nHigh++;
		}
		pi0= (double(nHigh)/(1.0f-fLambda))/N;
		if (pi0>1.0)
			pi0= 1.0;
	}
	//else
	//	printf("Skipping pi0 estimation\n");

	//if (GetVerbosity()>1)
	//	printf("Type of multiple-testing correction: False discovery rate (pi0=%.2f).\n", pi0);

	CVector<double> vQ;
	vQ.SetSize(N);
	double fmin= pi0;
	for (int i=N;--i>=0;)
	{
		int j= vSI[i].m_nIndex;
		double f= pi0*vSI[i].m_Value/(double(i+1)/N);

		// Constrained
		if (f<fmin)
			fmin= f;

		if (bSortQ)
			vQ.SetAt(i, fmin);
		else
			vQ.SetAt(j, fmin);
	}

	return vQ;
}

int FindInIndex_alphabetic(const CVector<Table_SortItem_int> &vSI, CString &sValue)
{
	// Assumes vSI is sorted alphabetically in ascending order.
	// Returns position of sValue in vSI.
	int i0= 0;
	int i1= vSI.GetSize();
	while (i0<i1)
	{
		int i_mid= (i0+i1)/2;

		if (vSI[i_mid].m_sKey==sValue)
			return i_mid;

		if (strcmp(vSI[i_mid].m_sKey, sValue)<0)
			i0= i_mid+1;
		else
			i1= i_mid;		
	}

	return -1;
}

//
// Actual program

bool CApp::CreateBarcodeLookup()
{
	printf("Creating barcode lookup table...\n");

	int cOligo= Table_FindCol(m_idx_barcodes, "Oligo");
	if (cOligo<0)
	{
		ReportError("Column 'Oligo' not found", 0);
		return false;
	}

	const int R= m_idx_barcodes.GetRowCount()-1;
	// printf("R=%d\n", R);
	m_vSI_barcodes.SetSize(R);
	Table_SortItem_int *pSI= m_vSI_barcodes.GetBuffer();

	// sort on OligoID
	for (int r=0;r<R;r++)
	{
		Table_SortItem_int si;
		si.m_nIndex= r;
		si.m_sKey= Table_GetAt(m_idx_barcodes, r+1, cOligo);
		si.m_Value= 0;
		pSI[r]= si;
	}
	Table_SortItems_int(pSI, R, true, false);

	// create lists of positions of the first and last barcodes for each oligo within the index
	for (int i0=0;i0<R;)
	{
		int i1= i0;
		while (i1<R && pSI[i0].m_sKey==pSI[i1].m_sKey)
			i1++;

		m_barcodes_vFirst.Add(i0); // position of first and last barcode for this oligo in m_vSI_barcodes 
		m_barcodes_vLast.Add(i1);
		i0= i1;
	}

	return true;
}

bool CApp::NormalizeBarcodeCounts()
{
	printf("Normalizing barcode counts...\n");

	const int scale= m_nNormalizationScaleFactor;

	for (int c=0;c<m_idx_barcodes_RNA.GetWidth();c++)
	{
		double s=0;
		for (int r=0;r<m_idx_barcodes_RNA.GetHeight();r++)
		{
			double rna= m_idx_barcodes_RNA.GetAt(r,c);
			if (rna<0)
			{
				ReportError("RNA counts cannot be negative", 0);
				return false;
			}
			s += rna;
		}
		if (s==0)
		{
			ReportError("Total RNA counts cannot be zero", 0);
			return false;
		}
		for (int r=0;r<m_idx_barcodes_RNA.GetHeight();r++)
			m_idx_barcodes_RNA.SetAt(r,c, m_idx_barcodes_RNA.GetAt(r,c)/s*scale);
	}

	for (int c=0;c<m_idx_barcodes_DNA.GetWidth();c++)
	{
		double s=0;
		for (int r=0;r<m_idx_barcodes_DNA.GetHeight();r++)
		{
			double dna= m_idx_barcodes_DNA.GetAt(r,c);
			if (dna<0)
			{
				ReportError("DNA counts cannot be negative", 0);
				return false;
			}
			s += dna;
		}	
		if (s==0)
		{
			ReportError("Total DNA counts cannot be zero", 0);
			return false;
		}

		for (int r=0;r<m_idx_barcodes_RNA.GetHeight();r++)
			m_idx_barcodes_DNA.SetAt(r,c, m_idx_barcodes_DNA.GetAt(r,c)/s*scale);
	}

	return true;
}

bool CApp::GetBarcodeStats()
{
	if (!NormalizeBarcodeCounts())
		return false;

	if (m_barcodes_vFirst.GetSize()!=m_barcodes_vLast.GetSize())
	{
		ReportError("Internal error: m_barcodes_vFirst.GetSize()!=m_barcodes_vLast.GetSize()", 0);
		return false;
	}

	printf("Calculating raw variances...\n");

	const int n_rep= m_idx_barcodes_RNA.GetWidth();
	m_vSmyth_DF0.SetSize(n_rep);
	m_vSmyth_Var0.SetSize(n_rep);

	// get score stats across each oligo
	CVector< CVector<double> > vVar(n_rep);
	CVector<double> *pVar= vVar.GetBuffer();
	CVector< CVector<double> > vDF(n_rep);
	CVector<double> *pDF= vDF.GetBuffer();

		CVector<double> v_s(n_rep);
		double *p_s= v_s.GetBuffer();
		CVector<double> v_ss(n_rep);
		double *p_ss= v_ss.GetBuffer();
		CVector<int> v_n(n_rep);
		int *p_n= v_n.GetBuffer();

	for (int i=0;i<m_barcodes_vFirst.GetSize();i++)
	{
		int b0= m_barcodes_vFirst[i];
		int b1= m_barcodes_vLast[i];

		for (int j=0;j<n_rep;j++)
		{
			p_s[j]= 0;
			p_ss[j]= 0;
			p_n[j]= 0;
		}

		for (int b=b0;b<b1;b++)
		{
			int r_barcode= m_vSI_barcodes[b].m_nIndex;
			/*
			if (r_barcode<0 || r_barcode>=m_idx_barcodes_DNA.GetHeight() || r_barcode>=m_idx_barcodes_RNA.GetHeight())
			{
				ReportError(::Format("r_barcode out of bounds: %d", r_barcode), 0);
				return false;
			}
			*/

			double mu_DNA= Stat_GetMean_double(m_idx_barcodes_DNA.GetPointer(r_barcode, 0), m_idx_barcodes_DNA.GetWidth());
			for (int rep=0;rep<n_rep;rep++)
			{
				double counts= m_idx_barcodes_RNA.GetAt(r_barcode, rep);
				double score_barcode= GetBarcodeScore(m_idx_barcodes_RNA.GetAt(r_barcode, rep), m_idx_barcodes_DNA.GetAt(r_barcode, rep));
				if (IsN_double(score_barcode))
				{
					p_s[rep] += score_barcode;
					p_ss[rep] += score_barcode*score_barcode;
					p_n[rep] ++;
				}
			}
		}

		for (int rep=0;rep<n_rep;rep++)	
		{
			if (p_n[rep]>2)
			{
				double var= Stat_GetVariance(p_s[rep], p_ss[rep], p_n[rep]);
				if (var>0)
				{
					pVar[rep].Add(var);
					pDF[rep].Add(p_n[rep]);
				}
			}
		}
	}

	printf("Calculating regularized variances...\n");
	bool bOk= true;
	for (int rep=0;bOk && rep<n_rep;rep++)	
	{
		CVector<double> v_out;
		double s0sqr, df0;
		bOk= GetVariances_Smyth(pVar[rep], pDF[rep], v_out, s0sqr, df0);

		if (GetVerbosity()>2)
		{
			printf("Smyth, d0[%d]= %g\n", rep, df0);
			printf("Smyth, var0[%d]= %g\n", rep, s0sqr);
		}
		m_vSmyth_DF0.SetAt(rep, df0);
		m_vSmyth_Var0.SetAt(rep, s0sqr);
	}

	return bOk;
}

double CApp::GetBarcodeScore(const double rna, const double dna)
{
	switch (m_nScoreType)
	{
		case SCORE_RNA_ONLY: return rna;
		case SCORE_LOG2_RNA_DNA: return log((1+rna)/(1+dna))/log(2.0);
		case SCORE_RNA_DNA: 
		default:
			return rna/(1+dna);
	}		
}

bool CApp::GetBarcodeScores(const int b0, const int b1, CMatrix<double> &m_scores)
{
	if (b0<0 || b1<0 || b1<=b0)
		return false;
	if (b0>=b1)
	{
		ReportError("Internal error", "b0>=b1");
		return false;
	}
	if (!m_vRNAcols.GetSize())
	{
		ReportError("Interal error", "!m_vRNAcols.GetSize()");
		return false;
	}

	m_scores.ReInit(m_vRNAcols.GetSize(), b1-b0);
	
	int b= b0;	
	while (b<b1)
	{
		int r_barcodes= m_vSI_barcodes[b].m_nIndex;
		double DNA_mu= 0;

		// new code
		for (int j=0;j<m_idx_barcodes_DNA.GetWidth();j++)
		{
			double dna_counts= m_idx_barcodes_DNA.GetAt(r_barcodes, j);
			DNA_mu += dna_counts;
		}
		DNA_mu /= m_vDNAcols.GetSize();

		for (int j=0;j<m_idx_barcodes_RNA.GetWidth();j++)
		{
			double rna_counts= m_idx_barcodes_RNA.GetAt(r_barcodes, j);

			// calculate per-barcode score
			double score= GetBarcodeScore(rna_counts, DNA_mu);
			m_scores.SetAt(j, b-b0, score);
		}

		b++;
	}
	return true;
}

void PermuteBarcodeScores(CMatrix<double> &m_a, CMatrix<double> &m_r)
{
	for (int i=__max(m_a.GetWidth(), m_r.GetWidth());i>=0;i--)
	{
		int c_a= rand_int31() % m_a.GetWidth();
		int c_r= rand_int31() % m_r.GetWidth();
		for (int r=0;r<m_a.GetHeight();r++)
		{
			double tmp= m_a.GetAt(r,c_a);
			m_a.SetAt(r,c_a,m_r.GetAt(r,c_r));
			m_r.SetAt(r,c_r,tmp);
		}
	}
}

void CApp::GetOligoScores(CMatrix<double> &mBarcodeScores_a, CMatrix<double> &mBarcodeScores_r, CMatrix<double> &mB, CMatrix<double> &mW)
{
	// will output unpermuted scores on first row, permuted scores on remaining rows, each column will represent a replicate
	if (mBarcodeScores_a.GetHeight()!=mBarcodeScores_r.GetHeight())
	{
		printf("%d\n", mBarcodeScores_a.GetHeight());
		printf("%d\n", mBarcodeScores_r.GetHeight());
		ReportError("Internal error", "mBarcodeScores_a.GetHeight()!=mBarcodeScores_r.GetHeight()");
		return;
	}

	mB.ReInit(1+m_nNullDistribution_permutations, mBarcodeScores_a.GetHeight());
	mW.ReInit(1+m_nNullDistribution_permutations, mBarcodeScores_a.GetHeight());

	const int n_a= mBarcodeScores_a.GetWidth();
	const int n_r= mBarcodeScores_r.GetWidth();

	for (int p=0;p<mB.GetHeight();p++)
	{
		for (int i=0;i<mBarcodeScores_a.GetHeight();i++)
		{
			double mu_a, var_a;
			Stat_GetMeanAndVariance_double(mBarcodeScores_a.GetPointer(i, 0), n_a, mu_a, var_a);

			double mu_r, var_r;
			Stat_GetMeanAndVariance_double(mBarcodeScores_r.GetPointer(i, 0), n_r, mu_r, var_r);

			// Smyth adjustment
			double n_0= m_vSmyth_DF0[i];
			double var0= m_vSmyth_Var0[i];
			var_a = (n_a*var_a + n_0*var0)/(n_a+n_0);
			var_r = (n_r*var_r + n_0*var0)/(n_r+n_0);

			mB.SetAt(p, i, mu_a-mu_r);
			mW.SetAt(p, i, 1.0/(var_a/n_a + var_r/n_r));
		}

		// permute afterwards to make sure first row is unpermuted
		PermuteBarcodeScores(mBarcodeScores_a, mBarcodeScores_r);
	}
}

bool CApp::GetStatistics_SNP2(Table_SortItem_int *pSI, const int i0, const int i1, const int c_RSID, const int c_Oligo,
	const int c_Ref, const int c_Alt, const int c_Allele_in_oligo, const int c_Strand, const int c_Offset,
	CVector<int> &vOut_r_variant,
	CVector<double> &vOut_fem,
	CVector<double> &vOut_fem_var,
	CVector<double> &vOut_fem_n,
	CVector<double> &vOut_fem_z,
	CVector<double> &vOut_ttest_t,
	CVector<double> &vOut_ttest_df,
	double *p_NULL, int &n_NULL, const int n_NULL_max)
{
	const int I= i1-i0;
	CVector<bool> vAlt(I);
	CVector<bool> vStrand(I);
	CVector<int> vOffset(I);
	CVector<int> v_b0(I);
	CVector<int> v_b1(I);

	// create inventory of oligos
	for (int i=i0;i<i1;i++)
	{
		int r_variants= pSI[i].m_nIndex;

		if (Table_GetAt(m_idx_variants, r_variants, c_Strand)=="+")
			vStrand.SetAt(i-i0, true);
		else if (Table_GetAt(m_idx_variants, r_variants, c_Strand)=="-")
			vStrand.SetAt(i-i0, false);
		else
		{
			ReportError("Illegal strand value", (const char *)Table_GetAt(m_idx_variants, r_variants, c_Strand));
			return false;
		}

		if (Table_CheckIntFormat(m_idx_variants.GetFirstPointer(r_variants, c_Offset)))
			vOffset.SetAt(i-i0, atoi(m_idx_variants.GetFirstPointer(r_variants, c_Offset)));
		else
		{
			ReportError("Illegal offset value", (const char *)Table_GetAt(m_idx_variants, r_variants, c_Offset));
			return false;
		}

		if (Table_GetAt(m_idx_variants, r_variants, c_Allele_in_oligo)==Table_GetAt(m_idx_variants, r_variants, c_Ref))
			vAlt.SetAt(i-i0, false);
		else if (Table_GetAt(m_idx_variants, r_variants, c_Allele_in_oligo)==Table_GetAt(m_idx_variants, r_variants, c_Alt))
			vAlt.SetAt(i-i0, true);
		else
		{
			ReportError("Allele_in_oligo allele must be either Ref or Alt", (const char *)Table_GetAt(m_idx_variants, r_variants, c_Allele_in_oligo));
			return false;
		}

		CString sOligo= Table_GetAt(m_idx_variants, r_variants, c_Oligo);
		int b= FindInIndex_alphabetic(m_vSI_barcodes, sOligo);
		if (b<0)
		{
			// ReportWarning("Oligo not represented in barcode file", sOligo);
			v_b0.SetAt(i-i0, -1);
			v_b1.SetAt(i-i0, -1);
		}
		else
		{
			int b0= b;
			int b1= b;
			while (b0-1>=0 && m_vSI_barcodes[b0-1].m_sKey==m_vSI_barcodes[b].m_sKey)
				b0--;
			while (b1<m_vSI_barcodes.GetSize() && m_vSI_barcodes[b1].m_sKey==m_vSI_barcodes[b].m_sKey)
				b1++;
			v_b0.SetAt(i-i0, b0);
			v_b1.SetAt(i-i0, b1);
		}
	}

	//
	// oligo-level scores for all ref-alt sequence pairs
	//
	// for fixed-effects-model (fem)
	CVector<double> v_s_DeltaMu_weighted(1+m_nNullDistribution_permutations);
	v_s_DeltaMu_weighted.Fill(0);
	double *p_s_DeltaMu_weighted= v_s_DeltaMu_weighted.GetBuffer();
	CVector<double> v_s_weight(1+m_nNullDistribution_permutations);
	v_s_weight.Fill(0);
	double *p_s_weight= v_s_weight.GetBuffer();
	CVector<int> v_n_DeltaMu(1+m_nNullDistribution_permutations);
	v_n_DeltaMu.Fill(0);
	int *p_n_DeltaMu= v_n_DeltaMu.GetBuffer();

	// for paired t
	CVector<double> v_s_DeltaMu_unweighted(1+m_nNullDistribution_permutations);
	v_s_DeltaMu_unweighted.Fill(0);
	double *p_s_DeltaMu_unweighted= v_s_DeltaMu_unweighted.GetBuffer();
	CVector<double> v_ss_DeltaMu_unweighted(1+m_nNullDistribution_permutations);
	v_ss_DeltaMu_unweighted.Fill(0);
	double *p_ss_DeltaMu_unweighted= v_ss_DeltaMu_unweighted.GetBuffer();

	for (int ir=0;ir<I;ir++)
	{
		if (!vAlt[ir])
		{
			for (int ia=0;ia<I;ia++)
				if (vAlt[ia] && vStrand[ia]==vStrand[ir] && vOffset[ia]==vOffset[ir])
				{
					if ((vStrand[ia]==false && m_bIncludeNegStrand) || (vStrand[ia]==true && m_bIncludePosStrand))
					{
						int r_variants_a= pSI[ia+i0].m_nIndex;
						int r_variants_r= pSI[ir+i0].m_nIndex;

						CMatrix<double> mBarcodeScores_a, mBarcodeScores_r; // one replicate per row
						if (GetBarcodeScores(v_b0[ia], v_b1[ia], mBarcodeScores_a) && GetBarcodeScores(v_b0[ir], v_b1[ir], mBarcodeScores_r))
						{
							if (GetVerbosity()>2)
							{
								printf("%s_%s", (const char *)Table_GetAt(m_idx_variants, r_variants_a, c_Oligo), (const char *)Table_GetAt(m_idx_variants, r_variants_r, c_Oligo));
								printf("\t%s", (const char *)Table_GetAt(m_idx_variants, r_variants_a, c_RSID));
								printf("\t%s", (const char *)Table_GetAt(m_idx_variants, r_variants_a, c_Ref));
								printf("\t%s", (const char *)Table_GetAt(m_idx_variants, r_variants_a, c_Alt));
								printf("\t%s", (const char *)Table_GetAt(m_idx_variants, r_variants_a, c_Strand));
								printf("\t%s", (const char *)Table_GetAt(m_idx_variants, r_variants_a, c_Offset));
							}

							CMatrix<double> mDeltaMu, mW; 
							GetOligoScores(mBarcodeScores_a, mBarcodeScores_r, mDeltaMu, mW);
							for (int p=0;p<mDeltaMu.GetHeight();p++) // loop across permutations
							{
								// first row of mB and mW will contain unpermuted scores, the succeeding rows permuted scores
								const double *pDeltaMu= mDeltaMu.GetPointer(p,0);
								const double *pW= mW.GetPointer(p,0);

								for (int i=0;i<mDeltaMu.GetWidth();i++) // loop across oligos (i.e., windows, strands, and replicates)
								{
									if (IsN_double(pDeltaMu[i]) && mBarcodeScores_a.GetWidth()>=m_nMinBarcodesPerOligo && mBarcodeScores_r.GetWidth()>=m_nMinBarcodesPerOligo)
									{
										if (p==0 && GetVerbosity()>2)
										{
											printf("\t%.3g", pDeltaMu[i]);
											printf("\t[%.3g]", pW[i]);
										}

										// for fem
										p_s_weight[p] += pW[i];
										p_s_DeltaMu_weighted[p] += pDeltaMu[i]*pW[i];

										// for paired t
										p_s_DeltaMu_unweighted[p] += pDeltaMu[i];
										p_ss_DeltaMu_unweighted[p] += pDeltaMu[i]*pDeltaMu[i];

										p_n_DeltaMu[p]++;
										//printf("\nsum p_s_weight[%d]= %g, p_n_beta[%d]= %d\n", p, p_s_weight[p], p, p_n_beta[p]);
									}
									else
									{
										if (p==0 && GetVerbosity()>2)
										{
											printf("\t---");
											printf("\t---");
										}
									}
								}

								if (p==0 && GetVerbosity()>2)
									printf("\n");
							}
						}
					}
				}
		}
	}

	// export pooled
	for (int p=0;p<v_s_DeltaMu_weighted.GetSize();p++)
	{
		if (p_n_DeltaMu[p]>0)
		{
			double fem_mu= p_s_DeltaMu_weighted[p]/p_s_weight[p];
			double fem_var= 1.0/p_s_weight[p];
			double z= fem_mu/sqrt(fem_var);

			if (p==0)
			{
				// export to unpermuted distribution
				vOut_r_variant.Add(pSI[i0].m_nIndex);
				vOut_fem.Add(fem_mu);
				vOut_fem_var.Add(fem_var);
				vOut_fem_n.Add(v_n_DeltaMu[p]);
				vOut_fem_z.Add(z);

				double t_mu= p_s_DeltaMu_unweighted[p]/p_n_DeltaMu[p];
				double t_sem= Stat_GetStandardError_double(p_s_DeltaMu_unweighted[p], p_ss_DeltaMu_unweighted[p], p_n_DeltaMu[p])/sqrt(double(p_n_DeltaMu[p]));
				vOut_ttest_t.Add(t_mu/t_sem);
				vOut_ttest_df.Add(p_n_DeltaMu[p]-1);
			}
			else
			{
				// export to null distribution
				p_NULL[n_NULL++]= z;
				if (n_NULL>n_NULL_max)
				{
					ReportError("n_NULL>n_NULL_max", 0);
					return false;
				}
			}
		}
	}

	return true;
}

bool CApp::ExportNullDistribution(const char *aFilename, const double *p_NULL, const int n_NULL)
{
	FILE *pf= OpenOutputFile(aFilename);
	if (!pf)
		return false;

	if (p_NULL!=0)
		for (int i=0;i<n_NULL;i++)
			fprintf(pf, "%g\n", p_NULL[i]);

	CloseFile(pf);
	return true;
}

void GetNullDiagnostics(const double *pNULL, const int nNULL, int &df_opt)
{
	CVector<double> v(nNULL);
	double *p= v.GetBuffer();
	for (int i=0;i<nNULL;i++)
		v.SetAt(i, __abs(pNULL[i]));
	Stat_Quicksort_double(p, nNULL, true);
	
	double r_opt= 0;

	for (int df=2;df<100;df++)
	{
		CVector<double> vqNULL, vqT;
		for (double q=0.00;q<1.00;q+=0.005) // since we start from the middle because of __abs, then q=0 corresponds to 50% cdf
		{
			double z= p[int(nNULL*q)];
			double qT= Stat_Tcdf(z, df);
			double qNULL= 1.0 - (0.5+q/2);
			vqNULL.Add(qNULL);
			vqT.Add(qT);
		}
		double r= Stat_GetCorrelation_double(vqNULL, vqT, vqNULL.GetSize());
		if (r>r_opt)
		{
			r_opt= r;
			df_opt= df;
		}
		// printf("%d\t%g\n", df, r);
	}
}

CVector<double> GetCorrectedP(const CVector<double> &v_p)
{
	double p_mid= Stat_GetMedian_double(v_p);
	double lambda= 0.5/p_mid;
	// printf("Assay inflation factor= %g\n", lambda);

	const int n=v_p.GetSize();
	CVector<double> v_out(n);
	for (int i=0;i<n;i++)
		v_out.SetAt(i, __min(1.0, v_p[i]*lambda));

	return v_out;
}

bool CApp::CreateOutputFile(const char *aFilename_out)
{
	int c_Oligo= Table_FindCol(m_idx_variants, "Oligo");
	if (c_Oligo<0)
	{
		ReportError("Column 'Oligo' not found", 0);
		return false;
	}
	/*
	int c_Sequence= Table_FindCol(m_idx_variants, "Sequence");
	if (c_Sequence<0)
	{
		ReportError("Column 'Sequence' not found", 0);
		return false;
	}
	*/
	int c_RSID= Table_FindCol(m_idx_variants, "RSID");
	if (c_RSID<0)
	{
		ReportError("Column 'RSID' not found", 0);
		return false;
	}
	int c_Ref= Table_FindCol(m_idx_variants, "Ref");
	if (c_Ref<0)
	{
		ReportError("Column 'Ref' not found", 0);
		return false;
	}
	int c_Alt= Table_FindCol(m_idx_variants, "Alt");
	if (c_Alt<0)
	{
		ReportError("Column 'Alt' not found", 0);
		return false;
	}
	int c_Allele_in_oligo= Table_FindCol(m_idx_variants, "Allele_in_oligo");
	if (c_Allele_in_oligo<0)
	{
		ReportError("Column 'Allele_in_oligo' not found", 0);
		return false;
	}
	int c_Strand= Table_FindCol(m_idx_variants, "Strand");
	if (c_Strand<0)
	{
		ReportError("Column 'Strand' not found", 0);
		return false;
	}
	int c_Offset= Table_FindCol(m_idx_variants, "baseOffset");
	if (c_Offset<0)
	{
		ReportError("Column 'baseOffset' not found", 0);
		return false;
	}

	const int R= m_idx_variants.GetRowCount()-1;
	CVector<Table_SortItem_int> vSI(R);
	Table_SortItem_int *pSI= vSI.GetBuffer();

	// sort on RSID
	for (int r=0;r<R;r++)
	{
		Table_SortItem_int si;
		si.m_nIndex= r+1;
		si.m_sKey= Table_GetAt(m_idx_variants, r+1, c_RSID)+"_"+Table_GetAt(m_idx_variants, r+1, c_Ref)+"_"+Table_GetAt(m_idx_variants, r+1, c_Alt);
		si.m_Value= 0;
		pSI[r]= si;
	}
	Table_SortItems_int(pSI, R, true, false);

	CVector<double> v_NULL;
	const int n_NULL_max = m_barcodes_vFirst.GetSize() * m_nNullDistribution_permutations;
	v_NULL.SetSize(n_NULL_max); // max size of null distribution (not all barcodes will succeed...)
	double *p_NULL= v_NULL.GetBuffer();
	int n_NULL= 0; // actual size of null distribution

	CVector<int> v_r_variant;
	CVector<double> v_fem;
	CVector<double> v_fem_var;
	CVector<double> v_fem_n;
	CVector<double> v_fem_z;
	CVector<double> v_ttest_t;
	CVector<double> v_ttest_df;
	CVector<double> v_fem_p, v_fem_p_t, v_fem_p_normal, v_ttest_p;
	CVector<double> v_fem_q, v_fem_q_t, v_fem_q_normal, v_ttest_q;

	const int n_report_frequency= __min(100, __max(5, (100*10000)/m_nNullDistribution_permutations));

	bool bOk= true;
	int n_processed= 0;
	for (int r=0;bOk && r<R;)
	{
		if (n_processed % n_report_frequency==0)
			printf("Processed %d variants...\n", n_processed);
		int r0= r;
		int r1= r0+1;
		while (r1<R && pSI[r1].m_sKey==pSI[r0].m_sKey)
			r1++;

		bOk= GetStatistics_SNP2(pSI, r0, r1, c_RSID, c_Oligo, c_Ref, c_Alt, c_Allele_in_oligo, c_Strand, c_Offset, v_r_variant, v_fem, v_fem_var, v_fem_n, v_fem_z, v_ttest_t, v_ttest_df, p_NULL, n_NULL, n_NULL_max);

		r= r1;
		n_processed++;
	}
	printf("Processed %d variants...\n", n_processed);

	if (bOk && m_bExportNull)
		bOk= ExportNullDistribution(g_FileSystem.GetNameWithoutExtension(aFilename_out)+".null.txt", p_NULL, n_NULL);

	if (bOk)
	{
		PQfromZ(v_fem_z, p_NULL, n_NULL, v_fem_p, v_fem_q);

		// printf("Adjusting p-values: fem, permutation testing\n");
		CVector<double> v_fem_p_corrected= GetCorrectedP(v_fem_p);
		CVector<double> v_fem_q_corrected= QfromP(v_fem_p_corrected, false, false);

		/*
		// p and q based on normal distribution
		v_fem_p_normal.SetSize(v_fem_z.GetSize());
		for (int i=0;i<v_fem_z.GetSize();i++)
		{
			double p= 2*(1.0 - Stat_Normcdf(__abs(v_fem_z[i])));
			v_fem_p_normal.SetAt(i, p);
		}
		printf("Adjusting p-values: fem, normal distribution\n");
		CVector<double> v_fem_p_normal_corrected= GetCorrectedP(v_fem_p_normal);
		v_fem_q_normal= QfromP(v_fem_p_normal, false, false);
		CVector<double>	v_fem_q_normal_corrected= QfromP(v_fem_p_normal_corrected, false, false);

		// p and q based on t distribution
		v_fem_p_t.SetSize(v_fem_z.GetSize());
		for (int i=0;i<v_fem_z.GetSize();i++)
		{
			double p= 2*Stat_Tcdf(__abs(v_fem_z[i]), v_fem_n[i]-1);
			v_fem_p_t.SetAt(i, p);
		}
		v_fem_q_t= QfromP(v_fem_p_t, false, false);
		printf("Adjusting p-values: fem, t distribution\n");
		CVector<double> v_fem_p_t_corrected= GetCorrectedP(v_fem_p_t);
		v_ttest_q= QfromP(v_fem_p_t, false, false);
		CVector<double>	v_fem_q_t_corrected= QfromP(v_fem_p_t_corrected, false, false);
		*/

		// p and q for paired t-test
		v_ttest_p.SetSize(v_fem_z.GetSize());
		for (int i=0;i<v_fem_z.GetSize();i++)
		{
			double p= 2*Stat_Tcdf(__abs(v_ttest_t[i]), v_ttest_df[i]);
			v_ttest_p.SetAt(i, p);
		}
		// printf("Adjusting p-values: paired t-test\n");
		CVector<double> v_ttest_p_corrected= GetCorrectedP(v_ttest_p);
		v_ttest_q= QfromP(v_ttest_p, false, false);
		CVector<double>	v_ttest_q_corrected= QfromP(v_ttest_p_corrected, false, false);

		// sort output
		CVector<Stat_SortItem_double> vSI(v_r_variant.GetSize());
		Stat_SortItem_double *pSI= vSI.GetBuffer();
		for (int i=0;i<v_r_variant.GetSize();i++)
		{
			Stat_SortItem_double si;
			si.m_nIndex= i;
			si.m_Value= v_fem_z[i]; 
			pSI[i]= si;
		}
		Stat_QuicksortIndexed_double(pSI, vSI.GetSize(), false);

		// write to file
		FILE *pf= OpenOutputFile(aFilename_out);
		if (!pf)
			return false;

		if (m_bExportHeader)
		{
			fprintf(pf, "RSID\tRef\tAlt\tn\tscore\tp_perm\tq_perm\tp_perm_adj\tq_perm_adj");
			if (m_bExportTTest)
				fprintf(pf, "\tp_t\tq_t\tp_t_adj\tq_t_adj");
			fprintf(pf, "\n");
		}

		for (int j=0;j<v_r_variant.GetSize();j++)
		{
			int i= pSI[j].m_nIndex;
			if (m_bExportHeader)
			{
				fprintf(pf, "%s\t", (const char *)Table_GetAt(m_idx_variants, v_r_variant[i], c_RSID));
				fprintf(pf, "%s\t", (const char *)Table_GetAt(m_idx_variants, v_r_variant[i], c_Ref));
				fprintf(pf, "%s\t", (const char *)Table_GetAt(m_idx_variants, v_r_variant[i], c_Alt));
			}
			fprintf(pf, "%g", v_fem_n[i]);
			fprintf(pf, "\t%g", v_fem[i]);
			// fprintf(pf, "\t%g", sqrt(v_fem_var[i]));
			// fprintf(pf, "\t%g", v_fem_z[i]);
			fprintf(pf, "\t%g", v_fem_p[i]);
			fprintf(pf, "\t%g", v_fem_q[i]);
			fprintf(pf, "\t%g", v_fem_p_corrected[i]);
			fprintf(pf, "\t%g", v_fem_q_corrected[i]);
			/*
			fprintf(pf, "\t%g", v_fem_p_normal[i]);
			fprintf(pf, "\t%g", v_fem_q_normal[i]);
			fprintf(pf, "\t%g", v_fem_p_normal_corrected[i]);
			fprintf(pf, "\t%g", v_fem_q_normal_corrected[i]);
			fprintf(pf, "\t%g", v_fem_p_t[i]);
			fprintf(pf, "\t%g", v_fem_q_t[i]);
			fprintf(pf, "\t%g", v_fem_p_t_corrected[i]);
			fprintf(pf, "\t%g", v_fem_q_t_corrected[i]);
			*/
			if (m_bExportTTest)
			{
				fprintf(pf, "\t%g", v_ttest_p[i]);
				fprintf(pf, "\t%g", v_ttest_q[i]);
				fprintf(pf, "\t%g", v_ttest_p_corrected[i]);
				fprintf(pf, "\t%g", v_ttest_q_corrected[i]);
			}
			fprintf(pf, "\n");
		}
		CloseFile(pf);
	}

	return bOk;
}

/******************************************************************************/
// Main application

int CApp::Main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRes= 100;
	SetDefaultParameters();
	if (ParseCommandLine(argc, argv, envp) && !m_bHelpShown)
	{
		if (m_vNonSwitchArgs.GetSize()<3)
		{
			ReportError("Not enough arguments", 0);
			return 1000;
		}

		// file 0: barcode map
		// file 1: barcode counts
		// file 2: output
		CString sErr;
		printf("Loading barcode map...\n");
		CLoadBuf *pLB0= Table_LoadTabDataset(m_vNonSwitchArgs[0], m_idx_variants, sErr);
		if (pLB0)
		{
			printf("Loading barcode counts...\n");
			CLoadBuf *pLB1= Table_LoadTabDataset(m_vNonSwitchArgs[1], m_idx_barcodes, sErr);
			if (pLB1)
			{
				if (CheckParameters())
				{
					if (CreateBarcodeLookup() && GetBarcodeStats())
					{
						if (CreateOutputFile(m_vNonSwitchArgs[2]))
							nRes= 0;
					}
				}
			}
			else
				ReportError(sErr, (const char *)m_vNonSwitchArgs[1]);

			if (pLB1)
				g_FileSystem.FreeBuf(pLB1);
		}
		else
			ReportError(sErr, (const char *)m_vNonSwitchArgs[0]);

		if (pLB0)
			g_FileSystem.FreeBuf(pLB0);
	}

	return nRes;
}

/******************************************************************************/

int main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	CApp app;
	return app.Main(argc, argv, envp);
}
