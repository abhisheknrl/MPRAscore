/*******************************************************************************
 *
 *   statfun.h -- 
 *
 *   Björn Nilsson, 2004-
 *
 */

#ifndef __STATFUN_H__
#define __STATFUN_H__

#include <math.h>
#include "types/vector.h"
#include "types/matrix.h"

inline float Stat_abs(float x) { return x>=0 ? x : -x; };
inline int Stat_fix(float x) { return x>=0 ? int(x) : -int(-x) ; };
inline float Stat_max(float x, float y) { return x>y ? x : y; };
inline float Stat_min(float x, float y) { return x<y ? x : y; };
inline double Stat_abs(double x) { if (x<0) return -x; else return x; };
inline int Stat_fix(double x) { return x>=0 ? int(x) : -int(-x) ; };
inline double Stat_max(double x, double y) { return x>y ? x : y; };
inline double Stat_min(double x, double y) { return x<y ? x : y; };

// double Stat_GetSum(const float *aData, const int nSize);
double Stat_GetSum_double(const double *aData, const int nSize);
double Stat_GetSum_int(const int *aData, const int nSize);
// void Stat_GetSums(const float *aData, const int nSize, double &sum, double &ssqr);
void Stat_GetSums_double(const double *aData, const int nSize, double &sum, double &ssqr);
// double Stat_GetMean(const float *aData, const int nSize);
double Stat_GetMean_double(const double *aData, const int nSize);
double Stat_GetMean_int(const int *aData, const int nSize);

// void Stat_Normalize(float *aData, int nSize);
// float Stat_GetMedian(float *aData, int nSize);
// float Stat_GetMedian(const CVector<float> &vData);
// void Stat_GetMeanAndVariance(const float *aData, const int nSize, double &mean, double &var);
// void Stat_GetMinAndMax(const float *aData, const int nSize, float &min, float &max);

double Stat_GetWeightedMean_double(const float *aData, const float *aWeight, const int nSize);
double Stat_GetMedian_double(double *aData, int nSize);
double Stat_GetMedian_double(const CVector<double> &vData);
void Stat_GetMeanAndVariance_double(const double *aData, const int nSize, double &mean, double &var);
void Stat_GetMinAndMax_double(const double *aData, const int nSize, double &min, double &max);
float Stat_GetMax(const float *aData, const int nSize);
float Stat_GetMin(const float *aData, const int nSize);
double Stat_GetMax_double(const double *aData, const int nSize);
double Stat_GetMin_double(const double *aData, const int nSize);
int Stat_GetMax(const int *aData, const int nSize);
int Stat_GetMin(const int *aData, const int nSize);
// float Stat_GetMeanSqrError(const float *aData0, const float *aData1, int nSize);
double Stat_GetMeanSqrError_double(const double *aData0, const double *aData1, int nSize);
void Stat_NormalizeL2_double(double *aData, int nSize);
void Stat_Normalize_CenterUnitVariance_double(double *aData, int nSize);
void Stat_QuantileNormalize_Normal_double(double *aData, int nSize);
void Stat_QuantileNormalize_Rank_double(double *aData, int nSize);
void Stat_LogTransform_double(double *aData, double bias, int nSize);
void Stat_GetLinreg_double(const double *aDataX, const double *aDataY, int nSize, double &fCoeff, double &fBias);
void Stat_GetWeightedLinreg(const float *aX, const float *aY, const float *aW, int nSize, float &fCoeff, float &fBias);
double Stat_GetInnerProduct(const double *aData0, const double *aData1, const int n);
double Stat_GetCorrelation(double s0, double s1, double s00, double s01, double s11, int n);
// double Stat_GetCorrelation(const float *aData0, const float *aData1, int nSize);
double Stat_GetCorrelation_double(const double *aData0, const double *aData1, int nSize);
// double Stat_GetCovariance(const float *aData0, const float *aData1, int nSize);
double Stat_GetCovariance_double(const double *aData0, const double *aData1, int nSize);
double Stat_GetCorrelation_Weighted(const float *aData0, const float *aData1, const float *aW, int nSize);
float Stat_GetCorrelation_Fast(const float *aData0, const float *aData1, int nSize);
// double Stat_GetRawCorrelation(const float *aData0, const float *aData1, int nSize);
double Stat_GetRawCorrelation_double(const float *aData0, const float *aData1, int nSize);
// double Stat_GetL1Norm(const float *aData, int nSize);
double Stat_GetL1Norm_double(const double *aData, int nSize);
// double Stat_GetL1Distance(const float *aData0, const float *aData1, int nSize);
// double Stat_GetL2Norm(const float *aData, int nSize);
double Stat_GetL2Norm_double(const double *aData, int nSize);
// double Stat_GetL2Distance(const float *aData0, const float *aData1, int nSize);
double Stat_GetL2Distance_double(const double *aData0, const double *aData1, int nSize);
// double Stat_GetL2Diff(const float *aData0, const float *aData1, int nSize);
double Stat_GetL2Diff_double(const double *aData0, const double *aData1, int nSize);
void Stat_GramSchmidt_double(const double *p0, const double *p1, double *pout, int nSize);

// float Stat_GetPooledVariance(const float s0, const float ss0, int n0, const float s1, const float ss1, int n1);
double Stat_GetPooledVariance_double(const double s0, const double ss0, int n0, const double s1, const double ss1, int n1);
// float Stat_T(const float s0, const float ss0, int n0, const float s1, const float ss1, int n1);
double Stat_T_double(const double s0, const double ss0, int n0, const double s1, const double ss1, int n1);
float Stat_SAM(const float s0, const float ss0, int n0, 
			   const float s1, const float ss1, int n1, float fudge);


// Non-indexed sorting
//void Stat_Quicksort(float *aData, int nSize, bool bAscending);
void Stat_Quicksort_int(int *aData, int nSize, bool bAscending);
void Stat_Quicksort_double(double *aData, int nSize, bool bAscending);

// Indexed sorting
/*
// deprecated 2017-08-29, use double versions/BN
struct Stat_SortItem
{
	float m_fValue;
	int m_nIndex;
};
void Stat_QuicksortIndexed(Stat_SortItem *aSI, int nSize, bool Ascending);
Stat_SortItem *Stat_InitSortItem(CVector<Stat_SortItem> &vSI, const CVector<float> &v, bool bAscending);
*/

struct Stat_SortItem_double
{
	double m_Value;
	int m_nIndex;
};
void Stat_QuicksortIndexed_double(Stat_SortItem_double *aSI, int nSize, bool bAscending);
void Stat_QuicksortIndexed_double(CVector<Stat_SortItem_double> &vSI, bool bAscending);
Stat_SortItem_double *Stat_InitSortItem_double(CVector<Stat_SortItem_double> &vSI, const CVector<double> &v, bool bAscending);

struct Stat_SortItem_int
{
	int m_nValue;
	int m_nIndex;
};
void Stat_QuicksortIndexed_int(Stat_SortItem_int *aSI, int nSize, bool Ascending);
Stat_SortItem_int *Stat_InitSortItem_int(CVector<Stat_SortItem_int> &vSI, const CVector<int> &v, bool bAscending);

//
// float Stat_GetQuantile(const float *aData, int nSize, float q);
// float Stat_GetQuantile(const CVector<float> &vData, float q);

double Stat_GetQuantile_double(const double *aData, int nSize, double q);

inline double Stat_GetVariance(const double s, const double ssqr, const int nSize)
{	
	if (nSize<=1) 
		return 0; 
	else
		return (ssqr-s*s/nSize)/(nSize-1); 
}

/*
inline double Stat_GetVariance(const float *aData, int nSize)
{ 
	double m, v; 
	Stat_GetMeanAndVariance(aData, nSize, m, v);
	return v;
}
*/

inline double Stat_GetVariance_double(const double *aData, int nSize)
{ 
	double m, v; 
	Stat_GetMeanAndVariance_double(aData, nSize, m, v);
	return v;
}

inline double Stat_GetStandardError_double(const double s, const double ssqr, const int nSize)
{ return (double)sqrt((ssqr-s*s/nSize)/(nSize-1)); }
inline double Stat_GetStandardError_double(const double *aData, const int nSize)
{ return (double) sqrt(Stat_GetVariance_double(aData, nSize)); }

float Stat_GetCohen_double(double m0, double v0, double n0, double m1, double v1, double n1);

// Inverse chi-square test
float Stat_Fisher(const float *aP, const int n);
double Stat_pochisq (double x, int df);
double Stat_poz (double z);

// Advanced functions related to the computation of distribution functions
double Stat_Gamma(double x);
double Stat_Gammaln(double x);
double Stat_Gammainc(double x, double a, bool lower);
double Stat_Beta(double x, double y);
double Stat_Betainc(double x, double a, double b);
double Stat_Betacore(double x, double a, double b);
double Stat_erfc(double x);
double Stat_Tcdf(double x, double df);
double Stat_Normcdf(double z);
double Stat_InvErf(double p);
double Stat_InvNormcdf(double p);

// columnwise pre-normalization
#define NORM_None (0)
#define NORM_Rank (1)
#define NORM_Normal (2)
#define NORM_Lognormal (3)
#define NORM_Last (NORM_Lognormal)
bool Stat_NormalizeMatrix(CMatrix<double> &m, int mode);
bool Stat_NormalizeMatrix(CMatrix<double> &m, int mode, int cmin, int cmax);

// NaN handling, float, deprecated 2017-08-29/BN
// #define float_Exp (0x7f800000)
// const long float_NaN= float_Exp | 1;
// const long float_PosInfty= float_Exp;
// const long float_NegInfty= float_Exp | 0x80000000;
// inline float Stat_GetNaN_float() { return *(float *)&float_NaN; }
// inline float Stat_GetPosInfty_float() { return *(float *)&float_PosInfty; }
// inline float Stat_GetNegInfty_float() { return *(float *)&float_NegInfty; }
// inline bool IsN_float(float x) { return ((*(long *)&x) & float_Exp)!=float_Exp; }
// inline bool IsNaN_float(float x) { return ((*(long *)&x) & float_Exp)==float_Exp; }

// NaN handling, double
double Stat_GetNaN_double();
double Stat_GetPosInfty_double();
double Stat_GetNegInfty_double();

#ifdef _WIN32
// isnan is not implemented in Windows distribution of math.h
#define double_Exp (0x7ff0000000000000)
inline bool IsN_double(double x) { return ((*(long long *)&x) & double_Exp)!=double_Exp; }
inline bool IsNaN_double(double x) { return ((*(long long *)&x) & double_Exp)==double_Exp; }
#else
#include <math.h>
inline bool IsNaN_double(double x) { return isnan(x)!=0; }
inline bool IsN_double(double x) { return isnan(x)==0; }
#endif

/*
#ifdef _WIN32
// Visual Studio
#define double_Exp (0x7ff0000000000000)
const unsigned long long double_NaN= 0x7ff0000000000001; // double_Exp | 1;
const unsigned long long double_PosInfty= double_Exp;
const unsigned long long double_NegInfty= double_Exp | 0x8000000000000000;
inline double Stat_GetNaN_double() { return *(double *)&double_NaN; }
inline double Stat_GetPosInfty_double() { return *(double *)&double_PosInfty; }
inline double Stat_GetNegInfty_double() { return *(double *)&double_NegInfty; }
inline bool IsN_double(double x) { return ((*(long long *)&x) & double_Exp)!=double_Exp; }
inline bool IsNaN_double(double x) { return ((*(long long *)&x) & double_Exp)==double_Exp; }
#else 
// GCC has problems with 64-bit longs -- TODO: Fix this properly!!!
#define double_Exp (0x7ff0000000000000)
//const unsigned long long double_NaN= 0x7ff0000000000001; // double_Exp | 1;
//const unsigned long long double_PosInfty= double_Exp;
//const unsigned long long double_NegInfty= double_Exp | 0x8000000000000000;
inline double Stat_GetNaN_double() { return double(Stat_GetNaN_float()); }
inline double Stat_GetPosInfty_double() { return double(Stat_GetPosInfty_float()); }
inline double Stat_GetNegInfty_double() { return double(Stat_GetNegInfty_float()); }
inline bool IsN_double(double x) { return IsN_float(x); }
inline bool IsNaN_double(double x) { return IsNaN_float(x); }
#endif
*/

bool Stat_HistConvolve(const CVector<double> &h0, double min0, double max0, 
					   const CVector<double> &h1, double min1, double max1, 
					   CVector<double> &h_out, double &min_out, double &max_out);
void Stat_HistNormalize(CVector<double> &h);

#ifdef _WIN32
inline double Acklam_stdnormal_pdf(double u);
double Acklam_stdnormal_cdf(double u);
double Acklam_stdnormal_inv(double p);
#endif

#endif

