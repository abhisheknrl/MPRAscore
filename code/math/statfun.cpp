/*******************************************************************************
 *
 *   statfun.cpp -- Basic statistics functions
 *
 *   Björn Nilsson, 2004--
 */

#include "statfun.h"
#include <stdlib.h> // For rand
#include "math/rand.h"
#include <float.h>
#include "types/nan.h"
#include "types/minmax.h"

// TODO: UNIX: Detect and handle big-endian format at run-time

// This trick avoids the 32-bit limit of the gcc preprocessor
unsigned char g_NaN_double[8]= { 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x7f };
unsigned char g_PosInfty_double[8]= { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x7f };
unsigned char g_NegInfty_double[8]= { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xff };

double Stat_GetNaN_double() { return *(double *)&g_NaN_double; }
double Stat_GetPosInfty_double() { return *(double *)&g_PosInfty_double; }
double Stat_GetNegInfty_double() { return *(double *)&g_NegInfty_double; }

/******************************************************************************/
// Basic functions for single and double precision floating point numbers
// TODO: Move to separate file?

#define pi (3.1415926535897932384)
#define Stat_eps (2.2204e-16)
/******************************************************************************/
// Basic statistics, including sums, means, norms, variances, quantiles

/*
double Stat_GetSum(const float *aData, const int nSize)
{
	// Get sum(x)
	double s= 0;
	for (int i=nSize;--i>=0;)
		s += *aData++;
	return s;
}
*/

double Stat_GetSum_double(const double *aData, const int nSize)
{
	// Get sum(x)
	double s= 0;
	for (int i=nSize;--i>=0;)
		s += *aData++;
	return s;
}

double Stat_GetSum_int(const int *aData, const int nSize)
{
	// Get sum(x)
	int s= 0;
	for (int i=nSize;--i>=0;)
		s += *aData++;
	return double(s);
}

/*
void Stat_GetSums(const float *aData, const int nSize, double &s, double &ss)
{
	// Compute sum and sum-of-squares of aData[nSize]
	s= 0;
	ss= 0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= *aData++;
		s += tmp;
		ss += tmp*tmp;
	}
}

void Stat_GetSums(const float *aData, const int nSize, double &s, double &ss)
{
	// Compute sum and sum-of-squares of aData[nSize]
	s= 0;
	ss= 0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= *aData++;
		s += tmp;
		ss += tmp*tmp;
	}
}
*/

void Stat_GetSums_double(const double *aData, const int nSize, double &s, double &ss)
{
	// Compute sum and sum-of-squares of aData[nSize]
	s= 0;
	ss= 0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= *aData++;
		s += tmp;
		ss += tmp*tmp;
	}
}

/*
double Stat_GetMean(const float *aData, const int nSize)
{
	// Returns the standard estimate of the mean.
	ASSERT(nSize>=0);
	if (nSize<=0)
		return 0;
	return Stat_GetSum(aData, nSize)/nSize;
}
*/

double Stat_GetMean_double(const double *aData, const int nSize)
{
	// Returns the standard estimate of the mean.
	if (nSize<=0)
		return 0;
	return Stat_GetSum_double(aData, nSize)/nSize;
}

double Stat_GetMean_int(const int *aData, const int nSize)
{
	// Returns the standard estimate of the mean.
	if (nSize<=0)
		return 0;
	return double(Stat_GetSum_int(aData, nSize))/nSize;
}

/*
float Stat_GetWeightedMean(const float *aData, const float *aWeight, const int nSize)
{
	double swx= 0;
	double sw= 0;
	for (int i=nSize;--i>=0;)
	{
		swx += aWeight[i]*aData[i];
		sw += aWeight[i];
	}

	if (sw>0)
		return float(swx/sw);
	else
		return 0;
}
*/

double Stat_GetWeightedMean_double(const double *aData, const double *aWeight, const int nSize)
{
	double swx= 0;
	double sw= 0;
	for (int i=nSize;--i>=0;)
	{
		swx += aWeight[i]*aData[i];
		sw += aWeight[i];
	}

	if (sw>0)
		return swx/sw;
	else
		return 0;
}

/*
float Stat_GetMedian(float *aData, int nSize)
{
	// Returns the median and a sorted vector.
	if (!aData || nSize<=0)
		return 0;
	Stat_Quicksort(aData, nSize, true);
	int i= nSize >> 1;
	if ((nSize & 1)==0)
		// Vector has even length. Interpolate.
		return (aData[i-1] + aData[i])/2;
	else
		// Vector has odd length.
		return aData[i];
}

float Stat_GetMedian(const CVector<float> &vData)
{
	CVector<float> v_tmp= vData;
	return Stat_GetMedian(v_tmp.GetBuffer(), v_tmp.GetSize());
}
*/

double Stat_GetMedian_double(double *aData, int nSize)
{
	// Returns the median and a sorted vector.
	if (!aData || nSize<=0)
		return 0;
	Stat_Quicksort_double(aData, nSize, true);
	int i= nSize >> 1;
	if ((nSize & 1)==0)
		// Vector has even length. Interpolate.
		return (aData[i-1] + aData[i])/2;
	else
		// Vector has odd length.
		return aData[i];
}

double Stat_GetMedian_double(const CVector<double> &vData)
{
	CVector<double> v_tmp= vData;
	return Stat_GetMedian_double(v_tmp.GetBuffer(), v_tmp.GetSize());
}

/*
void Stat_GetMeanAndVariance(const float *aData, const int nSize, double &mean, double &var)
{
	// Special case, small vector
	if (nSize<=1)
	{
		var= 0;
		if (nSize)
			mean= *aData;
		else
			mean= 0;
		return;
	}

	double s, ssqr;
	Stat_GetSums(aData, nSize, s, ssqr);

	mean= s/nSize;
	var= Stat_GetVariance(s, ssqr, nSize);
}
*/

void Stat_GetMeanAndVariance_double(const double *aData, const int nSize, double &mean, double &var)
{
	// Special case, small vector
	if (nSize<=1)
	{
		var= 0;
		if (nSize)
			mean= *aData;
		else
			mean= 0;
		return;
	}

	double s, ssqr;
	Stat_GetSums_double(aData, nSize, s, ssqr);

	mean= s/nSize;
	var= Stat_GetVariance(s, ssqr, nSize);
}

/*
float Stat_GetQuantile(const float *aData, int nSize, float q)
{
	// NOTE: Assumes that the data have been sorted in ascending order.
	if (!aData)
		return 0;
	if (nSize==1)
		return aData[0];
	if (q<0.0f)
		return aData[0]; // q>1.0 handled below

	float qi= float(nSize-1)*q;
	int i= int(qi);
	float t= qi - float(i);

	if (i>=nSize-1)
		return aData[nSize-1];
	else
		return (1-t)*aData[i] + t*aData[i+1];
}

float Stat_GetQuantile(const CVector<float> &vData, float q)
{
	CVector<float> v_tmp= vData;
	float *p_tmp= v_tmp.GetBuffer();
	Stat_Quicksort(p_tmp, v_tmp.GetSize(), true);
	return Stat_GetQuantile(p_tmp, v_tmp.GetSize(), q);
}
*/

double Stat_GetQuantile_double(const double *aData, int nSize, double q)
{
	// NOTE: Assumes that the data have been sorted in ascending order.
	if (!aData)
		return 0;
	if (nSize==1)
		return aData[0];
	if (q<0.0f)
		return aData[0]; // q>1.0 handled below

	double qi= double(nSize-1)*q;
	int i= int(qi);
	double t= qi - double(i);

	if (i>=nSize-1)
		return aData[nSize-1];
	else
		return (1.0-t)*aData[i] + t*aData[i+1];
}

/*
void Stat_GetMinAndMax(const float *aData, const int nSize, float &min, float &max)
{
	if (nSize<=0)
	{
		min= 0;
		max= 0;
		return;
	}

	min= *aData;
	max= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData<min)
			min= *aData;
		else if (*aData>max)
			max= *aData;
		aData++;
	}
}
*/

void Stat_GetMinAndMax_double(const double *aData, const int nSize, double &min, double &max)
{
	if (nSize<=0)
	{
		min= 0;
		max= 0;
		return;
	}

	min= *aData;
	max= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData<min)
			min= *aData;
		else if (*aData>max)
			max= *aData;
		aData++;
	}
}

/*
float Stat_GetMax(const float *aData, const int nSize)
{
	if (nSize<=0)
		return 0;

	float max= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData>max)
			max= *aData;
		aData++;
	}
	return max;
}

float Stat_GetMin(const float *aData, const int nSize)
{
	if (nSize<=0)
		return 0;

	float min= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData<min)
			min= *aData;
		aData++;
	}
	return min;
}
*/

double Stat_GetMax_double(const double *aData, const int nSize)
{
	if (nSize<=0)
		return 0;

	double max= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData>max)
			max= *aData;
		aData++;
	}
	return max;
}

double Stat_GetMin_double(const double *aData, const int nSize)
{
	if (nSize<=0)
		return 0;

	double min= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData<min)
			min= *aData;
		aData++;
	}
	return min;
}

int Stat_GetMax(const int *aData, const int nSize)
{
	if (nSize<=0)
		return 0;

	int max= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData>max)
			max= *aData;
		aData++;
	}
	return max;
}


int Stat_GetMin(const int *aData, const int nSize)
{
	if (nSize<=0)
		return 0;

	int min= *aData++;
	for (int i=nSize-1;--i>=0;)
	{
		if (*aData<min)
			min= *aData;
		aData++;
	}
	return min;
}

/*
// deprecated 2017-08-29, use double. /BN
void Stat_Normalize(float *aData, int nSize)
{
	double mu= Stat_GetMean(aData, nSize);
	for (int i=nSize;--i>=0;)
		aData[i] -= mu;

	double sigma= Stat_GetL2Norm(aData, nSize);
	for (int i=nSize;--i>=0;)
		aData[i] /= sigma;
}
*/

void Stat_NormalizeL2_double(double *aData, int nSize)
{
	// previously called Stat_Normalize_double, which was unclear.
	double mu= Stat_GetMean_double(aData, nSize);
	for (int i=nSize;--i>=0;)
		aData[i] -= mu;

	double sigma= Stat_GetL2Norm_double(aData, nSize);
	for (int i=nSize;--i>=0;)
		aData[i] /= sigma;
}

void Stat_Normalize_CenterUnitVariance_double(double *aData, int nSize)
{
	double mu= Stat_GetMean_double(aData, nSize);
	for (int i=nSize;--i>=0;)
		aData[i] -= mu;

	if (nSize>1)
	{
		double sd= sqrt(Stat_GetVariance_double(aData, nSize)); // Stat_GetL2Norm_double(aData, nSize)/(nSize-1));
		for (int i=nSize;--i>=0;)
			aData[i] /= sd;
	}
}

void Stat_QuantileNormalize_Normal_double(double *aData, int nSize)
{
	if (nSize<=0)
		return;
	CVector<Stat_SortItem_double> vSI(nSize);
	Stat_SortItem_double *pSI= vSI.GetBuffer();
	for (int i=0;i<nSize;i++)
	{
		Stat_SortItem_double si;
		si.m_nIndex= i;
		si.m_Value= aData[i];
		pSI[i]= si;
	}
	Stat_QuicksortIndexed_double(pSI, nSize, true);
	
	for (int i=0;i<nSize;i++)
	{
		double y0= aData[pSI[i].m_nIndex];
		double y1= Stat_InvNormcdf(double(i+1)/(nSize+1));
		// printf("%d\t%g\t%g\n", i, y0, y1);
		aData[pSI[i].m_nIndex]= y1;
	}
}

void Stat_QuantileNormalize_Rank_double(double *aData, int nSize)
{
	if (nSize<=0)
		return;
	CVector<Stat_SortItem_double> vSI(nSize);
	Stat_SortItem_double *pSI= vSI.GetBuffer();
	for (int i=0;i<nSize;i++)
	{
		Stat_SortItem_double si;
		si.m_nIndex= i;
		si.m_Value= aData[i];
		pSI[i]= si;
	}
	Stat_QuicksortIndexed_double(pSI, nSize, true);
	
	for (int i=0;i<nSize;i++)
		aData[pSI[i].m_nIndex]= double(i)/nSize;
}

void Stat_LogTransform_double(double *aData, double bias, int nSize)
{
	for (int i=0;i<nSize;i++)
		aData[i]= log(aData[i]+bias);
}

/*
double Stat_GetL1Norm(const float *aData, int nSize)
{
	double s= 0;
	for (int i=nSize;--i>=0;)
		aData[i]>=0 ? s += aData[i] : s -= aData[i];
	return s;
}
*/

double Stat_GetL1Norm_double(const double *aData, int nSize)
{
	double s= 0;
	for (int i=nSize;--i>=0;)
		aData[i]>=0 ? s += aData[i] : s -= aData[i];
	return s;
}

/*
double Stat_GetL1Distance(const float *aData0, const float *aData1, int nSize)
{
	double s= 0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= aData0[i]-aData1[i];
		tmp>=0 ? s += tmp : s -= tmp;
	}
	return s;
}

double Stat_GetL2Norm(const float *aData, int nSize)
{
	double ss=0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= aData[i];
		ss += tmp*tmp;
	}
	return sqrt(ss);
}
*/

double Stat_GetL2Norm_double(const double *aData, int nSize)
{
	double ss=0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= aData[i];
		ss += tmp*tmp;
	}
	return sqrt(ss);
}

/*
double Stat_GetL2Diff(const float *aData0, const float *aData1, int nSize)
{
	double ss= 0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= aData0[i]-aData1[i];
		ss += tmp*tmp;
	}
	return ss;
}
*/

double Stat_GetL2Diff_double(const double *aData0, const double *aData1, int nSize)
{
	double ss= 0;
	for (int i=nSize;--i>=0;)
	{
		double tmp= aData0[i]-aData1[i];
		ss += tmp*tmp;
	}
	return ss;
}

double Stat_GetScalarProduct_double(const double *p0, const double *p1, int nSize)
{
	double s=0;
	for (int i=0;i<nSize;i++)
		s += p0[i]*p1[i];
	return s;
}

void Stat_GramSchmidt_double(const double *p0, const double *p1, double *pout, int nSize)
{
	// Makes p0 orthogonal with respect to p1.
	// Stores the result in pout, which can be the same as p0 or p1.

	double p0p1= Stat_GetScalarProduct_double(p0, p1, nSize);
	double p1norm= Stat_GetL2Norm_double(p1, nSize);
	if (p1norm!=0)
	{
		p0p1 /= p1norm;

		// printf("p0p1= %g\n", p0p1);
		for (int i=0;i<nSize;i++)
			pout[i]= p0[i] -p1[i]*p0p1;
	}
	else
		for (int i=0;i<nSize;i++)
			pout[i]= p0[i];
}

/*
double Stat_GetL2Distance(const float *aData0, const float *aData1, int nSize)
{
	return sqrt(Stat_GetL2Diff(aData0, aData1, nSize));
}
*/

double Stat_GetL2Distance_double(const double *aData0, const double *aData1, int nSize)
{
	return sqrt(Stat_GetL2Diff_double(aData0, aData1, nSize));
}

/*
float Stat_GetMeanSqrError(const float *aData0, const float *aData1, int nSize)
{
	float s= 0;
	if (aData0 && aData1 && nSize>=0)
	{
		for (int i=nSize;--i>=0;)
		{
			float tmp= (*aData0++)-(*aData1++);
			s += tmp*tmp;
		}
		s /= nSize;
	}
	return s;
}
*/

double Stat_GetMeanSqrError_double(const double *aData0, const double *aData1, int nSize)
{
	double s= 0;
	if (aData0 && aData1 && nSize>=0)
	{
		for (int i=nSize;--i>=0;)
		{
			double tmp= (*aData0++)-(*aData1++);
			s += tmp*tmp;
		}
		s /= nSize;
	}
	return s;
}

/*
float Stat_GetPooledVariance(const float s0, const float ss0, int n0, 
							 const float s1, const float ss1, int n1)
{
	if (n0<=0 || n1<=0)
		return 0.0f;

	if (n0+n1>2)
	{
		double sigmasqr0= Stat_GetVariance(s0, ss0, n0);
		double sigmasqr1= Stat_GetVariance(s1, ss1, n1);
		return float((sigmasqr0*(n0-1) + sigmasqr1*(n1-1))/float(n0+n1-2));
	}
	return 1.0f;
}
*/

double Stat_GetPooledVariance_double(const double s0, const double ss0, int n0, 
									 const double s1, const double ss1, int n1)
{
	if (n0<=0 || n1<=0)
		return 0.0f;

	if (n0+n1>2)
	{
		double sigmasqr0= Stat_GetVariance(s0, ss0, n0);
		double sigmasqr1= Stat_GetVariance(s1, ss1, n1);
		return (sigmasqr0*(n0-1) + sigmasqr1*(n1-1))/float(n0+n1-2);
	}
	return 1.0f;
}

/*
float Stat_T(const float s0, const float ss0, int n0, 
			 const float s1, const float ss1, int n1)
{
	float nom= s0/n0 + s1/n1;
	float denom= (float) sqrt(Stat_GetPooledVariance(s0, ss0, n0, s1, ss1, n1)*(1/float(n0) + 1/float(n1)));
	return nom/denom;
}
*/

double Stat_T_double(const double s0, const double ss0, int n0, 
					 const double s1, const double ss1, int n1)
{
	double nom= s0/n0 + s1/n1;
	double denom= (double) sqrt(Stat_GetPooledVariance_double(s0, ss0, n0, s1, ss1, n1)*(1/double(n0) + 1/double(n1)));
	return nom/denom;
}

/*
float Stat_SAM(const float s0, const float ss0, int n0, 
			   const float s1, const float ss1, int n1,
			   const float fudge)
{
	float nom= s0/n0 + s1/n1;
	float denom= (float) sqrt(Stat_GetPooledVariance(s0, ss0, n0, s1, ss1, n1)*(1/float(n0) + 1/float(n1)));
	return nom/(denom + fudge);
}
*/

/******************************************************************************/
// Regression

/*
void Stat_GetLinreg(const float *aDataX, const float *aDataY, int nSize, double &fCoeff, double &fBias)
{
	// Get least square regression parameters.
	// Vector too small?
	if (nSize<2)
	{
		fCoeff= 0;
		fBias= 0;
		return; 
	}

	double fX, fXsqr;
	Stat_GetSums(aDataX, nSize, fX, fXsqr);
	double fY= Stat_GetSum(aDataY, nSize);

	double fXY= 0;
	for (int i=nSize;--i>=0;)
		fXY += aDataX[i]*aDataY[i];

	double SSxx= fXsqr - fX*fX/nSize;
	double SSxy= fXY - fX*fY/nSize;

	fCoeff= SSxy/SSxx;
	fBias= fY/nSize - fCoeff*fX/nSize;
}
*/

void Stat_GetLinreg_double(const double *aDataX, const double *aDataY, int nSize, double &fCoeff, double &fBias)
{
	// Get least square regression parameters.
	// Vector too small?
	if (nSize<2)
	{
		fCoeff= 0;
		fBias= 0;
		return; 
	}

	double fX, fXsqr;
	Stat_GetSums_double(aDataX, nSize, fX, fXsqr);
	double fY= Stat_GetSum_double(aDataY, nSize);

	double fXY= 0;
	for (int i=nSize;--i>=0;)
		fXY += aDataX[i]*aDataY[i];

	double SSxx= fXsqr - fX*fX/nSize;
	double SSxy= fXY - fX*fY/nSize;

	fCoeff= SSxy/SSxx;
	fBias= fY/nSize - fCoeff*fX/nSize;
}

void Stat_GetWeightedLinreg(const float *aX, const float *aY, const float *aW, int nSize, float &fCoeff, float &fBias)
{
	// Returns weighted least square linear fit.
	double w= 0;
	double wx= 0;
	double wy= 0;
	double wxx= 0;
	double wxy= 0;

	for (int i=nSize;--i>=0;)
	{
		w += aW[i];
		float wixi= aW[i]*aX[i];
		wx += wixi;
		wxx += wixi*aX[i];
		wxy += wixi*aY[i];
		wy += aW[i]*aY[i];
	}

	fCoeff= float((wy*wx-w*wxy)/( wx*wx - wxx*w ));
	fBias= float((wy*wxx-wx*wxy)/( w*wxx - wx*wx ));
}

double Stat_GetCorrelation(double s0, double s1, double s00, double s01, double s11, int n)
{
	double d= sqrt((s00 - s0*s0/n)*(s11 - s1*s1/n));
	if (d>0.0)
	{
		double corr= (s01-s0*s1/n)/d;
		if (corr>1.0) // protect against round-off errors
			return 1.0;
		if (corr<-1.0)
			return -1.0;
		return corr;
	}
	else
		return 0.0;
}

/*
double Stat_GetCorrelation(const float *aData0, const float *aData1, int nSize)
{
	float m0= float(Stat_GetMean(aData0, nSize));
	float m1= float(Stat_GetMean(aData1, nSize));

	double s00= 0;
	double s01= 0;
	double s11= 0;
	for (int i=0;i<nSize;i++)
	{
		float d0= aData0[i] - m0;
		float d1= aData1[i] - m1;
		s00 += d0*d0;
		s01 += d0*d1;
		s11 += d1*d1;
	}

	const double eps= 1.0e-10;
	if (s00*s11<eps)
		return 0.0;
	else
		return s01/sqrt(s00*s11);
}
*/

double Stat_GetCorrelation_double(const double *aData0, const double *aData1, int nSize)
{
	double m0= Stat_GetMean_double(aData0, nSize);
	double m1= Stat_GetMean_double(aData1, nSize);

	double s00= 0;
	double s01= 0;
	double s11= 0;
	for (int i=0;i<nSize;i++)
	{
		double d0= aData0[i] - m0;
		double d1= aData1[i] - m1;
		s00 += d0*d0;
		s01 += d0*d1;
		s11 += d1*d1;
	}

	const double eps= 1.0e-10;
	if (s00*s11<eps)
		return 0.0;
	else
		return s01/sqrt(s00*s11);
}

/*
double Stat_GetCovariance(const float *aData0, const float *aData1, int nSize)
{
	float m0= float(Stat_GetMean(aData0, nSize));
	float m1= float(Stat_GetMean(aData1, nSize));

	double s00= 0;
	double s01= 0;
	double s11= 0;
	for (int i=0;i<nSize;i++)
	{
		float d0= aData0[i] - m0;
		float d1= aData1[i] - m1;
		s01 += d0*d1;
	}

	return s01/(nSize-1);
}
*/

double Stat_GetCovariance_double(const double *aData0, const double *aData1, int nSize)
{
	double m0= Stat_GetMean_double(aData0, nSize);
	double m1= Stat_GetMean_double(aData1, nSize);

	double s00= 0;
	double s01= 0;
	double s11= 0;
	for (int i=0;i<nSize;i++)
	{
		double d0= aData0[i] - m0;
		double d1= aData1[i] - m1;
		s01 += d0*d1;
	}

	return s01/(nSize-1);
}

/*
double Stat_GetCorrelation(const float *aData0, const float *aData1, int nSize)
{
	// Returns Pearson's correlation coefficient
	double s0, s1, ss0, ss1;
	Stat_GetSums(aData0, nSize, s0, ss0);
	ASSERT(IsN_double(s0));
	ASSERT(IsN_double(ss0));
	Stat_GetSums(aData1, nSize, s1, ss1);
	ASSERT(IsN_double(s1));
	ASSERT(IsN_double(ss1));

	double d= sqrt((ss0 - s0*s0/nSize)*(ss1 - s1*s1/nSize));
	ASSERT(IsN_double(d));

	double s01= 0;
	for (int i=nSize;--i>=0;)
		s01 += aData0[i]*aData1[i];

	return (s01-s0*s1/nSize)/d;
}
*/

double Stat_GetCorrelation_Weighted(const float *aX, const float *aY, const float *aW, int n)
{
	// aX[], aY[] data
	// aW[] weights, sum assumed to be one 

	double mu_wx= 0;
	double mu_wy= 0;
	for (int i=0;i<n;i++)
	{
		mu_wx += aW[i]*aX[i];
		mu_wy += aW[i]*aY[i];
	}
	// no division as sum(aW[])==1

	double cov_xy= 0;
	double cov_xx= 0;
	double cov_yy= 0;
	for (int i=0;i<n;i++)
	{
		double dx= aX[i]-mu_wx;
		double dy= aY[i]-mu_wy;
		cov_xx += aW[i]*dx*dx;
		cov_xy += aW[i]*dx*dy;
		cov_yy += aW[i]*dy*dy;
	}

	if (cov_xx>0 && cov_yy>0)
		return cov_xy/sqrt(cov_xx*cov_yy);
	else
		return 0.0;
}

float Stat_GetCorrelation_Fast(const float *aData0, const float *aData1, int nSize)
{
	// Returns Pearson's correlation coefficient using all-float computations
	float s0= 0;
	float s1= 0;
	float ss0= 0;
	float ss1= 0;
	for (int i=nSize;--i>=0;)
	{
		float tmp= aData0[i];
		s0 += tmp;
		ss0 += tmp*tmp;
		tmp= aData1[i];
		s1 += tmp;
		ss1 += tmp*tmp;
	}

	float d= (float) sqrt((ss0 - s0*s0/nSize)*(ss1 - s1*s1/nSize));

	float s01= 0;
	for (int i=nSize;--i>=0;)
		s01 += aData0[i]*aData1[i];

	return (s01-s0*s1/nSize)/d;
}

/*
double Stat_GetRawCorrelation(const float *aData0, const float *aData1, int nSize)
{
	// Returns the raw correlation of the input vectors
	double s0, s1, ss0, ss1;
	Stat_GetSums(aData0, nSize, s0, ss0);
	Stat_GetSums(aData1, nSize, s1, ss1);

	double s01= 0;
	for (int i=nSize;--i>=0;)
		s01 += aData0[i]*aData1[i];

	return s01/sqrt(ss0)/sqrt(ss1);
}
*/

double Stat_GetRawCorrelation_double(const double *aData0, const double *aData1, int nSize)
{
	// Returns the raw correlation of the input vectors
	double s0, s1, ss0, ss1;
	Stat_GetSums_double(aData0, nSize, s0, ss0);
	Stat_GetSums_double(aData1, nSize, s1, ss1);

	double s01= 0;
	for (int i=nSize;--i>=0;)
		s01 += aData0[i]*aData1[i];

	return s01/sqrt(ss0)/sqrt(ss1);
}

double Stat_InverseChiSquare(const double chi, int nDF)
{
	// Returns P(chisq >= chi), with nDF (even) degrees of freedom.

	double m= chi / 2.0;
	double sum= exp(-m);
	double term= sum;
	for (int i=1;i<=nDF/2;i++)
	{
		term *= m/i;
		sum += term;
	}

	// Guard against round-off errors.
	if (sum<=1.0)
		return sum;
	else
		return 1.0;
}

/******************************************************************************/
// Sorting 

void Quicksort_Ascending(float *pLeft, float *pRight)
{
	if (pRight<=pLeft)
		return;

	float *p_i= pLeft-1;
	float *p_j= pRight+1;

	float v= *pLeft;
	float tmp;

	do 
	{
		while (*(++p_i)<v);
		while (*(--p_j)>v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Ascending(pLeft, p_j);
	Quicksort_Ascending(p_j+1, pRight);
}

void Quicksort_Descending(float *pLeft, float *pRight)
{
	if (pRight<=pLeft)
		return;

	float *p_i= pLeft-1;
	float *p_j= pRight+1;

	float v= *pLeft;
	float tmp;

	do 
	{
		while (*(++p_i)>v);
		while (*(--p_j)<v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Descending(pLeft, p_j);
	Quicksort_Descending(p_j+1, pRight);
}

void RandomizeItems(float *aData, int nSize)
{	
	if (nSize<100)
		return; // Vector is small, don't randomize

	float tmp;
	srand(1);
	int I= nSize/3; // 33% permut
	for (int i=0;i<I;i++)
	{
		int j0= (rand_int31() % nSize);
		int j1= (rand_int31() % nSize);
		tmp= aData[j0];
		aData[j0]= aData[j1];
		aData[j1]= tmp;
	}
}

void Stat_Quicksort(float *aData, int nSize, bool bAscending)
{
	RandomizeItems(aData, nSize);
	
	if (bAscending)
		Quicksort_Ascending(aData, aData + (nSize - 1));
	else
		Quicksort_Descending(aData, aData + (nSize - 1));
}


void Quicksort_double_Ascending(double *pLeft, double *pRight)
{
	if (pRight<=pLeft)
		return;

	double *p_i= pLeft-1;
	double *p_j= pRight+1;

	double v= *pLeft;
	double tmp;

	do 
	{
		while (*(++p_i)<v);
		while (*(--p_j)>v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_double_Ascending(pLeft, p_j);
	Quicksort_double_Ascending(p_j+1, pRight);
}

void Quicksort_double_Descending(double *pLeft, double *pRight)
{
	if (pRight<=pLeft)
		return;

	double *p_i= pLeft-1;
	double *p_j= pRight+1;

	double v= *pLeft;
	double tmp;

	do 
	{
		while (*(++p_i)>v);
		while (*(--p_j)<v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_double_Descending(pLeft, p_j);
	Quicksort_double_Descending(p_j+1, pRight);
}

void Stat_Quicksort_double(double *aData, int nSize, bool bAscending)
{
	if (bAscending)
		Quicksort_double_Ascending(aData, aData + (nSize - 1));
	else
		Quicksort_double_Descending(aData, aData + (nSize - 1));
}

void Quicksort_int_Ascending(int *pLeft, int *pRight)
{
	if (pRight<=pLeft)
		return;

	int *p_i= pLeft-1;
	int *p_j= pRight+1;

	int v= *pLeft;
	int tmp;

	do 
	{
		while (*(++p_i)<v);
		while (*(--p_j)>v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_int_Ascending(pLeft, p_j);
	Quicksort_int_Ascending(p_j+1, pRight);
}

void Quicksort_int_Descending(int *pLeft, int *pRight)
{
	if (pRight<=pLeft)
		return;

	int *p_i= pLeft-1;
	int *p_j= pRight+1;
	int v= *pLeft;
	int tmp;

	do 
	{
		while (*(++p_i)>v);
		while (*(--p_j)<v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_int_Descending(pLeft, p_j);
	Quicksort_int_Descending(p_j+1, pRight);
}

void Stat_Quicksort_int(int *aData, int nSize, bool bAscending)
{
	if (bAscending)
		Quicksort_int_Ascending(aData, aData + (nSize - 1));
	else
		Quicksort_int_Descending(aData, aData + (nSize - 1));
}

/*
void Quicksort_Indexed_Ascending_double(Stat_SortItem_double *pLeft, Stat_SortItem_double *pRight)
{
	if (pRight<=pLeft)
		return;

	Stat_SortItem_double *p_i= pLeft-1;
	Stat_SortItem_double *p_j= pRight+1;

	double v= pLeft->m_Value;
	Stat_SortItem_double tmp;

	do 
	{
		while ((++p_i)->m_Value<v);
		while ((--p_j)->m_Value>v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Indexed_Ascending_double(pLeft, p_j);
	Quicksort_Indexed_Ascending_double(p_j+1, pRight);
}

void Quicksort_Indexed_Descending_double(Stat_SortItem_double *pLeft, Stat_SortItem_double *pRight)
{
	if (pRight<=pLeft)
		return;

	Stat_SortItem_double *p_i= pLeft-1;
	Stat_SortItem_double *p_j= pRight+1;

	double v= pLeft->m_Value;
	Stat_SortItem_double tmp;

	do 
	{
		while ((++p_i)->m_Value>v);
		while ((--p_j)->m_Value<v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Indexed_Descending_double(pLeft, p_j);
	Quicksort_Indexed_Descending_double(p_j+1, pRight);
}

void Stat_QuicksortIndexed_double(Stat_SortItem_double *aSI, int nSize, bool bAscending)
{
	if (bAscending)
		Quicksort_Indexed_Ascending_double(aSI, aSI + (nSize - 1));
	else
		Quicksort_Indexed_Descending_double(aSI, aSI + (nSize - 1));
}

Stat_SortItem *Stat_InitSortItem(CVector<Stat_SortItem> &vSI, const CVector<float> &v, bool bAscending)
{
	vSI.SetSize(v.GetSize());
	for (int i=0;i<v.GetSize();i++)
	{
		Stat_SortItem si;
		si.m_fValue= v[i];
		si.m_nIndex= i;
		vSI.SetAt(i, si);
	}

	Stat_SortItem *pSI= vSI.GetBuffer();
	Stat_QuicksortIndexed(pSI, vSI.GetSize(), bAscending);
	return pSI;
}
*/

void Quicksort_Indexed_Ascending_double(Stat_SortItem_double *pLeft, Stat_SortItem_double *pRight)
{
	if (pRight<=pLeft)
		return;

	Stat_SortItem_double *p_i= pLeft-1;
	Stat_SortItem_double *p_j= pRight+1;

	double v= pLeft->m_Value;
	Stat_SortItem_double tmp;

	do 
	{
		while ((++p_i)->m_Value<v);
		while ((--p_j)->m_Value>v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Indexed_Ascending_double(pLeft, p_j);
	Quicksort_Indexed_Ascending_double(p_j+1, pRight);
}

void Quicksort_Indexed_Descending_double(Stat_SortItem_double *pLeft, Stat_SortItem_double *pRight)
{
	if (pRight<=pLeft)
		return;

	Stat_SortItem_double *p_i= pLeft-1;
	Stat_SortItem_double *p_j= pRight+1;

	double v= pLeft->m_Value;
	Stat_SortItem_double tmp;

	do 
	{
		while ((++p_i)->m_Value>v);
		while ((--p_j)->m_Value<v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Indexed_Descending_double(pLeft, p_j);
	Quicksort_Indexed_Descending_double(p_j+1, pRight);
}

void Stat_QuicksortIndexed_double(Stat_SortItem_double *aSI, int nSize, bool bAscending)
{
	if (bAscending)
		Quicksort_Indexed_Ascending_double(aSI, aSI + (nSize - 1));
	else
		Quicksort_Indexed_Descending_double(aSI, aSI + (nSize - 1));
}

void Stat_QuicksortIndexed_double(CVector<Stat_SortItem_double> &vSI, bool bAscending)
{
	Stat_SortItem_double *pSI= vSI.GetBuffer();
	Stat_QuicksortIndexed_double(pSI, vSI.GetSize(), bAscending);
}

Stat_SortItem_double *Stat_InitSortItem_double(CVector<Stat_SortItem_double> &vSI, const CVector<double> &v, bool bAscending)
{
	vSI.SetSize(v.GetSize());
	for (int i=0;i<v.GetSize();i++)
	{
		Stat_SortItem_double si;
		si.m_Value= v[i];
		si.m_nIndex= i;
		vSI.SetAt(i, si);
	}

	Stat_SortItem_double *pSI= vSI.GetBuffer();
	Stat_QuicksortIndexed_double(pSI, vSI.GetSize(), bAscending);
	return pSI;
}

void Quicksort_Indexed_Ascending_int(Stat_SortItem_int *pLeft, Stat_SortItem_int *pRight)
{
	if (pRight<=pLeft)
		return;

	Stat_SortItem_int *p_i= pLeft-1;
	Stat_SortItem_int *p_j= pRight+1;

	int v= pLeft->m_nValue;
	Stat_SortItem_int tmp;

	do 
	{
		while ((++p_i)->m_nValue<v);
		while ((--p_j)->m_nValue>v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Indexed_Ascending_int(pLeft, p_j);
	Quicksort_Indexed_Ascending_int(p_j+1, pRight);
}

void Quicksort_Indexed_Descending_int(Stat_SortItem_int *pLeft, Stat_SortItem_int *pRight)
{
	if (pRight<=pLeft)
		return;

	Stat_SortItem_int *p_i= pLeft-1;
	Stat_SortItem_int *p_j= pRight+1;

	int v= pLeft->m_nValue;
	Stat_SortItem_int tmp;

	do 
	{
		while ((++p_i)->m_nValue>v);
		while ((--p_j)->m_nValue<v);

		tmp= *p_i;
		*p_i= *p_j;
		*p_j= tmp;
	}
	while (p_i<p_j);

	tmp= *p_i;
	*p_i= *p_j;
	*p_j= tmp;

	Quicksort_Indexed_Descending_int(pLeft, p_j);
	Quicksort_Indexed_Descending_int(p_j+1, pRight);
}

void Stat_QuicksortIndexed_int(Stat_SortItem_int *aSI, int nSize, bool bAscending)
{
	if (bAscending)
		Quicksort_Indexed_Ascending_int(aSI, aSI + (nSize - 1));
	else
		Quicksort_Indexed_Descending_int(aSI, aSI + (nSize - 1));
}

Stat_SortItem_int *Stat_InitSortItem_int(CVector<Stat_SortItem_int> &vSI, const CVector<int> &v, bool bAscending)
{
	vSI.SetSize(v.GetSize());
	for (int i=0;i<v.GetSize();i++)
	{
		Stat_SortItem_int si;
		si.m_nValue= v[i];
		si.m_nIndex= i;
		vSI.SetAt(i, si);
	}

	Stat_SortItem_int *pSI= vSI.GetBuffer();
	Stat_QuicksortIndexed_int(pSI, v.GetSize(), bAscending);
	return pSI;
}


/******************************************************************************/
// FUNCTION poz: probability of normal z value
// Returns P(Z<'z'), Z a N(0,1) RV.
// Adapted from a polynomial approximation in:
//  Ibbetson D, Algorithm 209
//  Collected Algorithms of the CACM 1963 p. 616
// Note:
//  This routine has six digit accuracy, so it is only useful for absolute
//    z values < 6.  For z values >= to 6.0, poz() returns 0.0.
/******************************************************************************/

double Stat_poz (double z)
{
  double x = 0.0;
  
  if (z != 0.0)
  {
    double y = 0.5 * fabs (z);
    if (y >= 3.0)
      x = 1.0;
    else if (y < 1.0) 
	{
      double w = y*y;
      x = ((((((((0.000124818987 * w
        -0.001075204047) * w +0.005198775019) * w
        -0.019198292004) * w +0.059054035642) * w
        -0.151968751364) * w +0.319152932694) * w
        -0.531923007300) * w +0.797884560593) * y * 2.0;
    }
    else 
	{
      y -= 2.0;
      x = (((((((((((((-0.000045255659 * y
        +0.000152529290) * y -0.000019538132) * y
        -0.000676904986) * y +0.001390604284) * y
        -0.000794620820) * y -0.002034254874) * y
        +0.006549791214) * y -0.010557625006) * y
        +0.011630447319) * y -0.009279453341) * y
        +0.005353579108) * y -0.002141268741) * y
        +0.000535310849) * y +0.999936657524;
    }
  }

  return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

/******************************************************************************/
// FUNCTION pochisq: probability of chi square value 
// ALGORITHM Compute probability of chi square value.
//	Adapted from:
//		Hill, I. D. and Pike, M. C.  Algorithm 299
//		Collected Algorithms for the CACM 1967 p. 243
//	Updated for rounding errors based on remark in
//		ACM TOMS June 1985, page 185
/******************************************************************************/

// log (sqrt (pi))
#define	LOG_SQRT_PI     0.5723649429247000870717135 

// 1 / sqrt (pi) 
#define	I_SQRT_PI       0.5641895835477562869480795 

// max value to represent exp (x)
#define	BIGX           20.0         

#define	ex(x)             (((x) < -BIGX) ? 0.0 : exp (x))

double Stat_pochisq (double x, int df)
{
	double y= Stat_GetNaN_double(); // Avoids gcc warning 
	double	a, s;
	double	e, c, z;
	int 	even;     // true if df is an even number
	
	if (x <= 0.0 || df < 1)
		return (1.0);
	
	a = 0.5 * x;
	even = (2*(df/2)) == df;
	if (df > 1)
		y = ex(-a);
	s = (even ? y : (2.0 * Stat_poz (-sqrt (x))));
	if (df > 2)
	{
		x = 0.5 * (df - 1.0);
		z = (even ? 1.0 : 0.5);
		if (a > BIGX)
		{
			e = (even ? 0.0 : LOG_SQRT_PI);
			c = log (a);
			while (z <= x)
			{
				e = log (z) + e;
				s += ex(c*z-a-e);
				z += 1.0;
			}
			return (s);
		}
		else
		{
			e = (even ? 1.0 : (I_SQRT_PI / sqrt (a)));
			c = 0.0;
			while (z <= x)
			{
				e = e * (a / z);
				c = c + e;
				z += 1.0;
			}
			return (c * y + s);
		}
	}
	else
		return (s);
}

float Stat_Fisher(const float *aP, const int n)
{
	// Implements Fisher's inverse chi-square algorithm 
	// for combining p-values from independent tests.
	if (n<=0)
		return 0;

	double S= 0;
	for (int i=0;i<n;i++)
		if (aP[i]>=0 && aP[i]<=1)
			S -= log(aP[i]);
		else
			return 0;
	S *= 2;
	return (float)Stat_pochisq(S, 2*n);
}

//**********************************************************************

double Stat_erfcore(double x, int jint)
{
	// erfcore function
	// Results verified 2008-05-19

	double result;

	// evaluate  erf  for  |x| <= 0.46875
	double xbreak = 0.46875;
	if (Stat_abs(x)<=xbreak)
	{
		double a[5]= {
			3.16112374387056560e00, 1.13864154151050156e02,
			3.77485237685302021e02, 3.20937758913846947e03,
			1.85777706184603153e-1};
		double b[4]= {
			2.36012909523441209e01, 2.44024637934444173e02,
			1.28261652607737228e03, 2.84423683343917062e03};

		double y= Stat_abs(x);
		double z= y*y;
		double xnum= a[4]*z;
		double xden= z;
		for (int i=0;i<3;i++)
		{
			xnum= (xnum+a[i])*z;
			xden= (xden+b[i])*z;
		}
		result= x*(xnum+a[3])/(xden+b[3]);
		if (jint!=-0)
			result= 1 - result;
		else if (jint==2)
			result= exp(z)*result;
	}
	else if (Stat_abs(x)<=4.)
	{
		// evaluate  erfc  for 0.46875 <= |x| <= 4.0
		const double c[9]= {
			5.64188496988670089e-1, 8.88314979438837594e00,
			6.61191906371416295e01, 2.98635138197400131e02,
			8.81952221241769090e02, 1.71204761263407058e03,
			2.05107837782607147e03, 1.23033935479799725e03,
			2.15311535474403846e-8};
		const double d[8]= {
			1.57449261107098347e01, 1.17693950891312499e02,
			5.37181101862009858e02, 1.62138957456669019e03,
			3.29079923573345963e03, 4.36261909014324716e03,
			3.43936767414372164e03, 1.23033935480374942e03};

		double y = Stat_abs(x);
		double xnum= c[8]*y;
		double xden= y;
		for (int i=0;i<7;i++)
		{
		   xnum= (xnum + c[i]) * y;
		   xden= (xden + d[i]) * y;
		}
		result= (xnum + c[7])/(xden + d[7]);
		if (jint!=2)
		{
		   double z= Stat_fix(y*16.0)/16.0;
		   double del = (y-z)*(y+z);
		   result= exp(-z*z) * exp(-del) * result;
		}
	}
	else
	{
		// evaluate  erfc  for |x| > 4.0
		const double p[6] = {
			3.05326634961232344e-1, 3.60344899949804439e-1,
			1.25781726111229246e-1, 1.60837851487422766e-2,
			6.58749161529837803e-4, 1.63153871373020978e-2};
		const double q[5] = {
			2.56852019228982242e00, 1.87295284992346047e00,
			5.27905102951428412e-1, 6.05183413124413191e-2,
			2.33520497626869185e-3};

		double y= Stat_abs(x);
		double z= 1 / (y * y);
		double xnum= p[5]*z;
		double xden= z;
		for (int i=0;i<4;i++)
		{
		   xnum= (xnum + p[i]) * z;
		   xden= (xden + q[i]) * z;
		}
		result= z * (xnum + p[4]) / (xden + q[4]);
		result= (1.0/sqrt(pi) -  result) / y;
		if (jint != 2)
		{
		   z = Stat_fix(y*16)/16;
		   double del = (y-z)*(y+z);
		   result= exp(-z*z) * exp(-del) * result;
		   if (result>FLT_MAX)
			   result= 0;
		}
	}

	// fix up for negative argument, erf, etc.
	switch (jint)
	{
		case 0:
			if (x > xbreak)
				result= (0.5 - result) + 0.5;
			else if (x < -xbreak)
				result= (-0.5 + result) - 0.5;
			break;
		case 1:
			if (x < -xbreak)
				result= 2. - result;
			break;
		case 2:
			if (x < -xbreak)
			{
				double z = Stat_fix(x*16)/16;
				double del = (x-z)*(x+z);
				double y = exp(z*z) * exp(del);
				result= (y+y) - result;
			}
			break;
		default:
			result= Stat_GetNaN_double(); // Illegal jint value
	}

	return result;
}

double Stat_erfc(double x)
{
	return Stat_erfcore(x, 1);
}

double Stat_Gammainc(double x, double a, bool lower)
{
	// Incomplete gamma function
	// Input: x real, a non-negative real, lower=true for lower tail, lower=false for upper tail.

	/*
	%    gammainc(x,a) = 1 ./ gamma(a) .*
	%       integral from 0 to x of t^(a-1) exp(-t) dt
	%
	%   For any a>=0, as x approaches infinity, gammainc(x,a) approaches 1.
	%   For small x and a, gammainc(x,a) ~= x^a, so gammainc(0,0) = 1.
	%
	%   Y = GAMMAINC(X,A,TAIL) specifies the tail of the incomplete gamma
	%   function when X is non-negative.  Choices are 'lower' (the default)
	%   and 'upper'.  The upper incomplete gamma function is defined as
	%   1 - gammainc(x,a).
	%
	%   Warning: When X is negative, Y can be inaccurate for abs(X) > A+1.
	*/

	if (a<0)
		return Stat_GetNaN_double();

	// Upper limit for series and continued fraction
	double amax = exp(log(2.0)*20);

	// Approximation for a > amax.  Accurate to about 5.e-5.
	if (a>amax)
	{
	  x= Stat_max(amax-1.0/3 + sqrt(amax/a)*(x-(a-1.0/3)),0);
	  a= amax;
	}

	// Series expansion for x < a+1
	if (a!=0 && x!=0 && x<a+1)
	{
		double xk= x;
		double ak = a;
		double ap= ak;
		double del= 1;
		double sum= del;
		//while norm(del,'inf') >= 100*eps(norm(sum,'inf'))
		while (Stat_abs(del)>=100*Stat_eps*Stat_abs(sum))
		{
			ap= ap + 1;
			del = xk * del / ap;
			sum = sum + del;
		}
		double bk = sum * exp(-xk + ak*log(xk) - Stat_Gammaln(ak+1));
		
		// For very small a, the series may overshoot very slightly.
		// bk(xk > 0 & bk > 1) = 1;
		if (xk>0 && bk>1)
			bk= 1;
		return lower ? bk : 1-bk;
	}

	// Continued fraction for x >= a+1
	if (a!=0 && x>=a+1)
	{
		double xk= x;
		double a0 = 1;
		double a1 = x;
		double b0 = 0;
		double b1 = a0;
		double ak = a;
		double fac = 1.0/a1;
		int n = 1;
		double g = b1*fac;
		double gold = b0;
		while (Stat_abs(g-gold) >= 100*Stat_eps*Stat_abs(g))
		{
			gold = g;
			double ana = n - ak;
			a0 = (a1 + a0 *ana) * fac;
			b0 = (b1 + b0 *ana) * fac;
			double anf = n*fac;
			a1 = xk * a0 + anf * a1;
			b1 = xk * b0 + anf * b1;
			fac = 1.0 / a1;
			g = b1 * fac;
			n = n + 1;
		}
		double bk = exp(-xk + ak*log(xk) - Stat_Gammaln(ak)) * g;
		return lower ? 1-bk : bk;
	}

	if (x==0)
		return lower ? 0 : 1;
	if (a==0)
		return lower ? 1 : 0;

	return Stat_GetNaN_double();
}

double Stat_Gammaln(double x)
{
	// Gammaln function, defined by gammaln(x) = ln(gamma(x))
	double d1 = -5.772156649015328605195174e-1;
	double p1[8] = {4.945235359296727046734888e0, 2.018112620856775083915565e2,
					2.290838373831346393026739e3, 1.131967205903380828685045e4,
					2.855724635671635335736389e4, 3.848496228443793359990269e4,
					2.637748787624195437963534e4, 7.225813979700288197698961e3};
	double q1[8] = {6.748212550303777196073036e1, 1.113332393857199323513008e3,
					7.738757056935398733233834e3, 2.763987074403340708898585e4,
					5.499310206226157329794414e4, 6.161122180066002127833352e4,
					3.635127591501940507276287e4, 8.785536302431013170870835e3};

	double d2 = 4.227843350984671393993777e-1;
	double p2[8] = {4.974607845568932035012064e0, 5.424138599891070494101986e2,
					1.550693864978364947665077e4, 1.847932904445632425417223e5,
					1.088204769468828767498470e6, 3.338152967987029735917223e6,
					5.106661678927352456275255e6, 3.074109054850539556250927e6};
	double q2[8] = {1.830328399370592604055942e2, 7.765049321445005871323047e3,
					1.331903827966074194402448e5, 1.136705821321969608938755e6,
					5.267964117437946917577538e6, 1.346701454311101692290052e7,
					1.782736530353274213975932e7, 9.533095591844353613395747e6};

	double d4 = 1.791759469228055000094023e0;
	double p4[8] = {1.474502166059939948905062e4, 2.426813369486704502836312e6,
					1.214755574045093227939592e8, 2.663432449630976949898078e9,
					2.940378956634553899906876e10, 1.702665737765398868392998e11,
					4.926125793377430887588120e11, 5.606251856223951465078242e11};
	double q4[8] = {2.690530175870899333379843e3, 6.393885654300092398984238e5,
					4.135599930241388052042842e7, 1.120872109616147941376570e9,
					1.488613728678813811542398e10, 1.016803586272438228077304e11,
					3.417476345507377132798597e11, 4.463158187419713286462081e11};
	double  c[7] = {-1.910444077728e-03, 8.4171387781295e-04,
					-5.952379913043012e-04, 7.93650793500350248e-04,
					-2.777777777777681622553e-03, 8.333333333333333331554247e-02,
					5.7083835261e-03};

	if (x<=0) // Illegal value
		return 0.0f;
	
	if (x<=Stat_eps)
		return -log(x);
	
	if (x<=0.5)
	{
		double xden= 1.0;
		double xnum= 0.0;
		for (int i=0;i<8;i++)
		{
			xnum= xnum*x + p1[i];
			xden= xden*x + q1[i];
		}
		return -log(x)+ x*(d1+x*(xnum/xden));
	}

	if (x<=0.6796875)
	{
		double xm1= (x-0.5)-0.5;
		double xden= 1.0;
		double xnum= 0.0;
		for (int i=0;i<8;i++)
		{
			xnum= xnum*xm1 + p2[i];
			xden= xden*xm1 + q2[i];
		}
		return -log(x) + xm1*(d2 + xm1*(xnum/xden));
	}
	
	if (x<=1.5)
	{
		double xm1= (x-0.5)-0.5;
		double xden= 1.0;
		double xnum= 0.0;
		for (int i=0;i<8;i++)
		{
			xnum= xnum*xm1 + p1[i];
			xden= xden*xm1 + q1[i];
		}
		return xm1*(d1 + xm1*(xnum/xden));
	}

	if (x<=4.0)
	{
		double xm2= x-2.0;
		double xden= 1.0;
		double xnum= 0.0;
		for (int i=0;i<8;i++)
		{
			xnum= xnum*xm2 + p2[i];
			xden= xden*xm2 + q2[i];
		}
		return xm2*(d2 + xm2*(xnum/xden));
	}

	if (x<=12.0)
	{
		double xm4= x-4.0;
		double xden= -1.0;
		double xnum= 0.0;
		for (int i=0;i<8;i++)
		{
			xnum= xnum*xm4 + p4[i];
			xden= xden*xm4 + q4[i];
		}
		return d4 + xm4*(xnum/xden);
	}

	double r= c[6];
	double ysq= x*x;
	for (int i=0;i<6;i++)
		r= r/ysq+c[i];
	r= r/x;
	double corr= log(x);
	double spi= 0.9189385332046727417803297;
	return r+spi-0.5*corr+x*(corr-1.0);
}

double Stat_Gamma(double x)
{
	// Gamma function
	return exp(Stat_Gammaln(x));
}

double Stat_Beta(double x, double y)
{
	// Beta function
	return exp(Stat_Gammaln(x)+Stat_Gammaln(y)-Stat_Gammaln(x+y));
}

double Stat_Betaln(double x, double y)
{
	// Betaln function
	return Stat_Gammaln(x)+Stat_Gammaln(y)-Stat_Gammaln(x+y);
}

double Stat_Betacore(double x, double a, double b)
{
	/*
	function y = betacore(x, a, b)
	%BETACORE Core algorithm for the incomplete beta function.
	%   Y = BETACORE(X,A,B) computes a continued fraction expansion used by
	%   BETAINC.  Specifically,
	%
	%      BETAINC(X,A,B) = BETACORE(X,A,B) * (X^A * (1-X)^B) / (A*BETA(A,B)).
	%
	%   X must be strictly between 0 and 1.  Returns NaN if continued fraction
	%   does not converge.

	aplusb = a + b;
	aplus1 = a + 1;
	aminus1 = a - 1;
	C = 1;
	% When called from BETAINC, Dinv can never be zero unless (a+b) or (a+1)
	% round to a.
	Dinv = 1 - aplusb .* x ./ aplus1;
	y = C ./ Dinv;
	*/
	double aplusb= a+b;
	double aplus1= a+1;
	double aminus1= a-1;
	double C=1.0;
	double Dinv= 1.0-aplusb*x/aplus1;
	double y= C/Dinv;

	/*
	maxiter = 1000;
	for m = 1:maxiter
		yold = y;
		twom = 2 * m;
		d = m * (b - m) .* x ./ ((aminus1 + twom) .* (a + twom));
		C = 1 + d ./ C;
		% Using Dinv, not D, ensures that C = 1/D will be a stable fixed point
		Dinv = 1 + d ./ Dinv;
		y = y .* (C./Dinv);
		d = -(a + m) .* (aplusb + m) .* x ./ ((a + twom) .* (aplus1 + twom));
		C = 1 + d ./ C;
		Dinv = 1 + d ./ Dinv;
		y = y .* (C./Dinv);
		k = (abs(y(:)-yold(:)) > 1000*eps(y(:)));
		if ~any(k), break; end
		m = m + 1;
	end
	*/
	int maxiter= 1000;
	for (int m=1;m<=maxiter;m++)
	{
		double yold= y;
		int twom= 2*m;
		double d= m*(b-m)*x/((aminus1+twom)*(a+twom));
		C= 1+d/C;
		Dinv= 1+d/Dinv;
		y= y*(C/Dinv);
		d= -(a+m)*(aplusb+m)*x/((a+twom)*(aplus1+twom));
		C= 1+d/C;
		Dinv= 1+d/Dinv;
		y= y*(C/Dinv);
		if (Stat_abs(y-yold)<=1000*Stat_eps*y)
			return y;
	}
	/*
	if m > maxiter
		y(k) = NaN;
	end
	*/
	return Stat_GetNaN_double();
}

double Stat_expm1(double x)
{
	// Compute exp(x)-1 accurately (due to William Kahan)
	double u= exp(x);
	if (u == 1.)
		return x;
	if (u/1. == -1.)
		return -1.;
	return (u-1.)*x/log(u);
}

double Stat_log1p(double x)
{
	// Compute log(1+x) accurately (due to William Kahan)
	double u= 1.+x;
	if (u == 1.)
		return x;
	else
		return log(u)*x/(u-1.);
}

double Stat_Betainc(double x, double a, double b)
{
	// Incomplete beta function

	/*
	function y = betainc(x,a,b)
	%BETAINC Incomplete beta function.
	%   Y = BETAINC(X,Z,W) computes the incomplete beta function for
	%   corresponding elements of X, Z, and W.  The elements of X must be in
	%   the closed interval [0,1], and those of Z and W must be nonnegative.
	%   X, Z, and W must all be real and the same size (or any of them can be
	%   scalar).
	%
	%   The incomplete beta function is defined as
	%
	%     I_x(z,b) = 1./BETA(z,w) .*
	%                 integral from 0 to x of t.^(z-1) .* (1-t).^(w-1) dt
	%
	%   To compute the upper tail of the incomplete beta function, use
	%
	%     1 - BETAINC(X,Z,W) = BETAINC(1-X,W,Z).
	%
	%   If either Z or W is very large, BETAINC uses an approximation whose
	%   absolute accuracy is at least 5e-3 if Z+W > 6.
	%
	%   Class support for inputs X,Z,W:
	%      float: double, single
	%
	%   See also BETA, BETALN.

	%   Ref: Abramowitz & Stegun, Handbook of Mathematical Functions, sec. 26.5,
	%   especially 26.5.8, 26.5.20 and 26.5.21.

	%   Copyright 1984-2006 The MathWorks, Inc.
	%   $Revision: 5.16.4.9 $  $Date: 2006/04/03 17:11:13 $

	if nargin < 3
		error('MATLAB:betainc:NotEnoughInputs','Requires three input arguments.')
	elseif any(x(:) < 0 | x(:) > 1 | isnan(x(:))) || ~isreal(x)
		error('MATLAB:betainc:XoutOfRange','X must be in the interval [0,1].')
	elseif any(a(:) < 0 | isnan(a(:))) || ~isreal(a)
		error('MATLAB:betainc:PositiveZ','Z must be real and nonnegative.')
	elseif any(b(:) < 0 | isnan(b(:))) || ~isreal(b)
		error('MATLAB:betainc:PositiveW','W must be real and nonnegative.')
	end

	try
		% Preallocate y (using the size rules for plus)
		y = x + a + b;
	catch
		error('MATLAB:betainc:XZWsizeMismatch', ...
			  'X, Z and W must all the same size (or any of them can be scalar).')
	end
	% Initialize y(x==0) to 0, y(x==1) to 1. Everything else will be filled in.
	y(:) = (x==1);
	*/

	if (x<0.0 || x>1.0 || a<0.0 || b<0.0)
		return Stat_GetNaN_double();
	if (x==0)
		return 0;
	if (x==1)
		return 1;

	/*
	% Use the continued fraction unless either parameter is very large.
	approx = (a+b) > 1e7;
	*/
	bool approx= (a+b) > 1e7;
	double y= Stat_GetNaN_double();

	/*
	k = (0 < x & x < (a+1) ./ (a+b+2) & ~approx);
	if any(k(:))
		if isscalar(x), xk = x; else xk = x(k); end
		if isscalar(a), ak = a; else ak = a(k); end
		if isscalar(b), bk = b; else bk = b(k); end
		% This is x^a * (1-x)^b / (a*beta(a,b)), computed so that a==0 works.
		btk = exp(gammaln(ak+bk) - gammaln(ak+1) - gammaln(bk) + ...
									 ak.*log(xk) + bk.*log1p(-xk));
		y(k) = btk .* betacore(xk,ak,bk);
	end
	*/
	if (!approx && x>0 && x<(a+1.)/(a+b+2.))
	{
		double xk= x;
		double ak= a;
		double bk= b;
		double btk= exp(Stat_Gammaln(ak+bk)- Stat_Gammaln(ak+1) - Stat_Gammaln(bk) + ak*log(xk) + bk*Stat_log1p(-xk));
		y= btk*Stat_Betacore(xk,ak,bk);
	}

	/*
	k = ((a+1) ./ (a+b+2) <= x & x < 1 & ~approx);
	if any(k(:))
		if isscalar(x), xk = x; else xk = x(k); end
		if isscalar(a), ak = a; else ak = a(k); end
		if isscalar(b), bk = b; else bk = b(k); end
		% This is x^a * (1-x)^b / (b*beta(a,b)), computed so that b==0 works.
		btk = exp(gammaln(ak+bk) - gammaln(ak) - gammaln(bk+1) + ...
									 ak.*log(xk) + bk.*log1p(-xk));
		y(k) = 1 - btk .* betacore(1-xk,bk,ak);
	end
	*/
	if (!approx && x<1. && x>=((a+1.)/(a+b+2.)))
	{
		double xk= x;
		double ak= a;
		double bk= b;
		double btk= exp(Stat_Gammaln(ak+bk) - Stat_Gammaln(ak) - Stat_Gammaln(bk+1.) + ak*log(xk) + bk*Stat_log1p(-xk));
		y= 1.-btk*Stat_Betacore(1.-xk,bk,ak);
	}

	/*
	% NaNs may have come from a=b=0, leave those alone.  Otherwise if the
	% continued fraction in betacore failed to converge, or if we didn't use
	% it, use approximations.
	k = find((isnan(y) & (a+b>0)) | approx);
	if ~isempty(k)
		if isscalar(x), xk = x; else xk = x(k); end
		if isscalar(a), ak = a; else ak = a(k); end
		if isscalar(b), bk = b; else bk = b(k); end
		w1 = (bk.*xk).^(1/3);
		w2 = (ak.*(1-xk)).^(1/3);
		y(k) = 0.5*erfc(-3/sqrt(2)*((1-1./(9*bk)).*w1-(1-1./(9*ak)).*w2)./ ...
			   sqrt(w1.^2./bk+w2.^2./ak));

		k1 = find((ak+bk-1).*(1-xk) <= 0.8);
		if ~isempty(k1)
			if isscalar(x), xk = x; else xk = xk(k1); end
			if isscalar(a), ak = a; else ak = ak(k1); end
			if isscalar(b), bk = b; else bk = bk(k1); end
			s = 0.5*((ak+bk-1).*(3-xk)-(bk-1)).*(1-xk);
			y(k(k1)) = gammainc(s,bk,'upper');
		end
	end
	*/
	if ((IsNaN_double(y) && (a+b>0)) || approx)
	{
		double xk= x;
		double ak= a;
		double bk= b;
		double w1= exp(log(bk*xk)/3.); // = (bk*xk)^(1/3);
		double w2= exp(log(ak*(1-xk))/3.); // = (ak*(1-xk))^(1/3);
		y= 0.5*Stat_erfc(-3./sqrt(2.)*((1.-1.0/(9*bk))*w1-(1.-1.0/(9*ak))*w2)/sqrt(w1*w1/bk+w2*w2/ak));

		if ((ak+bk-1)*(1-xk)<0.8)
		{
			double xk= x;
			double ak= a;
			double bk= b;
			double s = 0.5*((ak+bk-1.0)*(3.0-xk)-(bk-1.0))*(1.0-xk);
			y= Stat_Gammainc(s,bk,false);
		}
	}

	return y;
}

// ******************************************************************
// Routines for retrieving common cdf:s

double Stat_Normcdf(double z)
{
	return 0.5 * Stat_erfc(-z / sqrt(2.0));
}

double Stat_Tcdf(double z, double df)
{
	if (df<Stat_eps)
		return Stat_GetNaN_double();
	if (df==1.0)
		return 0.5 + atan(z)/pi; // Cauchy distribution
	const double normcutoff= 1e7;
	if (df>normcutoff)
		return Stat_Normcdf(z); // Normal approximation
	double b= Stat_Betainc(df/(df+z*z), df/2, 0.5)/2;
	return b;
	// return z>=0.0 ? 1.-b : b;
}

// ******************************************************************
// Routines for convolving histograms

void Stat_HistNormalize(CVector<double> &h)
{
	// Normalize histogram to sum 1
	int n= h.GetSize();
	double *ph= h.GetBuffer(n); // optimization
	double s= 0;
	for (int i=0;i<n;i++)
		s += ph[i];
//	double s_inv= 1.0/s; // optimization: div --> mul
	for (int i=0;i<n;i++)
		ph[i] *= s;
}

bool Stat_HistConvolve(const CVector<double> &h0, double min0, double max0, 
					   const CVector<double> &h1, float min1, double max1, 
					   CVector<double> &h_out, double &min_out, double &max_out)
{
	// h0, h1 histograms, equally spaced between their respective min and max, normalized to sum 1.
	const int n= h0.GetSize();
	if (!n || n!=h1.GetSize())
		return false;
	h_out.SetSize(n);
	h_out.Fill(0);
	min_out= min0 + min1;
	max_out= max0 + max1;

	const int nm1= n-1;
	for (int i0=0;i0<n;i0++)
	{
		const double x0= double(i0)/nm1*(max0-min0)+min0;
		const double n0= h0[i0];
		for (int i1=0;i1<n;i1++)
		{
			double x1= double(i1)/nm1*(max1-min1)+min1;
			double x_out= x0+x1;
			int i_out= int((x_out-min_out)/(max_out-min_out)*nm1 + 0.5);
			h_out.SetAt(i_out, h_out[i_out] + n0*h1[i1]);
		}
	}
	return true;
}

double Stat_GetInnerProduct(const double *aData0, const double *aData1, const int n)
{
	double s=0;
	for (int i=n;--i>=0;)
		s += aData0[i]*aData1[i];
	return s;
}

double Stat_InvErf(double x)
{
	// From "A handy approximation of the error function and its inverse" (Eq. 7)
	// Sergei Winitzki, preprint, Feb 6 2008

	const double a= 0.147;
	const double a_inv= 1.0/a;
	const double b= 2.0/pi/a;

	double y1= log(1.0-x*x);
	double y2= b+y1/2.0;

	double y_out= sqrt(- b - y1/2.0 + sqrt( y2*y2 -  a_inv*y1 ) );
	if (x>=0)
		return y_out;
	else
		return -y_out;
}

double Stat_InvNormcdf(double p)
{
	double q= Stat_InvErf(2*p-1);
	if (IsN_double(q))
		return sqrt(2.0)*q;
	else
		return q;
}

bool Stat_NormalizeMatrix(CMatrix<double> &m, int mode, int cmin, int cmax)
{
	const int n= m.GetHeight();
	CVector<double> icdf(n);

	if (mode==NORM_Rank)
	{
		for (int i=0;i<n;i++)
			icdf.SetAt(i, double(i)/n); // relative rank
	}
	else if (mode==NORM_Normal)
	{
		for (int i=0;i<n;i++)
			icdf.SetAt(i, Stat_InvNormcdf(double(i+1)/(n+1)));
	}
	else if (mode==NORM_Lognormal)
	{
		for (int i=0;i<n;i++)
			icdf.SetAt(i, exp(Stat_InvNormcdf(double(i+1)/(n+1))));
	}
	else if (mode==NORM_None)
		return true;
	else
		return false;

	CVector<Stat_SortItem_double> vSI(n);
	Stat_SortItem_double *pSI= vSI.GetBuffer();
	for (int c=0;c<m.GetWidth();c++)
	{
		if (c>=cmin && c<=cmax)
		{
			for (int i=0;i<n;i++)
			{
				pSI[i].m_Value= m.GetAt(i,c);
				pSI[i].m_nIndex= i;
			}
		
			Stat_QuicksortIndexed_double(pSI, vSI.GetSize(), true);

			for (int i=0;i<n;i++)
				m.SetAt(pSI[i].m_nIndex, c, icdf[i]);
		}
	}

	return true;
}

bool Stat_NormalizeMatrix(CMatrix<double> &m, int mode)
{
	return Stat_NormalizeMatrix(m, mode, 0, m.GetWidth());
}

#ifdef _WIN32

//// 
// Acklam

#define M_PI (3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679)
#define M_SQRT2 (1.41421356237309504880)
#define M_SQRT_2 (0.7071067811865475244008443621048490392848359376887)
#define M_SQRT2PI (2.50662827463100050242)
#define M_1_SQRTPI (0.564189583547756286948)
#define M_SQRTPI (1.77245385090551602792981)
#define M_1_SQRT2PI (M_SQRT_2*M_1_SQRTPI)

/*
 * A normally distributed random number generator.  We avoid
 * the uniform rv's being 0.0 since this will result in infinte
 * values, and double count the 0 == 2pi.
 */
double Acklam_random_normal() {
 static int i = 1;
 static double u[2] = {0.0, 0.0};
 register double r[2];

 if (i == 1) {
  r[0] = sqrt(-2*log((double)(rand()+1)/(double)(RAND_MAX+1)));
  r[1] = 2*M_PI*(double)(rand()+1)/(double)(RAND_MAX+1);
  u[0] = r[0]*sin(r[1]);
  u[1] = r[0]*cos(r[1]);
  i = 0;
 } else {
  i = 1;
 }

 return u[i];
};

/*
 * The standard normal PDF, for one random variable.
 */
inline double Acklam_stdnormal_pdf(double u)
{
 return exp(-u*u/2)/M_SQRT2PI;
};

/*
 * An implementation of adaptive, recursive Newton-Cotes integration.
 * Based on the MATLAB implementation, but covered in a lot of books...
 *
 * This only does integration over the standard normal PDF.  It's just
 * here to check the error function approximations.
 */
#define LEVMAX 10
double Acklam_quad8_stdnormal_pdf(double a, double b, double Q = 1.0)
{
 /* The magic Newton-Cotes weights */
 const int w[9] = {3956, 23552, -3712, 41984, -18160, 41984, -3712, 23552,
3956};
 const int dw = 14175;
 static int level = -1;
 static double tol = 1e-30;
 register double h, Q1 = 0.0, Q2 = 0.0;
 register int i;

 level++;
 h = (b-a)/16.0;
 for (i = 0; i < 9; i++) {
  Q1 += h*w[i]*Acklam_stdnormal_pdf(a+i*h)/dw;
  Q2 += h*w[i]*Acklam_stdnormal_pdf(a+(i+8)*h)/dw;
 };
 /* This is the adaptive recursive bit.  We only recurse if we can
improve... */
 if (fabs(Q1+Q2-Q) > tol*fabs(Q1+Q2) && level <= LEVMAX) {
  tol = tol/2;
  Q1 = Acklam_quad8_stdnormal_pdf(a,(a+b)/2,Q1);
  Q2 = Acklam_quad8_stdnormal_pdf((a+b)/2,b,Q2);
  tol = tol*2;
 }
 level--;
 return Q1 + Q2;
}

/*
 * The standard normal CDF, for one random variable.
 *
 *   Author:  W. J. Cody
 *   URL:   http://www.netlib.org/specfun/erf
 *
 * This is the erfc() routine only, adapted by the
 * transform stdnormal_cdf(u)=(erfc(-u/sqrt(2))/2;
 */
double Acklam_stdnormal_cdf(double u)
{
 const double a[5] = {
  1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
  1.887426188426510e+002,3.209377589138469e+003
 };
 const double b[5] = {
  1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
  1.813893686502485e+003,8.044716608901563e+003
 };
 const double c[9] = {
  2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00,
  6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02,
  1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03
 };
 const double d[9] = {
  1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02,
  5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03,
  4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03
 };
 const double p[6] = {
  1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1,
  1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4
 };
 const double q[6] = {
  1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00,
  5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3
 };
 register double y, z;

 if (IsNaN_double(u))
  return Stat_GetNaN_double();
 if (!_finite(u))
  return (u < 0 ? 0.0 : 1.0);
 y = fabs(u);
    if (y <= 0.46875*M_SQRT2) {
  /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
  z = y*y;
  y = u*((((a[0]*z+a[1])*z+a[2])*z+a[3])*z+a[4])
       /((((b[0]*z+b[1])*z+b[2])*z+b[3])*z+b[4]);
  return 0.5+y;
 }
 z = exp(-y*y/2)/2;
 if (y <= 4.0) {
  /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
  y = y/M_SQRT2;
  y =
((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8])


/((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]);

  y = z*y;
    } else {
  /* evaluate erfc() for |u| > sqrt(2)*4.0 */
  z = z*M_SQRT2/y;
  y = 2/(y*y);
        y = y*(((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5])
    /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]);
        y = z*(M_1_SQRTPI-y);
    }
 return (u < 0.0 ? y : 1-y);
};

/*
 * The inverse standard normal distribution.
 *
 *   Author:      Peter John Acklam <pjacklam@online.no>
 *   URL:         http://home.online.no/~pjacklam
 *
 * This function is based on the MATLAB code from the address above,
 * translated to C, and adapted for our purposes.
 */
double Acklam_stdnormal_inv(double p)
{
 const double a[6] = {
  -3.969683028665376e+01,  2.209460984245205e+02,
  -2.759285104469687e+02,  1.383577518672690e+02,
  -3.066479806614716e+01,  2.506628277459239e+00
 };
 const double b[5] = {
  -5.447609879822406e+01,  1.615858368580409e+02,
  -1.556989798598866e+02,  6.680131188771972e+01,
  -1.328068155288572e+01
 };
 const double c[6] = {
  -7.784894002430293e-03, -3.223964580411365e-01,
  -2.400758277161838e+00, -2.549732539343734e+00,
   4.374664141464968e+00,  2.938163982698783e+00
 };
 const double d[4] = {
   7.784695709041462e-03,  3.224671290700398e-01,
   2.445134137142996e+00,  3.754408661907416e+00
 };

 register double q, t, u;

 if (IsNaN_double(p) || p > 1.0 || p < 0.0)
  return Stat_GetNaN_double();
 if (p == 0.0)
  return Stat_GetNegInfty_double(); 
 if (p == 1.0)
  return Stat_GetPosInfty_double();
 q = __min(p,1-p);
 if (q > 0.02425) {
  /* Rational approximation for central region. */
  u = q-0.5;
  t = u*u;
  u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
    /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
 } else {
  /* Rational approximation for tail region. */
  t = sqrt(-2*log(q));
  u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
   /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
 }
 /* The relative error of the approximation has absolute value less
    than 1.15e-9.  One iteration of Halley's rational method (third
    order) gives full machine precision... */
 t = Acklam_stdnormal_cdf(u)-q;    /* error */
 t = t*M_SQRT2PI*exp(u*u/2);   /* f(u)/df(u) */
 u = u-t/(1+u*t/2);     /* Halley's method */

 return (p > 0.5 ? -u : u);
};

#endif // _WIN32 Acklam