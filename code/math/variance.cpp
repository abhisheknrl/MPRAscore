/*******************************************************************************
 *
 *   variance.cpp -- Implements several shrinkage-based variance estimators
 *                   frequently used in the analysis of microarray data:
 *
 *   Björn Nilsson, 2006 
 *
 */

#include "types/vector.h"
#include "math/statfun.h"

/******************************************************************************/
// Helper functions

double Get_Trigamma(double x)
// 
// Computes the Trigamma funtion for x>0.
// Adopted from PROB.C by John Burkhardt.
//
{
	if (x<=0.0)
	{
		ASSERT(false);
		return 0.0;
	}

	const double a = 0.0001;
	const double b = 5.0;
	const double b2 =   1.0 / 6.0;
	const double b4 = - 1.0 / 30.0;
	const double b6 =   1.0 / 42.0;
	const double b8 = - 1.0 / 30.0;
	double y;
	double z;

	// If X is smaller than A, use a small value approximation.
	if (x<=a)
		return 1.0 / x / x;
  
	// Otherwise, increase the argument to B <= ( X + I ).
	z = x;
	double value = 0.0;
    while ( z < b )
	{
		// The recurrence relation for polygamma functions yields
		value = value + 1.0 / z / z;
		z = z + 1.0;
	}
	
	// ...and then apply an asymptotic formula.
	y = 1.0 / z / z;

	value = value + 0.5 * 
		y + ( 1.0 
		+ y * ( b2 
		+ y * ( b4 
		+ y * ( b6 
		+ y *   b8 )))) / z;
 
	return value;
}

double Get_PolynomialValue(double x, const double a[], int n)
{
	double s= a[0];
	for (int i=1;i<=n;i++)
	{
		s += a[i]*x;
		x *= x;
	}
	return s;
}

double Get_Digamma(double x) 
{ 
	// 
	// Computes Digamma (psi) function for x>0.
	// Adopted from PROB.C by John Burkhardt.
	//

	if (x<=0)
	{
		ASSERT(false);
		return 0.0;
	}

	// Euler-Mascheroni constant
	const double EM= 0.577215664901532860606512090082402431042;

	// Interpolation polynomial coefficients
	const double A[] = 
	{ 
		8.33333333333333333333E-2, 
		-2.10927960927960927961E-2, 
		7.57575757575757575758E-3, 
		-4.16666666666666666667E-3, 
		3.96825396825396825397E-3, 
	   -8.33333333333333333333E-3, 
		8.33333333333333333333E-2 
	};

	double s, w, y, z; 
	int i, n; 
	double nz = 0.0; 

	// Check for positive integer up to 10.
	if( (x <= 10.0) && (x == floor(x)) ) 
	{ 
	   y = 0.0; 
	   n = (int)x; 
	   for( i=1; i<n; i++ ) 
	   { 
		 w = i; 
		 y += 1.0/w; 
	   } 
	   y -= EM; 

	   return y; 
	} 

	s = x; 
	w = 0.0; 
	while(s<10.0) 
	{ 
	   w += 1.0/s; 
	   s += 1.0; 
	} 

	if(s<1.0e17) 
	{ 
	   z = 1.0/(s * s); 
	   y = z * Get_PolynomialValue( z, A, 6 ); 
	} 
	else 
	   y = 0.0; 

	y = log(s) -  (0.5/s)  -  y  -  w; 

	return(y); 
} 

double Get_Trigamma_Inverse(double x)
{
	if (x<=0.0)
	{
		ASSERT(false);
		return 0.0;
	}

	double y= 0.5+1.0/x;
	const double eps= 0.00000001;
	double delta= -2.0*y;
	while (-delta/y>eps)
	{
		// Newton-Raphson procedure, where we use left 
		// derivatives to ascertain that we never over-
		// estimate the step length. The termination 
		// criterion is the one proposed in Smyth (2004).
		double g3= Get_Trigamma(y);
		double dy= 0.001*y;
		double g3_deriv= (g3-Get_Trigamma(y-dy))/dy;

		delta= g3*(1.0-g3/x)/g3_deriv;
		y += delta;
	}
	return y;
}

/******************************************************************************/
// Exported functions

bool GetVariances_Standard(CVector<double> vVarIn, CVector<double> &vVarOut)
{
	// Dummy function to allow for a uniform code look
	vVarOut= vVarIn;
	return true;
}

bool GetVariances_Unit(CVector<double> vVarIn, CVector<double> &vVarOut)
{
	// Simply sets all variances to one.

	int N= vVarIn.GetSize();
	vVarOut.SetSize(N);
	for (int i=0;i<N;i++)
		vVarOut.SetAt(i, 1.0);

	return true;
}

bool GetVariances_Equal(CVector<double> vVarIn, CVector<double> &vVarOut)
{
	//
	// Implement the naïve shrunken variance estimator, i.e. all 
	// genes are assumed to have identical variance, estimated
	// as the mean variance across all genes.
	//

	int N= vVarIn.GetSize();
	vVarOut.SetSize(N);
	if (N>0)
	{
		double mu= 0;
		for (int i=0;i<N;i++)
			mu += vVarIn[i];
		mu /= N;

		for (int i=0;i<N;i++)
			vVarOut.SetAt(i, mu);
	}

	return true;
}

bool GetVariances_Tusher(CVector<double> vVarIn, double q, CVector<double> &vVarOut, double &sigma_0)
{
	//
	// Implements the SAM (Statistical Analysis of Microarray) 
	// variance estimator, proposed in: Tusher, V.G., Tibshirani, R., 
	// and Chu, G. (2001). "Statistical analysis of microarrays 
	// applied to the ionizing radiation response", Proc. Natl. 
	// Acad. Sci. U S A 98, 5116-5121.
	//
	// The additional parameter q indicates the variance quantile 
	// used for selecting the fudge factor. A common choice
	// for this parameter is 0.9 (a.k.a., Efron's regularization, 
	// see also GetVariances_Efron). The additional output parameter
	// sigma_0 will contain the value of the fudge factor itself.
	// Note that the "variances" output by SAM are not variances
	// in the traditional sense, implying that z-scores computed 
	// using these values will not be t-statistics.
	//

	if (q<0.0 || q>1.0)
		return false;
	int N= vVarIn.GetSize();
	vVarOut.SetSize(N);
	if (N==0)
		return true;

	// Compute fudge factor
	double *pV= vVarIn.GetBuffer(N);
	Stat_Quicksort_double(pV, N, true);
	sigma_0= sqrt(Stat_GetQuantile_double(pV, N, q));

	for (int i=0;i<N;i++)
	{
		double s= sqrt(vVarIn[i]) + sigma_0;
		vVarOut.SetAt(i, s*s);
	}

	return true;
}

bool GetVariances_Efron(CVector<double> vVarIn, CVector<double> &vVarOut)
{
	double sigma_0_dummy;
	return GetVariances_Tusher(vVarIn, 0.5, vVarOut, sigma_0_dummy);
}

bool GetVariances_Cui(CVector<double> vVarIn, CVector<double> &vVarOut)
{
	//
	// Implements the shrunken variance estimator proposed in
	// Cui, X. et al (2005) in "Improved statistical tests for differential
	// gene expression by shrinking variance components estimates",
	// Biostatistics 6(1), pp 59-75.
	//

	int N= vVarIn.GetSize();
	vVarOut.SetSize(N);
	if (N==0)
		return true;

	return true;
}

bool GetVariances_Smyth(CVector<double> vVarIn, CVector<double> vDF, CVector<double> &vVarOut, double &s0sqr, double &d0)
{
	// 
	// Implements the moderated variance estimator proposed by
	// Gordon Smyth in "Linear Models and Empirical Bayes Methods 
	// for Assessing Differential Expression in Microarray 
	// Experiments", Statistical Applications in Genetics
	// and Molecular Biology (2004), 3(1), Article 3.
	//
	// The additional vector vDF represents the number of degrees
	// of freedom used (zeros allowed) when computing each element 
	// in vVarIn. The length of vDF must match that of vVarIn. 
	//

	int N= vVarIn.GetSize();
	if (vDF.GetSize()!=N)
		return false;
	vVarOut.SetSize(N);
	if (N==0)
		return true;
	if (N<100) // Very small sample => non-intended use. Use standard variance estimator.
		return GetVariances_Standard(vVarIn, vVarOut);

	CVector<double> vE;
	vE.SetSize(N);
	for (int i=0;i<N;i++)
	{
		if (vDF[i]>0)
		{
			// Positive degrees-of-freedom
			if (vVarIn[i]<=0.0)
				return false; // Spurious negative variance detected. Fail.

			double df_half= vDF[i]/2.0;
			vE.SetAt(i, log(vVarIn[i]) - Get_Digamma(df_half) + log(df_half));
		}
		else if (vDF[i]==0)
			// Zero degrees-of-freedom
			vE.SetAt(i, 0.0);
		else
		{
			printf("Neg!");
			// Negative degrees-of-freedom detected. Fail.
			return false;
		}
	}

	double e_mean, e_var;
	Stat_GetMeanAndVariance_double(vE, vE.GetSize(), e_mean, e_var);
	double s= 0;
	for (int i=0;i<N;i++)
		s += Get_Trigamma(vDF[i]/2.0);

	double RHS= e_var-s/N;
	if (RHS>0)
	{
		// Solve Equation 1 (Smyth 2004, page 11) to obtain d0.
		d0= 2.0*Get_Trigamma_Inverse(RHS);

		// Solve for s0, given d0.
		s0sqr= exp(e_mean + Get_Digamma(d0/2.0) - log(d0/2.0));

		// Compute moderated variances
		double tmp= d0*s0sqr;
		for (int i=0;i<N;i++)
			vVarOut.SetAt(i, (tmp + vDF[i]*vVarIn[i])/(d0 + vDF[i]));
	}
	else
	{
		// No evidence that the underlying variances vary between genes. 
		d0= 1000000000000000000000.0; // infinity

		// Compute common moderated variance
		s0sqr= exp(e_mean);
		for (int i=0;i<N;i++)
			vVarOut.SetAt(i, s0sqr);
	}

	return true;
}

bool GetVariances_Smyth(CVector<double> vVarIn, CVector<double> vDF, CVector<double> &vVarOut)
{
	//
	// Simplfied call without hyperparameter output.
	//
	double s0sqr, d0;
	return GetVariances_Smyth(vVarIn, vDF, vVarOut, s0sqr, d0);
}

bool GetVariances_Smyth(CVector<double> vVarIn, double nDF, CVector<double> &vVarOut)
{
	//
	// Simplified call with equal degrees of freedom 
	// for all genes and no hyperparameter output, i.e.
	// when there is no missing data.
	//

	CVector<double> vDF;
	vDF.SetSize(vVarIn.GetSize());
	for (int i=0;i<vVarIn.GetSize();i++)
		vDF.SetAt(i, nDF);
	return GetVariances_Smyth(vVarIn, vDF, vVarOut);
}
