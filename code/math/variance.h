/*******************************************************************************
 *
 *   variance.h --  Implements commonly used shrinkage-based variance 
 *					estimators for the analysis of microarray data.
 *
 *   Björn Nilsson, 2006 
 *
 */

bool GetVariances_Unit(CVector<double> vVarIn, CVector<double> &vVarOut);
bool GetVariances_Equal(CVector<double> vVarIn, CVector<double> &vVarOut);
bool GetVariances_Standard(CVector<double> vVarIn, CVector<double> &vVarOut); // Dummy function
bool GetVariances_Tusher(CVector<double> vVarIn, double q, CVector<double> &vVarOut, double &sigma_0);
bool GetVariances_Efron(CVector<double> vVarIn, CVector<double> &vVarOut); // Simplified call
bool GetVariances_Cui(CVector<double> vVarIn, CVector<double> &vVarOut);
bool GetVariances_Smyth(CVector<double> vVarIn, CVector<double> vDF, CVector<double> &vVarOut, double &s0sqr, double &d0);
bool GetVariances_Smyth(CVector<double> vVarIn, CVector<double> vDF, CVector<double> &vVarOut);
bool GetVariances_Smyth(CVector<double> vVarIn, double nDF, CVector<double> &vVarOut); // Simplified call
