/*
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
an improved signal by accounting for experimentally observed patterns 
in probe behavior and handling error at the appropriately at low and 
high signal values.

Copyright (C) 2004 Affymetrix, Inc.

This program is free software; you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by 
the Free Software Foundation; either version 2 of the License, 
or (at your option) any later version.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details.

You should have received a copy of the GNU General Public License 
along with this program;if not, write to the 

Free Software Foundation, Inc., 
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

#include <malloc.h>
#include <memory>
#include <new>
#include <math.h>
#include "plieralg.h"
#include "affyheapsort.h"
#include "error.h"

//////////////////////////////////////////////////////////////////////

enum NormErrorCodes {
	NoError,
	MemoryError,
	NumErrors
};

//////////////////////////////////////////////////////////////////////

static int g_ErrorCode=NoError;

//////////////////////////////////////////////////////////////////////
#define FPMIN         0.00000001f
#define EPS           1e-2
#define IQR_DIVIDER   1.34f

//////////////////////////////////////////////////////////////////////

float logtwo(float value)
{
	if (value>FPMIN)
		return (float) (log(value)/log(2));
	else
		return(-50);
}

//////////////////////////////////////////////////////////////////////

float LN(float value, float zerolevel)
{
	if (value>EPS)
		return (float) logtwo(value);
	else
		return zerolevel;
}

//////////////////////////////////////////////////////////////////////

float Exp(float value)
{
	return (float) pow(2,value);
}

//////////////////////////////////////////////////////////////////////

inline
void TransferVector(double *A, double *B, long Length)
{
	for (long i=0; i<Length; i++)
		A[i] = B[i];
}

//////////////////////////////////////////////////////////////////////
inline
void ScrambleTransferVector(double *A, double *B, long *Twist, long Length)
{
	long i;
	for (i=0; i<Length; i++)
		A[i] = B[Twist[i]];
}

//////////////////////////////////////////////////////////////////////
inline
void InitializeVector(double *v, long n, double fValue)
{
	for (long i=0; i<n; i++)
		v[i] = fValue;
}

//////////////////////////////////////////////////////////////////////

void
LogVector(double *Vector, long Length)
{
	long i;
	for (i=0; i<Length; i++)
		Vector[i] = log(Vector[i]);
}

//////////////////////////////////////////////////////////////////////

void
ExpVector(double *Vector, long Length)
{
	long i;
	for (i=0; i<Length; i++)
		Vector[i] = exp(Vector[i]);
}

//////////////////////////////////////////////////////////////////////
inline
void ShrinkVector(double *Vector, double *ShiftVector, long Length, double dropmax)
{
	long i;
	for (i=0; i<Length; i++)
		if (Vector[i]+ShiftVector[i] < Vector[i]/dropmax)
			ShiftVector[i] = Vector[i]/dropmax - Vector[i];
}

//////////////////////////////////////////////////////////////////////

void StepVector(double *Vector, double *ShiftVector, long Length, double lambda)
{
	long i;
	for (i=0; i<Length; i++)
		Vector[i] += lambda*ShiftVector[i];
}

//////////////////////////////////////////////////////////////////////

inline
void AugmentData(double **Matrix, long NumExperiments, long NumProbes, double Safety)
{
	long i,j;
	// ensure data is positive by adding small constant everywhere
	for (i=0; i<NumExperiments; i++)
		for (j=0; j<NumProbes; j++)
			Matrix[i][j] += Safety; // small positive number
}

//////////////////////////////////////////////////////////////////////

void AugmentFloatData(float* data, long nSize, float safety)
{
	for (long i=0; i<nSize; i++)
		data[i] += safety;
}

//////////////////////////////////////////////////////////////////////

void AugmentDoubleData(double* data, long nSize, double safety)
{
	for (long i=0; i<nSize; i++)
		data[i] += safety;
}

//////////////////////////////////////////////////////////////////////

long DoHeapSeaRaw(plier_data *pData, double *Concentration, double *Affinity,
			   double **TestDataTable, double ssqlimit, bool bNoAffinity)
{
	long iter;
	long i,j, index, ti;
	long medxval, medxvallo, medpval, medpvallo;
	double ssq;


	// do median polish by heapindex
	double median;
	double *Val;
	long *Rank;
	long ValSize;

	long nSize = pData->m_nAnalyses*pData->m_nProbes;

	Val = new double [nSize]; //scratch space
	if (Val == 0)
		return NO_DATAMEM;

	Rank = new long [nSize]; // scratch space for rank
	if (Rank == 0)
	{
		delete[] Val;
		return NO_DATAMEM;
	}

	ssq = 10;
	long iterations = pData->m_algParams->seaiteration;
	for (iter=0; iter<iterations && ssq>ssqlimit; iter++)
	{
		ssq = 0;
		for (i=0; i<pData->m_nAnalyses;)
		{
			index = 0;
			for (ti=i; ti<pData->m_nReplicates[i]; ti++)
			{
				// if replicates, would repeat this operation to load all values into Val
				for (j=0; j<pData->m_nProbes; j++, index++)
					Val[index] = TestDataTable[ti][j]-Affinity[j]; // copy the values
					//Val[index] = TestDataTable[ti*pData->m_nProbes + j]-Affinity[j]; // copy the values
			}
			// Valsize is number of things loaded
			ValSize = index;
			HeapIndex(Val, 0, Rank, ValSize);

			// get median, subtract from column
			medpval = (ValSize-ValSize%2)/2;
			medpvallo = (ValSize-1-(ValSize-1)%2)/2;
			median = (Val[Rank[medpval]]+Val[Rank[medpvallo]])/2;

			// if replicates, would set concentration[i] for each replicate
			for (ti=i; ti<pData->m_nReplicates[i]; ti++)
			{
				ssq += (Concentration[ti]-median)*(Concentration[ti]-median);
				Concentration[ti] = median;
			}
			i = pData->m_nReplicates[i]; //jump forwards by an appropriate amount
		}

		for (j=0; j<pData->m_nProbes && !bNoAffinity; j++) // only do this if updating affinities
		{
			// no replicate probes
			for (i=0; i<pData->m_nAnalyses; i++)
				Val[i] = TestDataTable[i][j]- Concentration[i];
				//Val[i] = TestDataTable[i*pData->m_nProbes+j]- Concentration[i];
			ValSize = pData->m_nAnalyses;
			HeapIndex(Val, 0, Rank, ValSize);
			// get median, subtract from column
			medxval = (ValSize-ValSize%2)/2;
			medxvallo = (ValSize-1-(ValSize-1)%2)/2;
			median = (Val[Rank[medxval]]+Val[Rank[medxvallo]])/2;

			ssq += (Affinity[j]-median)*(Affinity[j]-median);
			Affinity[j] = median;
		}
	}

	delete[] Val;
	delete[] Rank;

	if (iter==iterations && ssq>ssqlimit)
		return MAXIT_SEA_REACHED;
	else
		return NO_PLIER_ERROR;
}

//////////////////////////////////////////////////////////////////////

void BalanceAffinity(double *Concentration, double *Affinity,
					 long NumExperiments, long NumProbes)
{
	double total;
	long i,j;

	// balance affinities - crude - should build in at top?
	total = 0;
	for (j=0; j<NumProbes; j++)
		total += Affinity[j];
	total /= NumProbes;

	for (i=0; i<NumExperiments; i++)
		Concentration[i] *= total;
	for (j=0; j<NumProbes; j++)
		Affinity[j] /=total;
}

//////////////////////////////////////////////////////////////////////

long doSEA(plier_data *pData, double* Concentration,
		   double* Affinity, double **IHash, bool bNoAffinity)
{
	long i,j;
	double pm, mm;
	double pm_minus_mm;
	double fFourTimesAttenuation = 4*pData->m_algParams->attenuation;

	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nProbes; j++)
		{
			pm = pData->m_fPM[i][j];
			mm = pData->m_fMM[i][j];
			pm_minus_mm = pm-mm;
			IHash[i][j] = log(((pm_minus_mm)+sqrt((pm_minus_mm)*(pm_minus_mm)+fFourTimesAttenuation*pm*mm))/2);
		}

	LogVector(Concentration, pData->m_nAnalyses);
	LogVector(Affinity, pData->m_nProbes);

	long lRes = DoHeapSeaRaw(pData, Concentration, Affinity, IHash, pData->m_algParams->seaconvergence, bNoAffinity);
	if (lRes != NO_PLIER_ERROR)
		return lRes;

	ExpVector(Concentration, pData->m_nAnalyses);
	ExpVector(Affinity, pData->m_nProbes);

	BalanceAffinity(Concentration, Affinity, pData->m_nAnalyses, pData->m_nProbes);

	return NO_PLIER_ERROR;
}

//////////////////////////////////////////////////////////////////////
inline
double GenerateReferenceForX(double *X, long Length)
{
	long i;
	double total=0;
	for (i=0; i<Length; i++)
		total += log(X[i]); // logarithm of X
	total /= Length;
	total = exp(total); // geometric mean
	return(total); // geometric mean
}

//////////////////////////////////////////////////////////////////////
double JustError(double Concentration,
		   double Affinity,
		   double Hash, double PM, double MM, bool bUseMM)
{
	double a, c, y, q, r, e;

	if (bUseMM)
	{
	   a = Affinity;
      c = Concentration;
      y = a * c;
      q = sqrt (y * y + Hash);
      r = (y + q) / (2 * PM);
      e = log (r);
	}
	else
	{
      a = Affinity;
      c = Concentration;
      y = a * c;
      r = (y + MM) / PM;
      e = log (r);
	}
	return(e);
}

//////////////////////////////////////////////////////////////////////

double JustFit(double Concentration, double Affinity, double Hash, double PM, double MM, double z, bool bUseMM)
{
	double e;

	e = JustError(Concentration, Affinity, Hash, PM, MM, bUseMM);
	return(e*e/(1+e*e/z));
}

//////////////////////////////////////////////////////////////////////
inline
void PLIERMMLikelihood(double *Likelihood,
				  double *CDeriv,
				  double *ADeriv,
				  double *CGrad,
				  double *AGrad,
				  double Concentration,
				  double Affinity,
				  double Hash,
				  double PM,
				  double MM,
				  double z, bool bUseMM)
{
	double a,c,y,q,r,e,x,df,dgC,dgA,ddf, xsq, esq;

	if (bUseMM)
	{
		a = Affinity;
		c = Concentration;
		y = a*c;
		q = sqrt(y*y+Hash);
		r = (y+q)/(2*PM);
		e = log(r);
		esq = e*e;
		x = 1+(esq)/z;
		xsq = x*x;
		// f = g^2/(1+g^2/z)
		*Likelihood = esq/x;
	// give me the likelihood I need
		df = (2*e)/xsq;

		// df = 2g dg
		//df = 2*e;

		// dg = 1/r dr
		// dr = (a/q)r
		dgC = a/q;
		dgA = c/q;
		*CDeriv = df*dgC;
		*ADeriv = df*dgA;

		// d^2f = 2g/(1+g^2/z)^2 d^2g + 2[3*(1+g^2/z)-4]/(1+g^2/z)^3 (dg)^2
		ddf = 2/xsq;

		*CGrad = ddf*(dgC*dgC);
		*AGrad = ddf*(dgA*dgA);
	}
	else
	{
		a = Affinity;
		c = Concentration;
		y = a*c;
		r = (y+MM)/PM;
		e = log(r);
		esq = e*e;
		x = 1+esq/z;
		xsq = x*x;

		*Likelihood = esq/x;
		df = (2*e)/xsq;

		// 1/r dr = PM/(y+MM)*a/PM = a/(y+MM)
		dgC = a/(y+MM);
		dgA = c/(y+MM);

		*CDeriv = df*dgC;
		*ADeriv = df*dgA;

		ddf = 2/xsq;
		*CGrad = ddf*(dgC*dgC);
		*AGrad = ddf*(dgA*dgA);
	}
}

//////////////////////////////////////////////////////////////////////
inline
void RoughnessPenaltyForX(double *LocalLikelihood, double *Deriv, double *Grad, double X, double Reference, double alpha, double alpha_times_2)
{
	// some function of a concentration against a reference level
	double green;
	//double alpha_times_2=alpha*2;
	green = X/Reference;
	double loggreen = log(green);
	*LocalLikelihood = alpha*loggreen*loggreen;
	*Deriv = alpha_times_2*loggreen/X;
	*Grad = alpha_times_2/(X*X);
}

//////////////////////////////////////////////////////////////////////
inline
double UpdateLikelihoodForRoughness(double *X, double *XDeriv, double *XGrad, long Length, double penalty, double Reference)
{
	double  LocalLikelihood, Deriv, Grad, LogLikelihood;
	long i;
	LogLikelihood = 0;

	// do not compute reference here, to uncouple likelihood
	double penalty_times_2 = penalty*2;
	for (i=0; i<Length; i++)
	{
		RoughnessPenaltyForX(&LocalLikelihood,&Deriv, &Grad,X[i],Reference,penalty,penalty_times_2);
		LogLikelihood += LocalLikelihood;
		XDeriv[i] += Deriv;
		XGrad[i] += Grad;
	}
	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////
inline
double ComputeGlobalLikelihood(plier_data* pData, double *Concentration, double *Affinity,
							   double *NewtonCDeriv, double *NewtonADeriv, double *NewtonCGrad, double *NewtonAGrad,
							   double **IHash)
{
	double LogLikelihood, LocalLikelihood;
	double CDeriv, ADeriv, CGrad, AGrad;
	long i,j;
	double ReferenceAffinity, ReferenceConcentration;

	LogLikelihood = 0;
	memset(NewtonCDeriv, 0, pData->m_nAnalyses*sizeof(double));
	memset(NewtonCGrad, 0,  pData->m_nAnalyses*sizeof(double));
	memset(NewtonADeriv, 0, pData->m_nProbes*sizeof(double));
	memset(NewtonAGrad, 0, pData->m_nProbes*sizeof(double));

	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nProbes; j++)
		{
			PLIERMMLikelihood(&LocalLikelihood, &CDeriv, &ADeriv, &CGrad, &AGrad,
				Concentration[i], Affinity[j], IHash[i][j],
				pData->m_fPM[i][j], pData->m_fMM[i][j],
				pData->m_algParams->gmcutoff, pData->m_algParams->usemm);
			LogLikelihood += LocalLikelihood;
			NewtonCDeriv[i] += CDeriv;
			NewtonADeriv[j] += ADeriv;
			NewtonCGrad[i] += CGrad;
			NewtonAGrad[j] += AGrad;
		}

	ReferenceAffinity = GenerateReferenceForX(Affinity, pData->m_nProbes);
	ReferenceConcentration = GenerateReferenceForX(Concentration, pData->m_nAnalyses);

	LogLikelihood += UpdateLikelihoodForRoughness(Affinity, NewtonADeriv, NewtonAGrad,
						pData->m_nProbes, pData->m_algParams->probepenalty, ReferenceAffinity);
	LogLikelihood += UpdateLikelihoodForRoughness(Concentration, NewtonCDeriv, NewtonCGrad,
						pData->m_nAnalyses, pData->m_algParams->concpenalty, ReferenceConcentration);

	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////

void Join_Replicates(long *Replicate, double *NewtonCDeriv, double *NewtonCGrad, long NumExperiments)
{
	long i, ti, NumRep;
	double CDeriv;
	double CGrad;

	for (i=0; i<NumExperiments;)
	{
		CDeriv = CGrad = 0;
		for (ti=i; ti<Replicate[i]; ti++)
		{
			CDeriv += NewtonCDeriv[ti];
			CGrad += NewtonCGrad[ti];
		}
		NumRep = Replicate[i]-i;
		CDeriv /= NumRep;
		CGrad /= NumRep;
		for (ti=i; ti<Replicate[i]; ti++)
		{
			NewtonCDeriv[ti] = CDeriv;
			NewtonCGrad[ti] = CGrad;
		}
		i = Replicate[i];
	}
}

//////////////////////////////////////////////////////////////////////

double ComputeExperimentLogLikelihood(plier_data *pData,
									  double *Concentration, double *Affinity,
									  double *NewtonCDeriv, double *NewtonCGrad,
									  double **IHash, long WhichExp)
{
	double ReferenceConcentration;
	double CDeriv, ADeriv, CGrad, AGrad, LocalLikelihood;
	double LogLikelihood;
	long j, ti;


	ReferenceConcentration = GenerateReferenceForX(Concentration, pData->m_nAnalyses);

	// for each experimental replicate
	LogLikelihood = 0;
	for (ti=WhichExp; ti<pData->m_nReplicates[WhichExp]; ti++)
	{
		for (j=0; j<pData->m_nProbes; j++)
		{
			PLIERMMLikelihood(&LocalLikelihood, &CDeriv, &ADeriv, &CGrad, &AGrad,
				Concentration[ti], Affinity[j], IHash[ti][j],
				pData->m_fPM[ti][j], pData->m_fMM[ti][j],
				pData->m_algParams->gmcutoff, pData->m_algParams->usemm);
			LogLikelihood += LocalLikelihood;
		}
	}
	// roughness penalty over all experiments
	LogLikelihood += UpdateLikelihoodForRoughness(Concentration, NewtonCDeriv, NewtonCGrad,
						pData->m_nAnalyses, pData->m_algParams->concpenalty, ReferenceConcentration);
	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////

double ComputeProbeLogLikelihood(plier_data* pData,
								 double *Concentration, double *Affinity,
								 double *NewtonADeriv, double *NewtonAGrad,
								 double **IHash, long WhichProbe)
{
	double ReferenceAffinity;
	double CDeriv, ADeriv, CGrad, AGrad, LocalLikelihood;
	double LogLikelihood;
	long i;

	ReferenceAffinity = GenerateReferenceForX(Affinity, pData->m_nProbes);
	// for each probe generate default values
	LogLikelihood = 0;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		PLIERMMLikelihood(&LocalLikelihood, &CDeriv, &ADeriv, &CGrad, &AGrad,
				Concentration[i], Affinity[WhichProbe], IHash[i][WhichProbe],
				pData->m_fPM[i][WhichProbe], pData->m_fMM[i][WhichProbe],
				pData->m_algParams->gmcutoff, pData->m_algParams->usemm);
		LogLikelihood += LocalLikelihood;
	}
	LogLikelihood += UpdateLikelihoodForRoughness(Affinity, NewtonADeriv, NewtonAGrad,
						pData->m_nProbes, pData->m_algParams->probepenalty, ReferenceAffinity);
	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////

long SearchGridOptimum(plier_data* pData, double *Concentration, double *Affinity,
					   double *NewtonCDeriv, double *NewtonCGrad, double *NewtonADeriv, double *NewtonAGrad,
					   double **IHash,
					   double epsilon, bool bNoAffinity)
{
	long i,j, ti, tti;
	double trialval, oldval;
	double oldLogLikelihood;
	double LogLikelihood;
	long converged;

	converged = 1; // assume at optimum
	// search plausible grid for better possible optima
	// obvious grid placement is determined by concentration values
	for (i=0; i<pData->m_nAnalyses;)
	{
		// for each replicate group compute best old likelihood
		oldLogLikelihood = ComputeExperimentLogLikelihood(pData, Concentration, Affinity,
								NewtonCDeriv, NewtonCGrad, IHash, i);
		for (ti=i; ti<pData->m_nReplicates[i]; ti++)
		{
			for (j=0; j<pData->m_nProbes; j++) // generate trial values for each replicate - denser attempts
			{
				if (Affinity[j]>0)
					trialval = (pData->m_fPM[ti][j]-pData->m_fMM[ti][j])/(Affinity[j]);
					//trialval = (pData->m_fPM[nIndex]-pData->m_fMM[nIndex])/(Affinity[j]);
				else
					trialval = -1;
				if (trialval>0)
				{
					// set each replicate within the group containing experiment i to the same trial value
					oldval = Concentration[i];
					for (tti=i; tti<pData->m_nReplicates[i]; tti++)
						Concentration[tti] = trialval;
					LogLikelihood = ComputeExperimentLogLikelihood(pData, Concentration, Affinity,
										NewtonCDeriv, NewtonCGrad, IHash, i);
					if (LogLikelihood<oldLogLikelihood)
					{
						converged = 0; // obviously found something better
						oldLogLikelihood = LogLikelihood; // the new champion
					}
					else // reset each replicate to the same old value
						for (tti=i; tti<pData->m_nReplicates[i]; tti++)
							Concentration[tti] = oldval;
				}
			}
		}
		i = pData->m_nReplicates[i]; // update to next value
	}

	// now do the same for affinities
	for (j=0; j<pData->m_nProbes && !bNoAffinity; j++)
	{
		oldLogLikelihood = ComputeProbeLogLikelihood(pData, Concentration, Affinity,
								NewtonADeriv, NewtonAGrad, IHash, j);
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			// run thru experiments generating trial values
			if (Concentration[i]>0)
				trialval = (pData->m_fPM[i][j]-pData->m_fMM[i][j])/(Concentration[i]);
				//trialval = (pData->m_fPM[nIndex]-pData->m_fMM[nIndex])/(Concentration[i]);
			else
				trialval = -1;
			if (trialval>0)
			{
				oldval = Affinity[j];
				Affinity[j] = trialval;
				LogLikelihood = ComputeProbeLogLikelihood(pData, Concentration, Affinity,
									NewtonADeriv, NewtonAGrad, IHash, j);
				if (LogLikelihood<oldLogLikelihood)
				{
					converged = 0; // obviously found something better
					oldLogLikelihood = LogLikelihood; // the new champion
				}
				else
					Affinity[j] = oldval;
			}
		}
	}

	return(converged);
}

//////////////////////////////////////////////////////////////////////

long CorrectReplicatesSlow(long *NewOrder, long *Replicate, long NumExperiments)
{
	// but, what about replicates?
	long *TOrder = 0;
	long *TRep = 0;
	long i,j, ti, tj;

	TOrder = new long [NumExperiments];
	if (TOrder == 0)
		return NO_DATAMEM;

	TRep = new long [NumExperiments];
	if (TRep == 0)
	{
		delete[] TOrder;
		return NO_DATAMEM;
	}

	// do n^2 correction for replicates
	// which basically is slow, but compared with PLIER iterations is still quick
	// quickest writing of code to get something functional
	// use OldOrder as scratch space here
	j = 0; // start off with nothing
	for (i=0; i<NumExperiments; i++)
	{
		// start with the current least element
		if (NewOrder[i]>-1)
		{
			// add to new order
			tj =j;
			TOrder[j] = NewOrder[i];
			j++;
			// search rest of list
			for (ti=i+1; ti<NumExperiments; ti++)
			{
				if (NewOrder[ti]>-1)
				{
					// if part of the same replicate group, add here
					if (Replicate[NewOrder[ti]]==Replicate[NewOrder[i]])
					{
						TOrder[j] = NewOrder[ti]; // next replicate in order
						j++;
						NewOrder[ti] = -1;
					}
				}
			}
			NewOrder[i] = -1; // done this, turn off
			for (; tj<j; tj++)
				TRep[tj] = j; // set replicate number
		}
	}
	// replace NewOrder
	for (i=0; i<NumExperiments; i++)
	{
		NewOrder[i] = TOrder[i];
		Replicate[i] = TRep[i];
	}

	delete[] TOrder;
	delete[] TRep;
	return NO_PLIER_ERROR;
}

//////////////////////////////////////////////////////////////////////

long SortInputs(plier_data* pData, long *OldOrder)
{
	// sort the input data into a more useful format, keeping the old order in the appropriate position
	// Treats PMdata and MMdata columns as a single, large number
	// must handle ties, due to quantile normalization and possible duplicate values
	// i.e. A vs B is handled by A[0]<>B[0], ties broken by A[1]<>B[1], etc
	double **TmpMatrix = 0; // hold both PM and MM together
	long *NewOrder = 0;
	long i,j;
	long TmpLength;
	long lRes = NO_PLIER_ERROR;

	TmpLength = 2*pData->m_nProbes;
	TmpMatrix = new double * [pData->m_nAnalyses];
	if (TmpMatrix == 0)
		return NO_DATAMEM;

	NewOrder = new long [pData->m_nAnalyses];
	if (NewOrder == 0)
	{
		delete[] TmpMatrix;
		return NO_DATAMEM;
	}

	long nFailedAt = -1;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		TmpMatrix[i] = new double [TmpLength];
		if (TmpMatrix[i] == 0)
		{
			nFailedAt = i;
			break;
		}
	}
	if (nFailedAt != -1)
	{
		for (i=0; i<nFailedAt; i++)
			delete[] TmpMatrix[i];
		delete[] TmpMatrix;
		delete[] NewOrder;
	}

	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nProbes; j++)
		{
			TmpMatrix[i][j] = pData->m_fPM[i][j];
			TmpMatrix[i][j+pData->m_nProbes] = pData->m_fMM[i][j];
		}

	for (i=0; i<pData->m_nAnalyses; i++)
		NewOrder[i] = i;

	// okay, have set up TmpMatrix
	HeapIndexMatrix(TmpMatrix, NewOrder, pData->m_nAnalyses, TmpLength);
	// NewOrder now contains a sort of the experiments

	// but the replicates need to be put in correct order!
	// i.e. adjacent so that the iterators work
	lRes = CorrectReplicatesSlow(NewOrder, pData->m_nReplicates, pData->m_nAnalyses);
	if (lRes == NO_PLIER_ERROR)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
			OldOrder[NewOrder[i]] = i; // inverse map

		for (i=0; i<pData->m_nAnalyses; i++)
			for (j=0; j<pData->m_nProbes; j++)
			{
				pData->m_fPM[i][j] = TmpMatrix[NewOrder[i]][j];
				pData->m_fMM[i][j] = TmpMatrix[NewOrder[i]][j+pData->m_nProbes];
			}
	}
	for (i=0; i<pData->m_nAnalyses; i++)
		delete[] TmpMatrix[i];
	delete[] TmpMatrix;
	delete[] NewOrder;
	return lRes;
}

///////////////////////////////////////////////////////////////////////

long UnScrambleMatrix(plier_data *pData, long *Twist)
{
	long lRes = NO_PLIER_ERROR;
	double *TmpValue;
	long i,j;

	TmpValue = new double [pData->m_nAnalyses];
	if (TmpValue==0)
		return(NO_DATAMEM);

	// unscramble PM/MM values to original order
	for (j=0; j<pData->m_nProbes; j++)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			TmpValue[i] = pData->m_fPM[Twist[i]][j];
		}
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			pData->m_fPM[i][j] = TmpValue[i];
		}
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			TmpValue[i] = pData->m_fMM[Twist[i]][j];
		}
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			pData->m_fMM[i][j] = TmpValue[i];
		}
	}

	delete[] TmpValue;
	return(lRes);
}

/////////////////////////////////////////////////////////////////////////

long UnScrambleReplicates(plier_data *pData, long *Twist)
{
	long lRes = NO_PLIER_ERROR;
	long *TmpValue;
	long i;

	TmpValue = new long [pData->m_nAnalyses];
	if (TmpValue==0)
		return(NO_DATAMEM);

	for (i=0; i<pData->m_nAnalyses; i++)
		TmpValue[i] = pData->m_nReplicates[Twist[i]];
	for (i=0; i<pData->m_nAnalyses; i++)
		pData->m_nReplicates[i] = TmpValue[i];

	delete[] TmpValue;
	return(lRes);
}

//////////////////////////////////////////////////////////////////////

long NewtonPlier(plier_data* pData, double& output)
{
	long i,j;
	double lambda;
	long count, icount;
	double epsilon;
	double dropmax;
	long converged;

	double LogLikelihood;
	double oldLogLikelihood;

	double *NewtonConcentration;
	double *NewtonOldConcentration;
	double *NewtonCDeriv;
	double *NewtonCGrad; // diagonal only
	double *NewtonCStep;
	double *NewtonAffinity;
	double *NewtonOldAffinity;
	double *NewtonADeriv;
	double *NewtonAGrad; // diagonal only
	double *NewtonAStep; // step size

	double **IHash;
	long *OldOrder;

	double lambdalimit;

	output = -1;

	NewtonConcentration = new double [pData->m_nAnalyses]; // experiments
	if (NewtonConcentration == 0)
		return NO_DATAMEM;

	NewtonOldConcentration = new double [pData->m_nAnalyses];
	if (NewtonOldConcentration == 0)
	{
		delete[] NewtonConcentration;
		return NO_DATAMEM;
	}

	NewtonCDeriv = new double [pData->m_nAnalyses]; // direction
	if (NewtonCDeriv == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		return NO_DATAMEM;
	}

	NewtonCGrad = new double [pData->m_nAnalyses]; // how far to go
	if (NewtonCGrad == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		return NO_DATAMEM;
	}

	NewtonCStep = new double [pData->m_nAnalyses];
	if (NewtonCStep == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		return NO_DATAMEM;
	}

	NewtonAffinity = new double [pData->m_nProbes];
	if (NewtonAffinity == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		return NO_DATAMEM;
	}

	NewtonOldAffinity = new double [pData->m_nProbes];
	if (NewtonOldAffinity == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		delete[] NewtonAffinity;
		return NO_DATAMEM;
	}

	NewtonADeriv = new double [pData->m_nProbes];
	if (NewtonADeriv == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		delete[] NewtonAffinity;
		delete[] NewtonOldAffinity;
		return NO_DATAMEM;
	}

	NewtonAGrad = new double [pData->m_nProbes];
	if (NewtonAGrad == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		delete[] NewtonAffinity;
		delete[] NewtonOldAffinity;
		delete[] NewtonADeriv;
		return NO_DATAMEM;
	}

	NewtonAStep = new double [pData->m_nProbes];
	if (NewtonAStep == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		delete[] NewtonAffinity;
		delete[] NewtonOldAffinity;
		delete[] NewtonADeriv;
		delete[] NewtonAGrad;
		return NO_DATAMEM;
	}

	long nSize = pData->m_nAnalyses*pData->m_nProbes;

	IHash = new double * [pData->m_nAnalyses];
	if (IHash == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		delete[] NewtonAffinity;
		delete[] NewtonOldAffinity;
		delete[] NewtonADeriv;
		delete[] NewtonAGrad;
		delete[] NewtonAStep;
		return NO_DATAMEM;
	}

	long nFailedAt = -1;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		IHash[i] = new double [pData->m_nProbes];
		if (IHash[i] == 0)
		{
			nFailedAt = i;
			break;
		}
	}
	if (nFailedAt != -1)
	{
		for (i=0; i<nFailedAt; i++)
			delete[] IHash[i];
		delete[] IHash;
		return NO_DATAMEM;
	}

	OldOrder = new long [pData->m_nAnalyses];
	if (OldOrder == 0)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		delete[] NewtonAffinity;
		delete[] NewtonOldAffinity;
		delete[] NewtonADeriv;
		delete[] NewtonAGrad;
		delete[] NewtonAStep;
		for (i=0; i<pData->m_nAnalyses; i++)
			delete[] IHash[i];
		delete[] IHash;
		return NO_DATAMEM;
	}

	// preprocess for stability
	SortInputs(pData, OldOrder);

	// read in affinities here? Or in previous function
	InitializeVector(NewtonConcentration, pData->m_nAnalyses, pData->m_algParams->defaultconcentration);
	if (pData->m_bUseModel == false)
		InitializeVector(NewtonAffinity, pData->m_nProbes, pData->m_algParams->defaultaffinity);
	else
		TransferVector(NewtonAffinity, pData->m_fAffinity, pData->m_nProbes);

	double augmentation = pData->m_algParams->augmentation;
	AugmentData(pData->m_fPM, pData->m_nAnalyses, pData->m_nProbes, augmentation);
	AugmentData(pData->m_fMM, pData->m_nAnalyses, pData->m_nProbes, augmentation);

	bool bNoAffinity = !pData->m_algParams->fitaffinity;

	// start with SEA for stability!
	long lRes = doSEA(pData, NewtonConcentration, NewtonAffinity, IHash, bNoAffinity);
	if (lRes != NO_PLIER_ERROR)
	{
		delete[] NewtonConcentration;
		delete[] NewtonOldConcentration;
		delete[] NewtonCDeriv;
		delete[] NewtonCGrad;
		delete[] NewtonCStep;
		delete[] NewtonAffinity;
		delete[] NewtonOldAffinity;
		delete[] NewtonADeriv;
		delete[] NewtonAGrad;
		delete[] NewtonAStep;
		for (i=0; i<pData->m_nAnalyses; i++)
			delete[] IHash[i];
		delete[] IHash;
		delete[] OldOrder;
		return lRes;
	}

	if (pData->m_algParams->optimization == FULL_MLE)
	{
		// now update conditions carefully
		for (i=0; i<pData->m_nAnalyses; i++)
			for (j=0; j<pData->m_nProbes; j++)
				IHash[i][j] = 4*pData->m_fPM[i][j]*pData->m_fMM[i][j]; // prepped data for distance

		//Now I can compute PLIER function as follows
		count=0;
		icount = 0;
		LogLikelihood=0;
		oldLogLikelihood = 1;
		epsilon = pData->m_algParams->plierconvergence;
		dropmax = pData->m_algParams->dropmax;
		lambdalimit = pData->m_algParams->lambdalimit;
		converged = 0;
		//beta = 0;

		long iterations = pData->m_algParams->plieriteration;
		while(count<iterations && converged<1)
		{
			LogLikelihood = ComputeGlobalLikelihood(pData, NewtonConcentration, NewtonAffinity,
				NewtonCDeriv, NewtonADeriv, NewtonCGrad, NewtonAGrad, IHash);

			// handle replicates:  Average Deriv and Grad for all replicates, substitute back
			// then everything steps the same amount
			// i.e. treat replicates as identical copies

			Join_Replicates(pData->m_nReplicates, NewtonCDeriv, NewtonCGrad, pData->m_nAnalyses);

			for (i=0; i<pData->m_nAnalyses; i++)
				NewtonCStep[i] = -NewtonCDeriv[i]/(NewtonCGrad[i]); // approximate inverse matrix by inverse diagonal
			for (j=0; j<pData->m_nProbes; j++)
				NewtonAStep[j] = -NewtonADeriv[j]/(NewtonAGrad[j]);

			ShrinkVector(NewtonConcentration, NewtonCStep, pData->m_nAnalyses, dropmax);
			ShrinkVector(NewtonAffinity, NewtonAStep, pData->m_nProbes, dropmax);
			if (bNoAffinity)
				memset(NewtonAStep, 0, pData->m_nProbes*sizeof(double)); // no change to affinities based on data

			oldLogLikelihood = LogLikelihood;
			TransferVector(NewtonOldConcentration, NewtonConcentration, pData->m_nAnalyses);
			TransferVector(NewtonOldAffinity, NewtonAffinity, pData->m_nProbes);

			// try decreasing moves in current direction until succeed
			LogLikelihood +=1; //not yet evaluated, but probably worse
			lambda = 2;
			while (oldLogLikelihood<LogLikelihood && lambda>lambdalimit)
			{
				lambda = lambda/2;
				TransferVector(NewtonConcentration, NewtonOldConcentration, pData->m_nAnalyses);
				TransferVector(NewtonAffinity, NewtonOldAffinity, pData->m_nProbes);
				StepVector(NewtonConcentration, NewtonCStep, pData->m_nAnalyses, lambda);
				StepVector(NewtonAffinity, NewtonAStep, pData->m_nProbes, lambda);

				LogLikelihood = ComputeGlobalLikelihood(pData, NewtonConcentration, NewtonAffinity,
					NewtonCDeriv, NewtonADeriv, NewtonCGrad, NewtonAGrad, IHash);
			}

			if (oldLogLikelihood<LogLikelihood)
			{
				// do no harm
				TransferVector(NewtonConcentration, NewtonOldConcentration, pData->m_nAnalyses);
				TransferVector(NewtonAffinity, NewtonOldAffinity, pData->m_nProbes);
				LogLikelihood = oldLogLikelihood;
			}

			count++;

			if (oldLogLikelihood<LogLikelihood+epsilon)
			{
				converged++;
				TransferVector(NewtonOldConcentration, NewtonConcentration, pData->m_nAnalyses);
				TransferVector(NewtonOldAffinity, NewtonAffinity, pData->m_nProbes);
				oldLogLikelihood = LogLikelihood;
				// found at least a local optimum - try other possible minima
				converged = SearchGridOptimum(pData, NewtonConcentration, NewtonAffinity,
					NewtonCDeriv, NewtonCGrad, NewtonADeriv, NewtonAGrad, IHash,
					epsilon, bNoAffinity);

				LogLikelihood = ComputeGlobalLikelihood(pData, NewtonConcentration, NewtonAffinity,
					NewtonCDeriv, NewtonADeriv, NewtonCGrad, NewtonAGrad, IHash);

				if (oldLogLikelihood<LogLikelihood+epsilon)
				{
					converged = 1; // didn't do any good to find another point
					// test again, because roundoff error can build up if epsilon too small
					if (oldLogLikelihood<LogLikelihood)
					{
						TransferVector(NewtonConcentration, NewtonOldConcentration, pData->m_nAnalyses);
						TransferVector(NewtonAffinity, NewtonOldAffinity, pData->m_nProbes);
					}
				}
				else
					converged =0; // it did work to find another point
			}
			else
				converged = 0;
			BalanceAffinity(NewtonConcentration, NewtonAffinity, pData->m_nAnalyses, pData->m_nProbes);
		}

		if (count == iterations && converged<1)
		{
			lRes = MAXIT_PLIER_REACHED;
		}
	}

	//TransferVector(pData->m_fConcentration, NewtonConcentration, pData->m_nAnalyses);
	ScrambleTransferVector(pData->m_fConcentration, NewtonConcentration, OldOrder, pData->m_nAnalyses);
	TransferVector(pData->m_fAffinity, NewtonAffinity, pData->m_nProbes);

	// clean up PM/MM matrix to be polite to further analysis
	// note that "data augmentation" is still in effect so everything will have augmentation added
	UnScrambleMatrix(pData, OldOrder);
	UnScrambleReplicates(pData, OldOrder);

	delete[] NewtonConcentration;
	delete[] NewtonOldConcentration;
	delete[] NewtonCDeriv;
	delete[] NewtonCGrad;
	delete[] NewtonCStep;
	delete[] NewtonAffinity;
	delete[] NewtonOldAffinity;
	delete[] NewtonADeriv;
	delete[] NewtonAGrad;
	delete[] NewtonAStep;
	for (i=0; i<pData->m_nAnalyses; i++)
		delete[] IHash[i];
	delete[] IHash;
	delete[] OldOrder;

	output = LogLikelihood/nSize;

	return lRes;
}

//////////////////////////////////////////////////////////////////////////////////////////

long Compute_Signed_Residuals(plier_data* pData, long ResidualType)
{
	long i,j;

	// just compute signed error for each PM/MM pair currently
	// do not count concentration penalty (goes with "concentration" related values)
	// do not count affinity penalty (goes with "affinity" related values)
	// do not do Geman-McClure transformation 
	// user can compute weights from raw error values if they like
	// signed residuals more useful than unsigned
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		for (j=0; j<pData->m_nProbes; j++)
		{
			pData->m_fRS[i][j] = JustError(pData->m_fConcentration[i],
				pData->m_fAffinity[j], 
				4*pData->m_fPM[i][j]*pData->m_fMM[i][j], 
				pData->m_fPM[i][j], 
				pData->m_fMM[i][j], 
				pData->m_algParams->usemm);

		}
	}
	return(NO_PLIER_ERROR);
}

