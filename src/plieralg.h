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

#if !defined(_PLIERALGORITHM_HEADER_)
#define _PLIERALGORITHM_HEADER_

#include "affyplier.h"

typedef struct
{
	long					m_nAnalyses;
	long					m_nProbes;
	long*					m_nReplicates;
	double*				m_fConcentration;
	double*				m_fAffinity;
	double**			m_fPM;
	double**			m_fMM;
	double**			m_fRS; // residuals
	bool					m_bUseModel;
	plier_params*	m_algParams;
} plier_data;

////////////////////////////////////////////////////////////////////////

long NewtonPlier(plier_data* pData, double& output);
long Compute_Signed_Residuals(plier_data* pData, long ResidualType);

////////////////////////////////////////////////////////////////////////

#endif // !defined(_PLIERALGORITHM_HEADER_)
