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

/*
 * plier.cpp: implementation of the caffyplier class.
 */

#include <memory.h>
#include <string.h>
#include "affyplier.h"
#include "plieralg.h"
#include "error.h"

/* default parameters for HG-U133A */
#define hgu133a_augmentation	0.1
#define hgu133a_defaultaffinity	1.0
#define hgu133a_defaultconcentration	1.0
#define hgu133a_attenuation	0.005f
#define hgu133a_seaconvergence	0.000001
#define hgu133a_seaiteration	2000
#define hgu133a_gmcutoff	0.15f
#define hgu133a_probepenalty	0.001f
#define hgu133a_concpenalty	0.000001f
#define hgu133a_usemm 1
#define hgu133a_usemodel 0
#define hgu133a_fitaffinity 1
#define hgu133a_plierconvergence	0.000001
#define hgu133a_plieriteration	3000
#define hgu133a_dropmax	3.0
#define hgu133a_lambdalimit	0.01
#define hgu133a_optimization 0


/***************************************************************/
/** caffyplier class implementation                           **/
/***************************************************************/

/*
 * Constructor
 */
caffyplier::caffyplier()
{
	memset(&params, 0, sizeof(plier_params));
	num_exp = 0;
	num_probe = 0;
	pm = NULL;
	mm = NULL;
	residuals = NULL;
	replicate = NULL;
	concentration = NULL;
	affinity = NULL;
	set_default();
}

/*
 * Destructor
 */
caffyplier::~caffyplier()
{
}

/*
 * set to HG-U133A default parameters
 */
void caffyplier::set_default(void* buffer)
{
	params.augmentation = hgu133a_augmentation;
	params.defaultaffinity = hgu133a_defaultaffinity;
	params.defaultconcentration = hgu133a_defaultconcentration;
	params.attenuation = hgu133a_attenuation;
	params.seaconvergence = hgu133a_seaconvergence;
	params.seaiteration = hgu133a_seaiteration;
	params.gmcutoff = hgu133a_gmcutoff;
	params.probepenalty = hgu133a_probepenalty;
	params.concpenalty = hgu133a_concpenalty;
	params.usemm = hgu133a_usemm;
	params.usemodel = hgu133a_usemodel;
	params.fitaffinity = hgu133a_fitaffinity;
	params.plierconvergence = hgu133a_plierconvergence;
	params.plieriteration = hgu133a_plieriteration;
	params.dropmax = hgu133a_dropmax;
	params.lambdalimit = hgu133a_lambdalimit;
	params.optimization = hgu133a_optimization;
}

// Like above but without the useless pointer
void caffyplier::set_default_u133() {
  set_default();
}

/*
 * validate the input parameters before calling PLIER algorithm
 */
long caffyplier::validate_params()
{
	if (params.augmentation < 0)
		return INVALID_AUGMENTATION;

	/* gmcutoff should not be 0, otherwise, division by zero will occur */
	if (params.gmcutoff == 0)
		return INVALID_GMCUTOFF;

	/* dropmax should be > 0, otherwise, division by zero will occur */
	if (params.dropmax <= 0)
		return INVALID_DROPMAX;

	/* concpenalty should be non-zero */
	if (params.concpenalty == 0)
		return INVALID_CONCPENALTY;

	/* probepenalty should be non-zero */
	if (params.probepenalty == 0)
		return INVALID_PROBEPENALTY;

	if (params.optimization != FULL_MLE && params.optimization != SEA)
		return INVALID_OPTIMIZATION;

	/* seaiteration should be > 0 */
	if (params.seaiteration <= 0)
		return INVALID_SEAITERATION;

	/* plieriteration should be > 0	*/
	if (params.optimization == FULL_MLE && params.plieriteration <= 0)
		return INVALID_PLIERITERATION;
	
	return NO_PLIER_ERROR;
}

/*
 * validate the inputs before calling PLIER algorithm
 */
long caffyplier::validate_inputs()
{
	if (num_exp <= 0)
		return INVALID_NUM_EXP;

	if (num_probe <= 0)
		return INVALID_NUM_PROBE;

	if (pm == 0)
		return INVALID_PM;

	if (mm == 0)
		return INVALID_MM;

	// if residuals null, assume no return required.

	if (concentration == 0)
		return INVALID_CONCENTRATION;

	if (affinity == 0)
		return INVALID_AFFINITY;

	return NO_PLIER_ERROR;
}

/*
 * run the plier algorithm
 */

void caffyplier::run(long* error_code)
{
	*error_code = validate_params();
	if (*error_code != NO_PLIER_ERROR)
		return;

	*error_code = validate_inputs();
	if (*error_code != NO_PLIER_ERROR)
		return;

	bool create_def_replicate = false;
	if (!replicate)
	{
		replicate = new long[num_exp];
		if (replicate == 0)
		{
			*error_code = NO_DATAMEM;
			return;
		}
		set_default_replicate(num_exp, replicate);
		create_def_replicate = true;
	}

	plier_data inputs;
	inputs.m_nAnalyses = num_exp;
	inputs.m_nProbes = num_probe;
	inputs.m_nReplicates = replicate;
	inputs.m_fPM = pm;
	inputs.m_fMM = mm;
	inputs.m_fRS = residuals;
	inputs.m_fAffinity = affinity;
	inputs.m_fConcentration = concentration;
	inputs.m_bUseModel = params.usemodel;
	inputs.m_algParams = &params;

	double output;
	*error_code = NewtonPlier(&inputs, output);

	// now we have affinities/concentrations fit, if we have residuals, fit them
	if (*error_code == NO_PLIER_ERROR && residuals!=0)
	{
		// construct signed residuals & return 
		// assumptions - ignore concentration/affinity penalties
		// want just the fit of the data
		// concentrations, affinities, >augmented< data, as returned from NewtonPlier
		// get raw signed residuals [no Geman-McClure] for each PM/MM value & return
		*error_code = Compute_Signed_Residuals(&inputs, 0);
	}

	if (create_def_replicate && replicate)
	{
		delete[] replicate;
		replicate = 0;
	}
}

/* 
 * set to default replicate groups
 */
void caffyplier::set_default_replicate(long size, long* replicate)
{
	if (replicate)
		for (long i=0; i<size; i++)
			replicate[i] = i;
}

/*
 * get the implementation name and version
 */
void caffyplier::get_impl_info(char* plier_impl_name, char* plier_impl_ver)
{
	try
	{
		strcpy(plier_impl_name, PLIER_IMPL_NAME1);
		strcpy(plier_impl_ver, PLIER_IMPL_VER1);
	}
	catch(...)
	{
	}
}
