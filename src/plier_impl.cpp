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
 * plier_impl.cpp: define PLIER implementation class
 */

#include <string.h>
#include "affyplier.h"
#include "error.h"

/* error strings */
static const char* szNO_PLIER_ERROR = "No error";
static const char* szNO_DATAMEM = "Failed to allocate memory for performing PLIER analysis. Corrective action is to close other applications and windows and retry.";
static const char* szINVALID_NUM_EXP = "Invalid number of input experiments.";
static const char* szINVALID_NUM_PROBE = "Invalid number of input probes";
static const char* szINVALID_PM = "Invalid perfect match intensities input data pointer";
static const char* szINVALID_MM = "Invalid mismatch intensities input data pointer";
static const char* szINVALID_CONCENTRATION = "Invalid concentration output data pointer";
static const char* szINVALID_AFFINITY = "Invalid affinity output data pointer";
static const char* szINVALID_AUGMENTATION = "Invalid parameter 'augmentation' - It should be greater than zero";
static const char* szINVALID_GMCUTOFF = "Invalid parameter 'gmcutoff' - It should not equal to zero";
static const char* szINVALID_DROPMAX = "Invalid parameter 'dropmax' - It should be greater than zero";
static const char* szINVALID_OPTIMIZATION = "Invalid parameter 'optimization' - It should be either 0 for Full Optimization or 1 for SEA.";
static const char* szINVALID_CONCPENALTY = "Invalid parameter 'concpenalty' - It should be a non-zero value";
static const char* szINVALID_PROBEPENALTY = "Invalid parameter 'probepenalty' - It should be a non-zero value";
static const char* szMAXIT_PLIER_REACHED = "Possible convergence issue in PLIER - maximum iterations hit before convergence detected. Convergence criteria may be too restrictive, maximum iterations too small, or data ill-conditioned.";
static const char* szMAXIT_SEA_REACHED = "Possible convergence issue in SEA - maximum iterations hit before convergence detected. Convergence criteria may be too restrictive, maximum iterations too small, or data ill-conditioned.";

/*
 * PLIER object creation function
 * Input parameter :
 *       plier_impl_name - a plier implementation name
 *       plier           - a pointer to a pointer of iplier object to be returned
 */
void create_plier_object(const char* plier_impl_name, iaffyplier** plier) 
{
	if (plier_impl_name == 0 || strcmp(plier_impl_name, PLIER_IMPL_NAME1) == 0)
	{
		*plier = new caffyplier;
		if (*plier)
			(*plier)->addref();
	}
	else
		*plier = 0;
}

/*
 * get the error string based on the error code
 */
void get_plier_error(long error_code, char* error)
{
	if (!error)
		return;

	switch (error_code)
	{
		case NO_PLIER_ERROR:
			strcpy(error, szNO_PLIER_ERROR);
			break;
		case NO_DATAMEM:
			strcpy(error, szNO_DATAMEM);
			break;
		case INVALID_NUM_EXP:
			strcpy(error, szINVALID_NUM_EXP);
			break;
		case INVALID_NUM_PROBE:
			strcpy(error, szINVALID_NUM_PROBE);
			break;
		case INVALID_PM:
			strcpy(error, szINVALID_PM);
			break;
		case INVALID_MM:
			strcpy(error, szINVALID_MM);
			break;
		case INVALID_CONCENTRATION:
			strcpy(error, szINVALID_CONCENTRATION);
			break;
		case INVALID_AFFINITY:
			strcpy(error, szINVALID_AFFINITY);
			break;
		case INVALID_AUGMENTATION:
			strcpy(error, szINVALID_AUGMENTATION);
			break;
		case INVALID_GMCUTOFF:
			strcpy(error, szINVALID_GMCUTOFF);
			break;
		case INVALID_DROPMAX:
			strcpy(error, szINVALID_DROPMAX);
			break;
		case INVALID_OPTIMIZATION:
			strcpy(error, szINVALID_OPTIMIZATION);
			break;
		case INVALID_CONCPENALTY:
			strcpy(error, szINVALID_CONCPENALTY);
			break;
		case INVALID_PROBEPENALTY:
			strcpy(error, szINVALID_PROBEPENALTY);
			break;
		case MAXIT_PLIER_REACHED:
			strcpy(error, szMAXIT_PLIER_REACHED);
			break;
		case MAXIT_SEA_REACHED:
			strcpy(error, szMAXIT_SEA_REACHED);
			break;
		default:
			strcpy(error, "Unknown error");
			break;
	}
}
