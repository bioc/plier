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
 * affyplier.h: interface for the caffyplier class.
 */

#if !defined(PLIERCLASS_HEADER)
#define PLIERCLASS_HEADER

#include "iaffyplier.h"

#define PLIER_IMPL_NAME1 "PLIER1"
#define PLIER_IMPL_VER1 "1.0"

/*
 * structures for storing PLIER parameters
 */

typedef struct
{
	double augmentation;
	double defaultaffinity;
	double defaultconcentration;
	double seaconvergence;
	double plierconvergence;
	double dropmax;
	double lambdalimit;
	float attenuation;
	float gmcutoff;
	float probepenalty;
	float concpenalty;
	bool usemm;
	bool usemodel;
	bool fitaffinity;
	long seaiteration;
	long plieriteration;
	long optimization;
} plier_params;

/*
 * PLIER implementation class which inherits from iplier interface class
 */

class caffyplier : public iaffyplier
{
	/* member variables */
	long num_exp;
	long num_probe;
	double** pm;
	double** mm;
	double** residuals;
	long* replicate;
	double* concentration;
	double* affinity;
	plier_params params;

public:
	caffyplier();
	virtual ~caffyplier();

public:

	/* PLIER parameters - Initialization */
	virtual void set_init_augmentation(double val) { params.augmentation=val; }
	virtual void get_init_augmentation(double* outval) { *outval=params.augmentation; }
	virtual void set_init_defaultaffinity(double val) { params.defaultaffinity=val; }
	virtual void get_init_defaultaffinity(double* outval) { *outval=params.defaultaffinity; }
	virtual void set_init_defaultconcentration(double val) {params.defaultconcentration=val; }
	virtual void get_init_defaultconcentration(double* outval) { *outval=params.defaultconcentration; }

	/* PLIER parameters - SEA */
	virtual void set_sea_attenuation(float val) { params.attenuation=val; }
	virtual void get_sea_attenuation(float* outval) { *outval=params.attenuation; }

	/* PLIER parameters - SEA Optimization */
	virtual void set_sea_opt_convergence(double val) { params.seaconvergence=val; }
	virtual void get_sea_opt_convergence(double* outval) { *outval=params.seaconvergence; }
	virtual void set_sea_opt_iteration(long val) { params.seaiteration=val; }
	virtual void get_sea_opt_iteration(long* outval) { *outval=params.seaiteration; }

	/* PLIER parameters - PLIER */
	virtual void set_plier_gmcutoff(float val) { params.gmcutoff=val; }
	virtual void get_plier_gmcutoff(float* outval) { *outval=params.gmcutoff; }
	virtual void set_plier_probepenalty(float val) { params.probepenalty=val; }
	virtual void get_plier_probepenalty(float* outval) { *outval=params.probepenalty; }
	virtual void set_plier_concpenalty(float val) { params.concpenalty=val; }
	virtual void get_plier_concpenalty(float* outval) { *outval=params.concpenalty; }
	virtual void set_plier_use_mm(bool val) { params.usemm=val; }
	virtual void get_plier_use_mm(bool* outval) { *outval=params.usemm; }
	virtual void set_plier_use_model(bool val) { params.usemodel=val; }
	virtual void get_plier_use_model(bool* outval) { *outval=params.usemodel; }
	virtual void set_plier_fit_affinity(bool val) { params.fitaffinity=val; }
	virtual void get_plier_fit_affinity(bool* outval) { *outval=params.fitaffinity; }
	
	/* PLIER parameters - PLIER Optimization */
	virtual void set_plier_opt_convergence(double val) { params.plierconvergence=val; }
	virtual void get_plier_opt_convergence(double* outval) { *outval=params.plierconvergence; }
	virtual void set_plier_opt_iteration(long val) { params.plieriteration=val; }
	virtual void get_plier_opt_iteration(long* outval) { *outval=params.plieriteration; }
	virtual void set_plier_opt_dropmax(double val) { params.dropmax=val; }
	virtual void get_plier_opt_dropmax(double* outval) { *outval=params.dropmax; }
	virtual void set_plier_opt_lambdalimit(double val) { params.lambdalimit=val; }
	virtual void get_plier_opt_lambdalimit(double* outval) { *outval=params.lambdalimit; }
	virtual void set_plier_opt_optimization(long val) { params.optimization=val; }
	virtual void get_plier_opt_optimization(long *outval) { *outval=params.optimization; }

	/* PLIER Inputs */
	virtual void set_num_exp(long val) { num_exp=val; }
	virtual void get_num_exp(long* val) { *val=num_exp; }
	virtual void set_num_probe(long val) { num_probe=val; }
	virtual void get_num_probe(long* val) { *val=num_probe; }
	virtual void set_replicate(long* val) { replicate=val; }
	virtual void get_replicate(long** val) { *val=replicate; }
	virtual void set_pm(double** val) { pm=val; }
	virtual void get_pm(double*** val) { *val=pm; }
	virtual void set_mm(double** val) { mm=val; }
	virtual void get_mm(double*** val) { *val=mm; }

	/* PLIER Inputs/Outputs */
	virtual void set_concentration(double* val) { concentration=val; }
	virtual void get_concentration(double** val) { *val=concentration; }
	virtual void set_affinity(double* val) { affinity=val; }
	virtual void get_affinity(double** val) { *val=affinity; }

	/* PLIER Parameters (RESERVED FOR FUTURE) */
	virtual void set_use_concentration(bool) {};
	virtual void get_use_concentration(bool*outval) {};
	virtual void set_use_bg(bool) {};
	virtual void get_use_bg(bool*outval) {};
	virtual void set_use_conc_sd(bool) {};
	virtual void get_use_conc_sd(bool*outval) {};
	virtual void set_use_affinity_sd(bool) {};
	virtual void get_use_affinity_sd(bool*outval) {};

	/* PLIER Inputs (RESERVED FOR FUTURE) */
	virtual void set_bg(double**) {};
	virtual void get_bg(double***) {};

	/* PLIER Inputs/Outputs (RESERVED FOR FUTURE) */
	virtual void set_conc_sd(double*) {};
	virtual void get_conc_sd(double**) {};
	virtual void set_affinity_sd(double*) {};
	virtual void get_affinity_sd(double**) {};

	// IMPLEMENTED FOR UTILITY
	virtual void set_residuals(double** val) {residuals=val;};
	virtual void get_residuals(double*** val) { *val=residuals;};
	/* PLIER Outputs (RESERVED FOR FUTURE) */
	virtual void set_call_pvalue(double*) {};
	virtual void get_call_pvalue(double**) {};
	virtual void set_likelihood(double*) {};
	virtual void get_likelihood(double**) {};

	/*
	 * Set to HG-U133A default PLIER parameters
	 * Input parameters :
	 *       buffer - [in] Reserved; must be 0
	 */
	virtual void set_default(void* buffer=0);
  virtual void set_default_u133();

	/*
	 * Interface to run the PLIER algorithm to obtain outputs
	 * Input parameter :
	 *       error_code - [out] Pointer to a long integer which receives
	 *                    the error code value, set to 0 if no error
	 */
	virtual void run(long* error_code);

	/*
	 * Get the plier object implementation name and version number
	 * Input parameters :
	 *       plier_impl_name - [out] Pointer to a buffer that receives
	 *                         a null-terminated string containing the implementation name.
	 *                         The buffer size should be large enough to contain
	 *                         MAX_IMPLNAME_LENGTH + 1 characters
	 *
	 *       plier_impl_ver  - [out] Pointer to a buffer that receives
	 *                         a null-terminated string containing the 
	 *                         The buffer size should be large enough to contain
	 *                         MAX_IMPLVER_LENGTH + 1 characters
	 */
	virtual void get_impl_info(char* plier_impl_name, char* plier_impl_ver);

protected:

	/* 
	 * if input replicate is NULL, set it to default replicate
	 * i.e. each experiment is in its own replicate group
	 */
	virtual void set_default_replicate(long size, long* replicate);

	/*
	 * validate the input parameters before running the algorithm
	 * return 0 if validation ok
	 */
	virtual long validate_params();

	/*
	 * validate the inputs before running the algorithm
	 * return 0 if validation ok
	 */
	virtual long validate_inputs();

};

#endif /* !defined(PLIERCLASS_HEADER) */
