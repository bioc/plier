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

/* ADDED By Crispin */

#define _GCC

/* --- */



 /*
 * iaffyplier.h: iaffyplier interface class.
 *
 */

#if !defined(PLIER_INTERFACE_HEADER)
#define PLIER_INTERFACE_HEADER

#ifdef _MSVC
#define plierdll __declspec( dllexport )
#else
#ifdef _GCC
#define plierdll
#else
#define plierdll __declspec( dllimport )
#endif /* _GCC */
#endif /* _MSVC */

#include "affy_ptr.h"

#define FULL_MLE 0
#define SEA 1
#define MAX_IMPLNAME_LENGTH 255
#define MAX_IMPLVER_LENGTH 32
#define MAX_ERROR_LENGTH 1024

/*
 * PLIER Interfaces
 * An abstract C++ class
 */

class iaffyplier
{
protected:
	iaffyplier() { refcnt = 0; }
	virtual ~iaffyplier() {}

public:

	/*
	 * Increment reference count
	 */
	void addref() { ++refcnt; }

	/*
	 * Decrement reference count, if zero, delete the current object
	 */
	void release() 
	{ 
		if (refcnt > 0 && --refcnt == 0)
			delete this;
	}

	/* PLIER parameters - Initialization */
	virtual void set_init_augmentation(double) = 0;          /* avoid zero values in input data */
	virtual void get_init_augmentation(double*) = 0;         /* return byref the value of augmentation */
	virtual void set_init_defaultaffinity(double) = 0;       /* default affinities if not supplied */
	virtual void get_init_defaultaffinity(double*) = 0;      /* return byref the value of default affinity */
	virtual void set_init_defaultconcentration(double) = 0;  /* default concentration if not supplied */
	virtual void get_init_defaultconcentration(double*) = 0; /* return byref the value of default concentation */

	/* PLIER parameters - SEA */
	virtual void set_sea_attenuation(float) = 0;             /* set up attenuation of background */
	virtual void get_sea_attenuation(float*) = 0;            /* return byref the value of SEA attenuation */

	/* PLIER parameters - SEA Optimization */
	virtual void set_sea_opt_convergence(double) = 0;        /* change in log-value at which to stop */
	virtual void get_sea_opt_convergence(double*) = 0;       /* return byref the value of SEA convergence */
	virtual void set_sea_opt_iteration(long) = 0;            /* max number of SEA iteration to avoid infinite loops in SEA */
	virtual void get_sea_opt_iteration(long*) = 0;           /* return byref the value of SEA iteration */

	/* PLIER parameters - PLIER */
	virtual void set_plier_gmcutoff(float) = 0;              /* set up robustness */
	virtual void get_plier_gmcutoff(float*) = 0;             /* return byref the value of GMCutoff */
	virtual void set_plier_probepenalty(float) = 0;          /* Bayes penalty for peculiar probes */
	virtual void get_plier_probepenalty(float*) = 0;         /* return byref the value of probe penalty */
	virtual void set_plier_concpenalty(float) = 0;           /* Bayes penalty for really peculiar concentrations */
	virtual void get_plier_concpenalty(float*) = 0;          /* return byref the value of concentration penalty */
	virtual void set_plier_use_mm(bool) = 0;                 /* use mm or background based likelihood */
	virtual void get_plier_use_mm(bool*) = 0;                /* return byref the value of usemm */
	virtual void set_plier_use_model(bool) = 0;              /* do not update affinities but use fixed values */
	virtual void get_plier_use_model(bool*) = 0;             /* return byref the value of usemodel */
	virtual void set_plier_fit_affinity(bool) = 0;           /* update affinities */
	virtual void get_plier_fit_affinity(bool*) = 0;          /* return byref the value of fitaffinity */

	/* PLIER parameters - PLIER Optimization */
	virtual void set_plier_opt_convergence(double) = 0;      /* change in likelihood */
	virtual void get_plier_opt_convergence(double*) = 0;     /* return byref the value of PLIER convergence */
	virtual void set_plier_opt_iteration(long) = 0;          /* max number of PLIER iteration to avoid infinite loops in PLIER */
	virtual void get_plier_opt_iteration(long*) = 0;         /* return byref the value of PLIER iteration */
	virtual void set_plier_opt_dropmax(double) = 0;          /* "logarithmic" barrier near zero */
	virtual void get_plier_opt_dropmax(double*) = 0;         /* return byref the value of dropmax */
	virtual void set_plier_opt_lambdalimit(double) = 0;      /* minimum step multiplier in method */
	virtual void get_plier_opt_lambdalimit(double*) = 0;     /* return byref the value of lambdalimit */
	virtual void set_plier_opt_optimization(long) = 0;       /* optimization method. either FULL_OPTIMIZATION or SEA */
	virtual void get_plier_opt_optimization(long*) = 0;      /* return byref the value of optimization */

	/* PLIER Inputs */
	virtual void set_num_exp(long) = 0;                      /* number of input experiments */
	virtual void get_num_exp(long*) = 0;                     /* return byref the number of input experiments */
	virtual void set_num_probe(long) = 0;                    /* number of input probes in each experiment */
	virtual void get_num_probe(long*) = 0;                   /* return byref the number of input probes in each experiment */
	virtual void set_replicate(long*) = 0;                   /* replicate groups */
	virtual void get_replicate(long**) = 0;                  /* return byref the pointer that set by using "set_replicate" method */
	virtual void set_pm(double**) = 0;                       /* matrix of num_exp by num_probe stored by experiments (i.e. by row) */
	virtual void get_pm(double***) = 0;                      /* return byref the pointer that set by using "set_pm" method */
	virtual void set_mm(double**) = 0;                       /* matrix of mismatch or background values, corresponds to structure of pm */
	virtual void get_mm(double***) = 0;                      /* return byref the pointer that set by using "set_mm" method */

	/* PLIER Inputs/Outputs */
	virtual void set_concentration(double*) = 0;             /* array of length num_exp for estimated intensities */
	virtual void get_concentration(double**) = 0;            /* return byref the pointer that set by using "set_est_intnsity" method */
	virtual void set_affinity(double*) = 0;                  /* array of length num_probe for estimated probe effects */
	virtual void get_affinity(double**) = 0;                 /* return byref the pointer that set by using "set_est_probe_effect" method */

	/* PLIER Parameters (RESERVED FOR FUTURE) */
	virtual void set_use_concentration(bool) = 0;            /* 1 if initial estimates of concentration are provided */
	virtual void get_use_concentration(bool*) = 0;           /* return byref the value of the use_concentration flag */
	virtual void set_use_bg(bool) = 0;                       /* 0 if the bg data is not initialized */
	virtual void get_use_bg(bool*) = 0;                      /* return byref the value of usebg */
	virtual void set_use_conc_sd(bool) = 0;                  /* 1 if initial estimates of concentrations' standard deviation are provided */
	virtual void get_use_conc_sd(bool*) = 0;                 /* return byref the value of the use_conc_sd flag */
	virtual void set_use_affinity_sd(bool) = 0;              /* 1 if estimates of affinities' standard deviation are provided */
	virtual void get_use_affinity_sd(bool*) = 0;             /* return byref the value of the use_affinity_sd flag */

	/* PLIER Inputs (RESERVED FOR FUTURE) */
	virtual void set_bg(double**) = 0;                       /* currently NULL, intended to allow for provision of background, corresponds to structure of pm */
	virtual void get_bg(double***) = 0;                      /* return byref the pointer that set by using "set_bg" method */

	/* PLIER Inputs/Outputs (RESERVED FOR FUTURE) */
	virtual void set_conc_sd(double*) = 0;                   /* array of length num_exp for standard deviation of estimated intensities */
	virtual void get_conc_sd(double**) = 0;                  /* return byref the pointer that set by using "set_signal_sd" method */
	virtual void set_affinity_sd(double*) = 0;               /* array of length num_probe for standard deviation of estimated probe effects */
	virtual void get_affinity_sd(double**) = 0;              /* return byref the pointer that set by using "set_affinity_sd" method */

	/* PLIER Outputs (RESERVED FOR FUTURE) */
	virtual void set_residuals(double**) = 0;                /* Used to return residuals, structure corresponds to pm */
	virtual void get_residuals(double***) = 0;               /* return byref the pointer that set by using "set_residuals" method */
	virtual void set_call_pvalue(double*) = 0;               /* array of length num_exp to return detection pvalues */
	virtual void get_call_pvalue(double**) = 0;              /* return byref the pointer that set by using "set_call_pvalue" method */
	virtual void set_likelihood(double*) = 0;                /* The likelihood, or some similarly-interpretable quantity */
	virtual void get_likelihood(double**) = 0;               /* return byref the pointer that set by using "set_likelihood" method */

	/*
	 * Set to HG-U133A default PLIER parameters
	 * Input parameters :
	 *       buffer - [in] Reserved; must be 0
	 */
	virtual void set_default(void* buffer=0) = 0;

	/*
	 * Interface to run the PLIER algorithm to obtain outputs
	 * Input parameter :
	 *       error_code - [out] Pointer to a long integer which receives
	 *                    the error code value, set to 0 if no error
	 */
	virtual void run(long* error_code) = 0;

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
	virtual void get_impl_info(char* plier_impl_name, char* plier_impl_ver) = 0;


private:
	int refcnt;  // reference count
};

/*
 * PLIER object creation function
 * Input parameters :
 *       plier_impl_name - [in] Pointer to a buffer that stores the plier implementation name.
 *                         It can be null and then default to PLIER v1 implementation
 *                         The buffer size should be at most MAX_IMPLNAME_LENGTH + 1 characters
 *
 *       plier           - [out] Pointer to a pointer of iaffyplier object to be returned
 */
plierdll void create_plier_object(const char* plier_impl_name, iaffyplier** plier);

/*
 * Get the error string
 * Input parameters :
 *       error_code - [in] Specifies the error code returned in run method
 *
 *       error      - [out] Pointer to a buffer that receives
 *                    a null-terminated string containing the error string.
 *                    The buffer size should be large enough to contain
 *                    MAX_ERROR_LENGTH + 1 characters
 */
plierdll void get_plier_error(long error_code, char* error);

#endif /* defined(PLIER_INTERFACE_HEADER) */
