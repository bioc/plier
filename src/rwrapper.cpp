/*
  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
  an improved signal by accounting for experimentally observed patterns 
  in probe behavior and handling error at the appropriately at low and 
  high signal values.

  Copyright (C) 2004 Affymetrix, Inc.

  This program is free software; you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 7
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
 * test.cpp : defines the entry point for the plier interfaces test program
 */

#include <string.h>
#include <vector>
#include "iaffyplier.h"
#include "error.h"
#include <R.h>
#define BUF_LEN	1024
#define MAX_PROBESET_SIZE 2048

/*
 * a plier wrapper function to use the iaffyplier interfaces
 */

void initialise_plier_wrapper(
			      iaffyplier* sp,
			      bool   use_hgu133a_default,
			      double augmentation,
			      double gmcutoff,
			      double probepenalty,
			      double concpenalty,
			      double defaultaffinity,
			      double defaultconcentration,
			      double attenuation,
			      double seaconvergence,
			      long   seaiteration,
			      double plierconvergence,
			      long   plieriteration,
			      bool   usemm,
			      bool   usemodel,
			      bool   fitaffinity,
			      double dropmax,
			      double lambdalimit,
			      long   optimization)
{
  if (sp)
    {
      if (use_hgu133a_default)
	sp->set_default();
      sp->set_init_augmentation(augmentation);
      sp->set_init_defaultaffinity(defaultaffinity);
      sp->set_init_defaultconcentration(defaultconcentration);

      sp->set_sea_attenuation(attenuation);
      sp->set_sea_opt_convergence(seaconvergence);
      sp->set_sea_opt_iteration(seaiteration);

      sp->set_plier_gmcutoff(gmcutoff);
      sp->set_plier_probepenalty(probepenalty);
      sp->set_plier_concpenalty(concpenalty);
      sp->set_plier_use_mm(usemm);
      sp->set_plier_use_model(usemodel);
      sp->set_plier_fit_affinity(fitaffinity);
		
      sp->set_plier_opt_convergence(plierconvergence);
      sp->set_plier_opt_iteration(plieriteration);
      sp->set_plier_opt_dropmax(dropmax);
      sp->set_plier_opt_lambdalimit(lambdalimit);
      sp->set_plier_opt_optimization(optimization);
    }
}

void do_one_probeset_internal(iaffyplier *sp,
			      int      num_exp,
  			      int      num_probe,
			      long*    replicate,
			      double** pm,
			      double** mm,
			      double*  concentration,
			      double*  affinity,
			      int*     error_code)
{
  *error_code = 0;
  sp->set_num_exp(num_exp);
  sp->set_num_probe(num_probe);
  sp->set_replicate(replicate);
  sp->set_pm(pm);
  sp->set_mm(mm);
  sp->set_concentration(concentration);
  sp->set_affinity(affinity);
  long ec;
  sp->run(&ec);
  
  *error_code = (int) ec;
 }


extern "C" {
  void one_probeset(bool    *use_hgu133a_default,
		    double  *augmentation,
		    double  *gmcutoff,
		    double  *probepenalty,
		    double  *concpenalty,
		    double  *defaultaffinity,
		    double  *defaultconcentration,
		    double  *attenuation,
		    double  *seaconvergence,
		    int     *seaiteration,
		    double  *plierconvergence,
		    int     *plieriteration,
		    bool    *usemm,
		    bool    *usemodel,
		    bool    *fitaffinity,
		    double  *dropmax,
		    double  *lambdalimit,
		    int     *optimization,
		    int     *num_exp,
		    int     *num_probes,
		    long    *replicate,
		    double  *pm,
		    double  *mm,
		    double  *concentration,
		    double  *affinity,
		    int     *error_code)
  {

    affy_ptr<iaffyplier> sp;
    create_plier_object(0, &sp);
    /* Create a double** array for pm and mm pointing into the pm and mm vectors for the underlying plier function */
    double **pm_ptr = new double *[*num_exp];
    double **mm_ptr = new double *[*num_exp];
    int i;
    for(i = 0; i < *num_exp; i++) {
      pm_ptr[i] = &pm[*num_probes * i];
      mm_ptr[i] = &mm[*num_probes * i];
    }

    long    *rpct; /* the API expects a long array that it wants to delete itself (replicate). We have an int array that R owns so we need to create the long one here, but let the API handle the clean up. Which is nice. */

    for(int tmpidx=0;tmpidx < *num_exp;tmpidx++) {
      rpct[tmpidx] = (long) replicate[tmpidx];
    }

    initialise_plier_wrapper(sp,
			     *use_hgu133a_default,
			     *augmentation,
			     *gmcutoff,
			     *probepenalty,
			     *concpenalty,
			     *defaultaffinity,
			     *defaultconcentration,
			     *attenuation,
			     *seaconvergence,
			     (long) *seaiteration,
			     *plierconvergence,
			     (long) *plieriteration,
			     *usemm,
			     *usemodel,
			     *fitaffinity,
			     *dropmax,
			     *lambdalimit,
			     (long) *optimization);

    do_one_probeset_internal(sp, *num_exp, *num_probes, rpct, pm_ptr, mm_ptr, concentration, affinity, error_code);
    if (*error_code != 0)
      {
	char error_str[MAX_ERROR_LENGTH];
	get_plier_error(*error_code, error_str); 
	fprintf(stderr,"Error in running plier: %s\n", error_str);
      }
    delete pm_ptr;  
    delete mm_ptr;
  }


  void an_experiment(bool    *use_hgu133a_default,
		     double  *augmentation,
		     double  *gmcutoff,
		     double  *probepenalty,
		     double  *concpenalty,
		     double  *defaultaffinity,
		     double  *defaultconcentration,
		     double  *attenuation,
		     double  *seaconvergence,
		     int     *seaiteration,
		     double  *plierconvergence,
		     int     *plieriteration,
		     bool    *usemm,
		     bool    *usemodel,
		     bool    *fitaffinity,
		     double  *dropmax,
		     double  *lambdalimit,
		     int     *optimization,
		     int     *num_exp,              /* how many experiments? */
		     int     *num_probes,           /* length of pm, mm is num_exp * this value */
		     int     *replicate,            /* replicate structure of the experiment */
		     double  *pm,                   /* the pm and mm data for the whole experiment */
		     double  *mm,                
		     char    **probenames,          /* there should be num-probes of these */
		     double  *concentration,        /* the result end up in here */
		     double  *affinity,             /* and here */
		     int     *error_code)           /* of course, this should only ever be '0' :-)) */
  {
    affy_ptr<iaffyplier> sp;
    create_plier_object(0, &sp);
    initialise_plier_wrapper(sp,
			     *use_hgu133a_default,
			     *augmentation,
			     *gmcutoff,
			     *probepenalty,
			     *concpenalty,
			     *defaultaffinity,
			     *defaultconcentration,
			     *attenuation,
			     *seaconvergence,
			     (long) *seaiteration,
			     *plierconvergence,
			     (long) *plieriteration,
			     *usemm,
			     *usemodel,
			     *fitaffinity,
			     *dropmax,
			     *lambdalimit,
			     (long) *optimization);
    /* The API expects a two-dimenional array of probes with the x axis corresponding to each probe, the y axis corresponding to each sample */
    /* R gives us a one dimensional array with the probes contigous and each sample following on from the last */
    /* most of this function is spent re-arranging the R data structure to give us the 2D array for each probeset that the API needs */


    /* a cache into which to gather the pm and mm probe data for each probeset */
    double *pms = 0;
    double *mms = 0;
    double **pm_ptr = 0;
    double **mm_ptr = 0;
    long    *rpct = new long[*num_exp]; /* the API expects a long array that it wants to delete itself (replicate). We have an int array that R owns so we need to create the long one here, but let the API handle the clean up. Which is nice. */

    for(int tmpidx=0;tmpidx < *num_exp;tmpidx++) {
      rpct[tmpidx] = (long) replicate[tmpidx];
    }

    pms = (double *) R_alloc(MAX_PROBESET_SIZE * *num_exp,sizeof(double));
    if(!pms) {
      fprintf(stderr,"Not enough memory to allocate pm cache."); 
      exit(1); 
    }
    mms = (double *) R_alloc(MAX_PROBESET_SIZE * *num_exp,sizeof(double));
    if(!mms) {
      fprintf(stderr,"Not enough memory to allocate mm cache."); 
      exit(1); 
    }

    /* Create a double** array for pm and mm pointing into the pm and mm vectors for the underlying plier function */
    pm_ptr = (double **) R_alloc(*num_exp,sizeof(double *));
    if(!pm_ptr) {
      fprintf(stderr,"Not enough memory to allocate pm_ptr cache."); 
      exit(1); 
    }
    mm_ptr = (double **) R_alloc(*num_exp,sizeof(double *));
    if(!mm_ptr) {
      fprintf(stderr,"Not enough memory to allocate mm_ptr cache."); 
      exit(1); 
    }

    /* run through the probenames array to identify each probeset one
       by one. As we do that, fill up pms and mms and count the number
       of probes in the probeset so far. Also stitch together the
       pm_ptr and mm_ptr arrays (which exist because the plier SDK we
       call expects a 2D array) */

    /* Do this initially for the first probeset */

    int start = 0;
    int i; /* used to count through all the probes in the experiment */
    int j; /* used to count through the experiments as we shuffle probes around */
    int k; /* used to count the number of probes we've seen in the current probeset as we shuffle probes around */
    int current = 0;

    /* do the first probe */
    for(j = 0; j < *num_exp; j++) {
      pms[j * MAX_PROBESET_SIZE] = pm[j * *num_probes]; 
      mms[j * MAX_PROBESET_SIZE] = mm[j * *num_probes];

      pm_ptr[j] = &pms[j * MAX_PROBESET_SIZE];
      mm_ptr[j] = &mms[j * MAX_PROBESET_SIZE];
    }
    k = 1;

    /* run through the array of probes looking at the names associated with them */
    /* each time we see a new name we've reached a new probeset */

    /* So... for each probeset build the 2D arrays of pm and mm  probes ready to pass through to the API */

    for(i = 1; i < *num_probes; i++) {

      if(strcmp(probenames[i],probenames[start]) == 0) {
        /* same probeset - copy over the next set of probes */  
	for(j = 0; j < *num_exp; j++) {
	  pms[k + j * MAX_PROBESET_SIZE] = pm[i + j * *num_probes]; 
	  mms[k + j * MAX_PROBESET_SIZE] = mm[i + j * *num_probes];
	}
        k++; /* get ready for the next probes, unless... */
        if(k > MAX_PROBESET_SIZE) fprintf(stderr,"Error in running plier: MAX_PROBESET_SIZE exceeded %d\n",k);
      }
      else {	
        /* onto a new probeset */
        
        /* so get the expression calls for the last one */
	do_one_probeset_internal(sp, *num_exp, k, rpct, pm_ptr, mm_ptr, &concentration[current * *num_exp], &affinity[start], error_code);
 
        /* and copy over the zero'th probeset into the temporary data structures */ 
	for(j = 0; j < *num_exp; j++) {
	  pms[j * MAX_PROBESET_SIZE] = pm[i + j * *num_probes]; 
	  mms[j * MAX_PROBESET_SIZE] = mm[i + j * *num_probes];
	}

	k=1;
        current ++; /* next result */
	start = i;  /* remember the start point of this probeset */
        /* comfort dots ... */
	if(current % 1000 == 0) fprintf(stderr,".",probenames[i]);
        if(k > MAX_PROBESET_SIZE) fprintf(stderr,"Error in running plier: MAX_PROBESET_SIZE exceeded %d\n",k);
      }
    }

    /* Mop up the last probeset left as we fall out of the loop */
    do_one_probeset_internal(sp, *num_exp, k, rpct, pm_ptr, mm_ptr, &concentration[current * *num_exp], &affinity[start], error_code);

    if (*error_code != 0)
      {
	char error_str[MAX_ERROR_LENGTH];
	get_plier_error(*error_code, error_str); 
	fprintf(stderr,"Error in running plier: %s\n", error_str);
      }

    fprintf(stderr," done.\n");
  }
}


