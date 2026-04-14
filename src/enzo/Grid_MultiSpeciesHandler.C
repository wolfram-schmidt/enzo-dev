/***********************************************************************
/
/  GRID CLASS (HANDLE CALLING AND SOLVING COOLING/CHEMISTRY)
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: Move logic for chemistry/cooling module selection here
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::MultiSpeciesHandler()
{
  if ((!MultiSpecies) && (!RadiativeCooling) && (!UseSGSDiffusion)) return SUCCESS; 
  if (GadgetEquilibriumCooling != 0) return SUCCESS;

  LCAPERF_START("grid_MultiSpeciesHandler");

#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == TRUE) {
    grackle_data->radiative_transfer_intermediate_step = FALSE;
    if (this->GrackleWrapper() == FAIL) {
      ENZO_FAIL("Error in GrackleWrapper.\n");
    }
  }
#else
  if (MultiSpecies && RadiativeCooling ) {
    int RTCoupledSolverIntermediateStep = FALSE;
    this->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep);
  } else {
    if (MultiSpecies)
      this->SolveRateEquations();
    if (RadiativeCooling)
      this->SolveRadiativeCooling();
  }
#endif

  if (HydroMethod != HD_RK && HydroMethod != MHD_RK && UseSGSModel && UseSGSDiffusion) {

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);

    // internal energy gradient
    if (DualEnergyFormalism) {
	    if (this->SGSUtil_ComputeGradient(GradEint,BaryonField[GENum]) == FAIL) {
	      fprintf(stderr, "grid::MultiSpeciesHandler: Error in SGSUtil_ComputeGradient(Eint)).\n");
	      return FAIL;
	    }
	  } else {
      // compute internal energy and store in AuxField
      if (this->SGSUtil_InternalEnergy() == FAIL) {
        fprintf(stderr, "grid::MultiSpeciesHandler: Error in SGSUtil_InternalEnergy.\n");
	      return FAIL;
      }
	    if (this->SGSUtil_ComputeGradient(GradEint,AuxField) == FAIL) {
	      fprintf(stderr, "grid::MultiSpeciesHandler: Error in SGSUtil_ComputeGradient(Eint)).\n");
	      return FAIL;
	    }
	  }

    int DeNum = FindField(ElectronDensity, FieldType, NumberOfBaryonFields);    
    if (debug1 && DeNum >= 0)
      printf("Free electron field: %"ISYM"\n",DeNum);

    int HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields);
    if (debug1 && HINum >= 0)
      printf("HI field: %"ISYM"\n",HINum);

    int MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
    if (debug1 && MetalNum >= 0)
      printf("Metal field: %"ISYM"\n",MetalNum);

    int ns_max = NEQ_HYDRO + NSpecies;
    if (MetalNum >= 0)
      ns_max++;

    if (debug1)
      printf("Number of baryon fields: %"ISYM", %"ISYM" hydro, %"ISYM" species, %"ISYM" all\n",
             NumberOfBaryonFields,NEQ_HYDRO,NSpecies,ns_max);

    // species gradients (excluding free electron field and colors)
	  for (int ns = NEQ_HYDRO, s = 0; ns < ns_max; ns++, s++) {

        if (ns == DeNum)
          continue;

        if (debug1)
          printf("Computing gradient of species %"ISYM", %"ISYM"\n",ns,s);

        // change species from density to mass fraction and store in AuxField
        if (this->SGSUtil_MassFraction(ns) == FAIL) {
          fprintf(stderr, "grid::MassFraction: Error in SGSUtil_MassFraction.\n");
	        return FAIL;
        } 
	      if (this->SGSUtil_ComputeGradient(GradSpec[s],AuxField) == FAIL) {
	        fprintf(stderr, "grid::MultiSpeciesHandler: Error in SGSUtil_ComputeGradient(Spec)).\n");
	        return FAIL;
	      }
    }

    if (this->SGS_AddDiffusionTermsDE() == FAIL) {
	    fprintf(stderr, "grid::MultiSpeciesHandler: Error in SGS_AddDiffusionTermsDE.\n");
	    return FAIL;
    }
  }

  if (ProblemType == 62)
    this->CoolingTestResetEnergies();

  LCAPERF_STOP("grid_MultiSpeciesHandler");
  return SUCCESS;
}
