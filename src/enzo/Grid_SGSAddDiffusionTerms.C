/***********************************************************************
 *
 * INFORMATION This file is part of a subgrid-scale (SGS) modeling 
 * framework in order to conduct explicit large eddy simulations (LES).
 *
 * The functions in this file concern models for the turbulent 
 * stress tensor in the momentum equation.
 * It consists of the turbulent (or SGS) Reynolds stress, the SGS
 * Maxwell stress, and the SGS magnetic pressure.
 *
 * The models have been verified "a priori", i.e. in comparison to
 * reference data, in 
 * Grete et al 2015 New J. Phys. 17 023070 doi: 10.1088/1367-2630/17/2/023070
 * Grete et al 2016 Phys. Plasmas 23 062317 doi: 10.1063/1.4954304 (Grete2016a)
 * and "a posteriori", i.e. used in simulations of decaying MHD turbulence, in
 * Grete et al 2017 Phys. Rev. E 05 033206 (Grete2017)
 *
 * WRITTEN BY Philipp Grete (mail@pgrete.de)
 *
 * DATE 2016
 *
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* 
 * This function adds to the SGS stress tensor (Reynolds stress component)
 * the pure (unscaled) nonlinear model
 * TauU = 1/12 * Delta^2 rho u_i,k u_j,k
 *
 * See equation (35) in Grete2016a for details (such as coefficient values)
 * or Vlaykov et al 2016 Phys. Plasmas 23 062316 doi: 10.1063/1.4954303 for
 * the derivation.
 */
void grid::SGS_AddDiff_nonlinear_energy(float **Flux) {
  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddDiff_nonlinear_energy start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  // if an explicit filter should be used
  // (at this point the fields are already filtered, 
  // see hydro_rk/Grid_MHDSourceTerms.C and the SGSNeedJacobians switch)
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  // if the model should be calculated based on grid-scale quantities
  // (not recommended, see Grete2017)
  } else {
    rho = BaryonField[DensNum];
  }

  int size = 1;
  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];

    /* we need Tau in the first ghost zone as well
     * as we'll take another derivative later on */
    StartIndex[dim] = GridStartIndex[dim] - 1;
    EndIndex[dim] = GridEndIndex[dim] + 1;
  }


  // the combined prefactor
  float CDeltaSqr = 1./12. * SGScoeffNLe * POW(SGSFilterWidth,2.) *
    POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        for (int l = 0; l < MAX_DIMENSION; l++) {
          Flux[SGSX][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSX][l][igrid] * GradEint[l][igrid];
          Flux[SGSY][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSY][l][igrid] * GradEint[l][igrid];
          Flux[SGSZ][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSZ][l][igrid] * GradEint[l][igrid];
        }
      }

}


/*
 * This function initializes a zero flux vector and calls the individual
 * functions that add the different terms to it.
 * Finally, the divergence of the flux is added to the dU vector used by
 * the MUSCL framework in hydro_rk/Grid_MHDSourceTerms.C 
 */
int grid::SGS_AddDiffusionTerms(float **dU) {
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  if (Time == 0.)
    return SUCCESS;

  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddDiffusionTerms start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  // if an explicit filter should be used
  // (at this point the fields are already filtered, 
  // see hydro_rk/Grid_MHDSourceTerms.C and the SGSNeedJacobians switch)
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  } else {
  // if the model should be calculated based on grid-scale quantities
  // (not recommended, see Grete2017)
    rho = BaryonField[DensNum];
  }

  int size = 1;
  float *Flux[6];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];
  }

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    Flux[dim] = new float[size];
    for (int i = 0; i < size; i++)
      Flux[dim][i] = 0.;
  }

  // the individual terms are added/activated by a non-zero coefficient
  if (SGScoeffNLe != 0.) 
    SGS_AddDiff_nonlinear_energy(Flux);

  int n = 0;
  int igrid, ip1, im1, jp1, jm1, kp1, km1;
  float EIncr;

  float facX = 1. / (2. * CellWidth[0][0]);
  float facY = 1. / (2. * CellWidth[1][0]);
  float facZ = 1. / (2. * CellWidth[2][0]);

  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];
        ip1 = i+1 + (j+k*GridDimension[1])*GridDimension[0];
        im1 = i-1 + (j+k*GridDimension[1])*GridDimension[0];
        jp1 = i + (j+1+k*GridDimension[1])*GridDimension[0];
        jm1 = i + (j-1+k*GridDimension[1])*GridDimension[0];
        kp1 = i + (j+(k+1)*GridDimension[1])*GridDimension[0];
        km1 = i + (j+(k-1)*GridDimension[1])*GridDimension[0];

        EIncr = - dtFixed * (
            (Flux[SGSX][ip1] - Flux[SGSX][im1])*facX + 
            (Flux[SGSY][jp1] - Flux[SGSY][jm1])*facY + 
            (Flux[SGSZ][kp1] - Flux[SGSZ][km1])*facZ);

        dU[iEtot][n] += EIncr;

	if (DualEnergyFormalism)
	  dU[iEint][n] += EIncr;
      }

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] Flux[dim];
  }

  return SUCCESS;
}
