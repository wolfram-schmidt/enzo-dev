/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE PRESENCE OF ENSTROPHY)
/
/  written by: Wolfram Schmidt
/  date:       February 2022
/
/  PURPOSE: refinement by vorticity using same normalization as 
/           refinement by shear; works for different solvers and
/           grid ranks
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

using namespace std;

int grid::FlagCellsToBeRefinedByRateOfStrainNormSqr()
{
  /* declarations */
 
  int i, j, k, index, dim;

  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* compute velocity derivative */

  if (this->ComputeJacobianVelocity(0) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeJacobianVelocity.");
  }

  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

 
  /* compute rate of strain squared */

  float *strain_sqr = new float[size];
	  
  if (this->ComputeJacobianVelocityNormSqr(strain_sqr) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeJacobianVelocityNormSqr.");
  }

  /* Combine all finite differences to determine if criterion is met */

  for (index = 0; index < size; index++) {
    FlaggingField[index] +=
      (strain_sqr[index] > (MinimumShearForRefinement*MinimumShearForRefinement)) ? 1 : 0;
  }
  
  /* clean up */
  delete [] strain_sqr;
 
  /* Count number of flagged Cells. */
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i])
      NumberOfFlaggedCells++;
 
  return NumberOfFlaggedCells;
}
