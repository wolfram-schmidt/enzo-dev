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
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int grid::FlagCellsToBeRefinedByVorticity()
{
  /* declarations */
 
  int i, j, k, index, dim, LeftOffset, RightOffset;
  float Vorticity, Denom, c_sound2;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate temp arrays and fill with zeroes */
 
  float *DelVel1 = new float[size];
  float *DelVel2 = new float[size];
  float *DelVel3 = new float[size];

  for (i = 0; i < size; i++) {
    DelVel1[i] = 0.0;
    DelVel2[i] = 0.0;
    DelVel3[i] = 0.0;
  }

  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  /* Create pressure field for calculating local sound speed */

  float *pressure = new float[size];

  if (this->ComputePressure(Time, pressure) == FAIL){
    ENZO_FAIL("Error in grid->ComputePressure.");
  }

  /* Set up stencil for Zeus vs PPM data, since PPM uses cell-centered 
     velocities and Zeus uses edge-centered velocities */
  if (HydroMethod == 2) { // Zeus
    LeftOffset = 1;
    RightOffset = 0;
    Denom = 1.0;
  }
  else {                 // PPM
    LeftOffset = 1;
    RightOffset = 1;
    Denom = 2.0;
  }

  /* loop over active dimensions */
 
  for (dim = 0; dim < GridRank; dim++) {
    if (GridDimension[dim] > 1) {
 
  /* Absolute Vorticity: [(dvy/dz - dvz/dy)^2 +
                          (dvz/dx - dvx/dz)^2 +
                          (dvx/dy - dvy/dx)^2]
     where:  dvx/dy = [vx(j-1) - vx(j+1)]/[2dy],
     in units of s^(-1)
    
     In order to use this as a refinement criterion, we need to compare
     it against the local cell width (ie dx).  So we refine when
     the shearing rate is larger than a user-defined threshold
     divided by the sound-crossing time of the cell:

     Vorticity > epsilon * c_s / dx

     We can save some extra calculations by removing dx from 
     the denominators of both sides of this equation:

     ((dvx - dvy)^2 + (dvz - dvx)^2 + (dvy - dvz)^2) / c_sound^2 > epsilon^2

     N.B.: DelVel1 = dvx(in dy) - dvy(in dx)
           DelVel2 = dvz(in dx) - dvx(in dz)
           DelVel3 = dvy(in dz) - dvz(in dy)
  */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    index = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];

	    switch (dim) {
	    case 0: 
	      if (GridRank > 1) // dvy (in dx)
		DelVel1[index] -= (BaryonField[Vel2Num][index - LeftOffset] -
				   BaryonField[Vel2Num][index + RightOffset]);
	      if (GridRank > 2) // dvz (in dx)
		DelVel2[index] += (BaryonField[Vel3Num][index - LeftOffset] -
				   BaryonField[Vel3Num][index + RightOffset]);
	      break;
	      
	    case 1:
	      if (GridRank > 1) // dvx (in dy)
		DelVel1[index] += (BaryonField[Vel1Num][index - LeftOffset] -
				   BaryonField[Vel1Num][index + RightOffset]);
	      if (GridRank > 2) // dvz (in dy)
		DelVel3[index] -= (BaryonField[Vel3Num][index - LeftOffset] -
				   BaryonField[Vel3Num][index + RightOffset]);
	      break;
	      
	    case 2:
	      if (GridRank > 2) { // dvx (in dz) and dvy (in dz)
		DelVel2[index] -= (BaryonField[Vel1Num][index - LeftOffset] -
				   BaryonField[Vel1Num][index + RightOffset]);
		DelVel3[index] += (BaryonField[Vel2Num][index - LeftOffset] -
				   BaryonField[Vel2Num][index + RightOffset]);
	      }
	      break;
	      
	    default:
	      break;
	    }
	  }
      LeftOffset *= GridDimension[dim];
      RightOffset *= GridDimension[dim];
    } 
  } // end loop over dimension

  /* Combine all finite differences to determine if criterion is met */
  
  for (index = 0; index < size; index++) {
    
    // If ideal gas, then calculate local sound speed
    if (EOSType == 0) 
      c_sound2 = Gamma*pressure[index]/BaryonField[DensNum][index];
    
    // If nonideal gas, then assume sound speed is 1
    else 
      c_sound2 = 1.0;
    
    Vorticity = (DelVel1[index]*DelVel1[index] + 
		 DelVel2[index]*DelVel2[index] + 
		 DelVel3[index]*DelVel3[index]) / 
      (Denom*Denom*c_sound2);
    
    FlaggingField[index] +=
      (Vorticity > (MinimumShearForRefinement*MinimumShearForRefinement)) ? 1 : 0;
  }
  
  /* clean up */
  delete [] DelVel1;
  delete [] DelVel2;
  delete [] DelVel3;
  delete [] pressure;
 
  /* Count number of flagged Cells. */
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i])

      NumberOfFlaggedCells++;
 
  return NumberOfFlaggedCells;
}
