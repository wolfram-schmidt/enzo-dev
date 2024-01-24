/***********************************************************************
/
/  COMPUTE CONTROL VARIABLES
/
/  written by: Wolfram Schmidt
/  date:       Feburary 2022
/
/  PURPOSE: computes control variables for refinement by variability
/           and increments sums over all cells at given refinment level
/
/  RETURNS: error status
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <iostream>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
  
using namespace std;

int CommunicationAllSumValues(float *Values, int Number);
 
int ComputeControlVariables(HierarchyEntry *Grid, int level, int &active_zones, float* sum, float* sum_of_sqrs)
{
  grid *CurrentGrid = Grid->GridData;
 
  /* If this is the lowest allowed level, then return. */
 
  if (level >= MaximumRefinementLevel)
    return SUCCESS;
 
  /* If this grid is not on this processor, then return. */
 
  if (MyProcessorNumber != CurrentGrid->ReturnProcessorNumber())
    return SUCCESS;
 
  int grid_active_zones;
  float grid_sum[MAX_FLAGGING_METHODS];
  float grid_sum_of_sqrs[MAX_FLAGGING_METHODS];
 
  /* Compute control variables and grid sums. */

  if (CurrentGrid->ComputeControlVariables(grid_active_zones, grid_sum, grid_sum_of_sqrs) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeControlVariables.");
  }

  /* Increment sums for refinement level. */

  active_zones += grid_active_zones;
  for (int method = 0; method < MAX_FLAGGING_METHODS; method++)
    if (CellFlaggingMethod[method] >= MIN_METHOD_VARIABILITY && CellFlaggingMethod[method] <= MAX_METHOD_VARIABILITY) {                      
      //cout << "[" << MyProcessorNumber << "] method " << CellFlaggingMethod[method] << ", " 
      //	   << grid_active_zones << " " << grid_sum[method] << " " << grid_sum_of_sqrs[method] << endl;
      sum[method] += grid_sum[method];
      sum_of_sqrs[method] += grid_sum_of_sqrs[method];
    }

  return SUCCESS; 
}
