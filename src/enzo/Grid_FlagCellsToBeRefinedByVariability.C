// LI implemented in this Enzo version the modifications performed originally by WS, Jan. 2012

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

using namespace std;

int grid::FlagCellsToBeRefinedByVariability(int level)
{
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Error check. */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  int i, j, k, index, dim;
  float ave, stdev, thresh;

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* If no control variable is defined, return INT_UNDEFINED. */

  int NumberOfFlaggedCells = INT_UNDEFINED;

  /* Flag cells for which control variable is above threshold. */

  for (int method = 0; method < MAX_FLAGGING_METHODS; method++)

    // check if method has control variable
    if (ControlVariable[method] != NULL) {

      NumberOfFlaggedCells = 0;

      ave = ControlVariableAve[method];
      stdev = sqrt(ControlVariableVar[method]);

      // threshold is given by a mulitiple of the maximum of grid average and standard deviation
      thresh = fabs(ave) + ThreshFct[method]*((stdev > fabs(ave)) ? stdev : fabs(ave));
      thresh = max(thresh, ThreshMin[method]);

      if (debug) 
	cout << "[" << ProcessorNumber << "] Refinement by variability, method " 
	     << CellFlaggingMethod[method] << ", level " << level
	     << "\n[" << ProcessorNumber << "]     ave     = " << ave
	     << "\n[" << ProcessorNumber << "]     std dev = " << stdev
	     << "\n[" << ProcessorNumber << "]     thresh  = " << thresh << endl;

      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    index = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
	    FlaggingField[index] += (fabs(ControlVariable[method][index]) >= thresh) ? 1 : 0;
	  }
    }

  /* Count number of flagged cells. */

  if (NumberOfFlaggedCells != INT_UNDEFINED)
    for (i = 0; i < size; i++)
      NumberOfFlaggedCells += FlaggingField[i];

  return NumberOfFlaggedCells;
}
