/***********************************************************************
/
/  GRID CLASS (COMPUTE CONTROL VARIABLES)
/
/  written by: Wolfram Schmidt
/  date:       April 2006
/  modified:   Feburary 2022
/
/  PURPOSE: computes control variables and sums over all grid cells
/           for refinement by variability
/
/  RETURNS: error status
/
************************************************************************/

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

int grid::ComputeControlVariables(int &active_zones, float* sum, float* sum_of_sqrs)
{
    if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

    int i, j, k, index, dim;

    int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);

    int size = 1;
    active_zones = 1;
    for (dim = 0; dim < GridRank; dim++) {
        size *= GridDimension[dim];
        active_zones *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    }

    //cout << "[" << MyProcessorNumber << "]" << " Computing control variables, size = " << size 
    //   << ", active zones = " << active_zones << endl;

    /* Compute velocity derivative */

    if (this->ComputeJacobianVelocity(0) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeJacobianVelocity.\n");
      return FAIL;
    }

    for (int method = 0; method < MAX_FLAGGING_METHODS; method++) {

      sum[method] = 0.0;
      sum_of_sqrs[method] = 0.0;

      if (CellFlaggingMethod[method] >= MIN_METHOD_VARIABILITY && CellFlaggingMethod[method] <= MAX_METHOD_VARIABILITY) {

        if (ControlVariable[method] == NULL)
            ControlVariable[method] = new float[size];

        for (i = 0; i < size; i++)
            ControlVariable[method][i] = 0.0;

        /* Select and control variable */

        switch (CellFlaggingMethod[method]) {

        /* ==== METHOD 20: BY DENSITY ==== */

        case 20: 
        
	    for (int i = 0; i < size; ++i)
            ControlVariable[method][i] = BaryonField[DensNum][i];
        break;

        /* ==== METHOD 21: BY ENSTROPHY ==== */

        case 21:

        if (this->ComputeVorticityNormSqr(ControlVariable[method]) == FAIL) {
            fprintf(stderr, "Error in grid->ComputeVorticityNormSqr.\n");
            return FAIL;
        }
        break;

        /* ==== METHOD 22: BY RATE OF STRAIN SQUARED ==== */

        case 22:

        if (this->ComputeRateOfStrainNormSqr(ControlVariable[method], 0) == FAIL) {
            fprintf(stderr, "Error in grid->ComputeRateOfStrainNormSqr.\n");
            return FAIL;
        }
        break;

        /* ==== METHOD 23: BY RATE OF COMPRESSION ==== */

        case 23:

        if (this->ComputeRateOfCompression(ControlVariable[method]) == FAIL) {
            fprintf(stderr, "Error in grid->ComputeRateOfCompression.\n");
            return FAIL;
        }
        break;

        default:
            ENZO_VFAIL("CellFlaggingMethod[%"ISYM"] = %"ISYM" not implemented for refinement by variability\n",
                       method, CellFlaggingMethod[method]);
        }

        for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
            for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
                for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
                    index = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
                    sum[method] += ControlVariable[method][index];
                    sum_of_sqrs[method] += ControlVariable[method][index]*ControlVariable[method][index];
                }

        //cout << "[" << MyProcessorNumber << "] method " << CellFlaggingMethod[method]
        //   << ", sum =  " << sum[method] << " sum of sqrs =  " << sum_of_sqrs[method] << endl;
      }
    }

    return SUCCESS;
}
